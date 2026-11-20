"""Stream complete SRA biological reads into native program paired FIFOs."""
from pathlib import Path
import os,threading,subprocess,hashlib,shutil,json,time,fcntl,queue
from common import OUT,now
BIN=next((OUT/'environments/sratoolkit').glob('sratoolkit.*/bin'))
def stage_and_validate(master,sha256,work,log):
 work=Path(work);work.mkdir(parents=True,exist_ok=True);sra=work/Path(master).name;stamp=work/'ARCHIVE_VALIDATED.json'
 cache=OUT/'tmp/posix_work/DOWNLOADS'/Path(master).name
 if cache.exists() and cache.stat().st_size==Path(master).stat().st_size:sra=cache
 if stamp.exists() and sra.exists() and json.loads(stamp.read_text()).get('sha256')==sha256:return sra
 if not sra.exists() or sra.stat().st_size!=Path(master).stat().st_size:shutil.copyfile(master,sra)
 # Reuse byte checksum from the verified master; native archive CRC/MD5 validation
 # detects local stage corruption and is retained alongside the command log.
 cmd=[str(BIN/'vdb-validate'),str(sra)];print('VALIDATE_COMMAND',cmd,file=log,flush=True)
 code=subprocess.run(cmd,stdout=log,stderr=log).returncode
 if code:raise RuntimeError('vdb-validate failed '+str(master))
 stamp.write_text(json.dumps({'master':str(master),'sha256':sha256,'exit_code':code,'UTC':now()}))
 return sra
class Streamer:
 def __init__(self,archives,layout,work,log):
  self.archives=archives;self.layout=layout;self.work=Path(work);self.work.mkdir(parents=True,exist_ok=True);self.log=log;self.counts=[];self.error=None;self.process=None
  self.paths=[self.work/'read1.fifo']+([self.work/'read2.fifo'] if layout=='PAIRED' else [])
  for p in self.paths:
   if p.exists():p.unlink()
   os.mkfifo(p)
  # RDWR keeps each pipe open while the consuming program initializes both mates.
  self.keep=[os.open(p,os.O_RDWR) for p in self.paths]
  for fd in self.keep:
   try:fcntl.fcntl(fd,fcntl.F_SETPIPE_SZ,1024*1024)
   except OSError:pass
 def start(self):self.thread=threading.Thread(target=self._feed,daemon=True);self.thread.start()
 def _feed(self):
  try:
   fs=[os.fdopen(os.dup(fd),'wb',buffering=0) for fd in self.keep]
   # Independent mate writers are required: native FASTQ readers may prefetch
   # a complete block from mate1 before opening/draining mate2. Alternating
   # blocking writes into two small FIFOs can deadlock on perfectly valid PE
   # input. Bounded queues keep at most 128 MiB/mate and preserve every byte.
   queues=[queue.Queue(maxsize=128) for _ in fs]
   writer_errors=[]
   def writer(j):
    try:
     while True:
      b=queues[j].get()
      if b is None:break
      view=memoryview(b)
      while view:
       n=fs[j].write(view);view=view[n:]
    except Exception as e:writer_errors.append(repr(e))
    finally:
     fs[j].close()
     # Release this mate's RDWR startup keeper immediately. Waiting until both
     # writers finish can withhold EOF on mate1 while a native decoder is
     # waiting to finish mate1 before draining mate2's final buffered block.
     fd=self.keep[j];self.keep[j]=None
     if fd is not None:os.close(fd)
   writers=[threading.Thread(target=writer,args=(j,),daemon=True) for j in range(len(fs))]
   for t in writers:t.start()
   buffers=[bytearray() for _ in fs]
   for archive in self.archives:
    cmd=[str(BIN/'fastq-dump'),'--split-spot','--skip-technical','--readids','--stdout',str(archive)]
    print('FULL_READ_STREAM_COMMAND',cmd,file=self.log,flush=True);proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=self.log,bufsize=4*1024*1024,start_new_session=True);self.process=proc;count=[0]*len(fs)
    while True:
     h=proc.stdout.readline()
     if not h:break
     record=[h,proc.stdout.readline(),proc.stdout.readline(),proc.stdout.readline()]
     if not record[3] or not record[2].startswith(b'+'):raise RuntimeError('Invalid streamed FASTQ record')
     if self.layout=='PAIRED':
      key=h.split()[0].rsplit(b'.',1)[-1]
      if key not in (b'1',b'2'):raise RuntimeError('Unresolved biological read index '+repr(h))
      mate=int(key)-1
      record[0]=h.split()[0].rsplit(b'.',1)[0]+b'/'+key+b'\n'
     else:mate=0
     buffers[mate].extend(b''.join(record));count[mate]+=1
     if len(buffers[mate])>=1024*1024:
      queues[mate].put(bytes(buffers[mate]));buffers[mate].clear()
    proc.stdout.close();code=proc.wait()
    if code:raise RuntimeError('fastq-dump exit '+str(code))
    if len(count)==2 and count[0]!=count[1]:raise RuntimeError('Unbalanced full biological read mates')
    self.counts.append({'archive':str(archive),'reads_by_mate':count,'exit_code':code})
   for j,b in enumerate(buffers):
    if b:queues[j].put(bytes(b))
    queues[j].put(None)
   for t in writers:t.join()
   if writer_errors:raise RuntimeError(';'.join(writer_errors))
  except Exception as e:self.error=repr(e)
  finally:
   # Close only after every byte is flushed so the consumer observes EOF.
   for fd in self.keep:
    if fd is not None:
     try:os.close(fd)
     except OSError:pass
 def join(self):
  self.thread.join()
  for p in self.paths:
   if p.exists():p.unlink()
  if self.error:raise RuntimeError(self.error)
  return self.counts
