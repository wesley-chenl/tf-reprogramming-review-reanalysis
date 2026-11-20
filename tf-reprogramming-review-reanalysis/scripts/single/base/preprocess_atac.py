"""Uniform, complete-read mm10 reconstruction; no literature labels used."""
from common import *
from stream_sra import stage_and_validate,Streamer
import subprocess,gzip,shutil,collections,concurrent.futures,argparse,fcntl
p=argparse.ArgumentParser();p.add_argument('--watch',action='store_true');a=p.parse_args()
SCRIPT_SHA_AT_IMPORT=digest(Path(__file__))
plans=list(csv.DictReader(open(OUT/'manifests/ATAC_DOWNLOAD_PLAN.tsv'),delimiter='\t'));samples=collections.defaultdict(list)
for r in plans:samples[r['sample_id']].append(r)
genome=OUT/'inputs/annotation/GRCm38.p6.genome.fa';gzg=Path(str(genome)+'.gz');gtf=OUT/'inputs/annotation/gencode.vM24.annotation.gtf';index=OUT/'inputs/bowtie2_mm10/genome';index.parent.mkdir(exist_ok=True)
def run(cmd,log,**kwargs):
 print('COMMAND',cmd,file=log,flush=True);return subprocess.run(cmd,stdout=log,stderr=log,**kwargs).returncode
def durable_copy(src,dst,log):
 # Windows native output avoids slow DrvFS 9P writes to the SATA HDD.
 shutil.copyfile(src,dst)

while not gzg.exists():
 if not a.watch:raise SystemExit('Genome download not yet complete')
 time.sleep(30)
for gz,plain in [(gzg,genome),(Path(str(gtf)+'.gz'),gtf)]:
 if not plain.exists():
  partial=Path(str(plain)+'.partial')
  with gzip.open(gz,'rb') as g,open(partial,'wb') as f:shutil.copyfileobj(g,f,8*1024*1024)
  partial.replace(plain)
key={'genome_sha256':digest(gzg),'bowtie2_build_version':subprocess.check_output(['bowtie2-build','--version'],text=True).splitlines()[0],'options':['--threads','4','--bmaxdivn','8','--dcv','1024'],'index_build_storage':'native SSD temporary prefix, then durable HDD copy; same genome and alignment algorithm'}
stamp=index.parent/'BUILD_COMPLETE.json'
if not stamp.exists() or json.loads(stamp.read_text())!=key:
 native_index=OUT/'tmp/posix_work/BOWTIE2_INDEX/genome';native_index.parent.mkdir(parents=True,exist_ok=True);native_fasta=native_index.parent/'GRCm38.p6.genome.fa'
 if not native_fasta.exists() or native_fasta.stat().st_size!=genome.stat().st_size:shutil.copyfile(genome,native_fasta)
 with open(OUT/'logs/bowtie2_index.log','w') as log:
  products=[p for p in native_index.parent.glob('genome.*') if p.suffix in ('.bt2','.bt2l')]
  if len(products)!=6:
   code=run(['bowtie2-build','--threads','4','--bmaxdivn','8','--dcv','1024',str(native_fasta),str(native_index)],log)
   if code:raise RuntimeError('Bowtie2 index build failed')
  assert run(['bowtie2-inspect','-s',str(native_index)],log)==0
  products=[p for p in native_index.parent.glob('genome.*') if p.suffix in ('.bt2','.bt2l')];assert len(products)==6
  for q in products:durable_copy(q,index.parent/q.name,log)
 stamp.write_text(json.dumps(key,indent=2))
 for q in products:q.unlink()
 native_fasta.unlink()
# Read-only native index shared by memory mapping across concurrent aligners.
# Byte-identical to the retained durable index; computational storage adapter.
align_index=OUT/'tmp/posix_work/BOWTIE2_SHARED_INDEX/genome';align_index.parent.mkdir(exist_ok=True)
shared_stamp=align_index.parent/'BYTE_IDENTITY.json'
if not shared_stamp.exists():
 copied=[]
 for src in index.parent.glob('genome.*'):
  if src.suffix not in ('.bt2','.bt2l'):continue
  dst=align_index.parent/src.name;h=hashlib.sha256()
  with open(src,'rb') as f,open(dst,'wb') as g:
   while b:=f.read(8*1024*1024):h.update(b);g.write(b)
  assert dst.stat().st_size==src.stat().st_size;copied.append({'source':str(src),'native_copy':str(dst),'bytes':dst.stat().st_size,'sha256_of_copy_stream':h.hexdigest()})
 assert len(copied)==6
 shared_stamp.write_text(json.dumps({'genome_index_key':key,'files':copied,'storage_only_adapter':True},indent=2))
toolkit=next((OUT/'environments/sratoolkit').glob('sratoolkit.*/bin'))
subprocess.run([sys.executable,str(OUT/'code/validate_atac_cache.py')],check=True)
def lpath(s):return Path('/mnt/'+s[0].lower()+s[2:].replace('\\','/')) if len(s)>2 and s[1]==':' else Path(s)
failed=set()
def ready(sid,rr):
 if rr[0]['group']=='source':return all(lpath(r['local_path']).exists() for r in rr)
 if rr[0]['group']=='liver':return all((OUT/'manifests/file_hashes'/(r['filename']+'.json')).exists() and json.loads((OUT/'manifests/file_hashes'/(r['filename']+'.json')).read_text()).get('status')=='verified' for r in rr)
 return all((OUT/'manifests/atac_sra_cloud'/(r['run']+'.json')).exists() and json.loads((OUT/'manifests/atac_sra_cloud'/(r['run']+'.json')).read_text()).get('status')=='verified' for r in rr)
def one(sid):
 rr=samples[sid];dest=OUT/'inputs/atac'/sid;dest.mkdir(parents=True,exist_ok=True);done=dest/'PREPROCESSING_COMPLETE.json'
 if done.exists():return True
 start=now();f1=[];f2=[];derived=[];archives=[];stream=None;hashes={};layout=rr[0]['layout'];logp=OUT/'logs'/('ATAC_'+sid+'.log')
 with open(logp,'a') as log:
  if rr[0]['group'] in ('source','liver'):
   for r in rr:
    fp=lpath(r['local_path']);cp=OUT/'manifests/file_hashes'/(r['filename']+'.json')
    if cp.exists() and json.loads(cp.read_text()).get('status')=='verified':rec=json.loads(cp.read_text())
    else:
     md=hashlib.md5();sh=hashlib.sha256()
     with open(fp,'rb') as f:
      while b:=f.read(8*1024*1024):md.update(b);sh.update(b)
     assert fp.stat().st_size==int(r['expected_bytes']) and md.hexdigest()==r['expected_md5']
     rec={'status':'verified','sha256':sh.hexdigest(),'md5':md.hexdigest(),'bytes':fp.stat().st_size,'path':str(fp)};dump('manifests/file_hashes/'+r['filename']+'.json',rec)
    hashes[r['filename']]=rec['sha256']
   if layout=='PAIRED':
    if rr[0]['group']=='liver':f1=[lpath(r['local_path']) for r in rr if r.get('paired_end')=='1'];f2=[lpath(r['local_path']) for r in rr if r.get('paired_end')=='2']
    else:f1=[lpath(r['local_path']) for r in rr if '_1.fastq' in r['filename']];f2=[lpath(r['local_path']) for r in rr if '_2.fastq' in r['filename']]
   else:f1=[lpath(r['local_path']) for r in rr]
  else:
   for runid in dict.fromkeys(r['run'] for r in rr):
    master_sr=OUT/'inputs/raw_atac_sra'/(runid+'.sra');srrec=json.loads((OUT/'manifests/atac_sra_cloud'/(runid+'.json')).read_text());hashes[runid]=srrec['sha256'];work=OUT/'tmp/posix_work'/('ATAC_'+runid)
    archives.append(stage_and_validate(master_sr,srrec['sha256'],work,log))
   stream=Streamer(archives,layout,OUT/'tmp/posix_work'/('STREAM_ATAC_'+sid),log);f1=[stream.paths[0]];f2=[stream.paths[1]] if layout=='PAIRED' else []
  if layout=='PAIRED':assert len(f1)==len(f2)>0
  cmd=['bowtie2','--very-sensitive-local','--mm','--seed','0','-p','12','-X','2000','-x',str(align_index)]
  if layout=='PAIRED':cmd+=['--no-mixed','--no-discordant']
  cmd+=['-1',','.join(map(str,f1)),'-2',','.join(map(str,f2))] if layout=='PAIRED' else ['-U',','.join(map(str,f1))]
  align_command=cmd.copy()
  native_temp=OUT/'tmp/posix_work'/('ALIGN_'+sid);native_temp.mkdir(parents=True,exist_ok=True)
  names=native_temp/'names.tmp.bam';coord=native_temp/'deduplicated.tmp.bam';filtered=dest/'filtered.bam';native_filtered=native_temp/'filtered.bam'
  if not filtered.exists():
   print('ALIGN_COMMAND',cmd,file=log,flush=True)
   align=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log)
   sort=subprocess.Popen(['samtools','sort','-n','-@','1','-m','256M','-o',str(names),'-'],stdin=align.stdout,stdout=log,stderr=log);align.stdout.close()
   if stream:stream.start()
   assert sort.wait()==0 and align.wait()==0
   stream_counts=stream.join() if stream else []
   fix=subprocess.Popen(['samtools','fixmate','-m',str(names),'-'],stdout=subprocess.PIPE,stderr=log)
   sort=subprocess.Popen(['samtools','sort','-@','1','-m','512M','-T',str(native_temp/'coordinate_sort'),'-'],stdin=fix.stdout,stdout=subprocess.PIPE,stderr=log);fix.stdout.close()
   mark=subprocess.Popen(['samtools','markdup','-r','-s','-',str(coord)],stdin=sort.stdout,stdout=log,stderr=log);sort.stdout.close();assert mark.wait()==0 and sort.wait()==0 and fix.wait()==0
   assert run(['samtools','index',str(coord)],log)==0
   chrs=[line.split('\t')[0] for line in subprocess.check_output(['samtools','idxstats',str(coord)],text=True).splitlines() if line.split('\t')[0] not in ('chrM','M','*')]
   cmd=['samtools','view','-b','-q','30','-F','3852']+(['-f','2'] if layout=='PAIRED' else [])+['-o',str(native_filtered),str(coord)]+chrs

   assert run(['samtools','index',str(native_filtered)],log)==0
   durable_copy(native_filtered,filtered,log);durable_copy(Path(str(native_filtered)+'.bai'),Path(str(filtered)+'.bai'),log)
   native_filtered.unlink();Path(str(native_filtered)+'.bai').unlink()
   for tmpfile in [names,coord,Path(str(coord)+'.bai')]:tmpfile.unlink()
  assert run(['samtools','quickcheck','-v',str(filtered)],log)==0
  if stream and not hasattr(stream,'thread'):
   for fd in stream.keep:os.close(fd)
   for pipe in stream.paths:pipe.unlink()
  stream_counts=stream.counts if stream else []
  with open(dest/'FLAGSTAT.txt','w') as f:subprocess.run(['samtools','flagstat',str(filtered)],stdout=f,stderr=log,check=True)
  peaks=dest/'peaks_peaks.narrowPeak'
  peak_cmd=['macs3','callpeak','-t',str(filtered),'-f','BAMPE' if layout=='PAIRED' else 'BAM','-g','mm','-q','.01','--keep-dup','all','-n','peaks','--outdir',str(dest)]
  if layout!='PAIRED':peak_cmd+=['--nomodel','--shift','-100','--extsize','200']
  if not peaks.exists():
   assert run(peak_cmd,log)==0
  rec={'sample_id':sid,'group':rr[0]['group'],'start_UTC':start,'end_UTC':now(),'input_hashes':hashes,'alignment_command':align_command,'peak_command':peak_cmd,'layout':layout,'bam':str(filtered),'bam_sha256':digest(filtered),'peaks':str(peaks),'peaks_sha256':digest(peaks),'status':'completed','biological_replicate':rr[0]['biological_replicate'],'all_raw_biological_reads_used':True,'complete_SRA_stream_counts':stream_counts,'technical_SRA_runs_merged_before_dedup':len(archives),'n_input_files':len(rr),'blacklist':'not applied; native Taiji default; chrM removed before peak calling','uniform_pipeline_code_sha256':SCRIPT_SHA_AT_IMPORT,'alignment_qualification':'uniform local soft clipping, no mixed/discordant rescue retained for PE; same-input shared index, seed0; QC-only decision before any C reference scoring'}
  done.write_text(json.dumps(rec,indent=2));print(sid,'ATAC real preprocessing complete',flush=True)
  for p in derived:p.unlink()
  for archive in archives:archive.unlink()
 return True
with concurrent.futures.ThreadPoolExecutor(2) as ex:
 futures={}
 while True:
  running=set(futures.values());pending=[sid for sid in samples if not (OUT/'inputs/atac'/sid/'PREPROCESSING_COMPLETE.json').exists() and sid not in failed and sid not in running]
  available=[sid for sid in pending if ready(sid,samples[sid])]
  for sid in available[:max(0,2-len(futures))]:futures[ex.submit(one,sid)]=sid
  if futures:
   done,_=concurrent.futures.wait(futures,timeout=10,return_when=concurrent.futures.FIRST_COMPLETED)
   for f in done:
    sid=futures.pop(f)
    try:f.result();progress('Real complete-read ATAC preprocessing completed: '+sid)
    except Exception as e:failed.add(sid);dump('runs/preprocessing/ATAC_'+sid+'_FAILED.json',{'sample_id':sid,'error':repr(e),'UTC':now()});print(sid,'failed',repr(e),flush=True)
  elif not pending or not a.watch:break
  else:time.sleep(30)
