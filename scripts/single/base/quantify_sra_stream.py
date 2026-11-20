from common import *
from stream_sra import stage_and_validate,Streamer
import subprocess,collections,argparse,signal,shutil
p=argparse.ArgumentParser();p.add_argument('--watch',action='store_true');a=p.parse_args()
plan=list(csv.DictReader(open(OUT/'manifests/RNA_DOWNLOAD_PLAN.tsv'),delimiter='\t'));rr={r['run']:r for r in plan if r['group'] not in ('source','liver')};index=OUT/'inputs/salmon_M24_k31';failed=set()
# Native cached index is byte-identical to the common frozen transcriptome
# index; placing it on the POSIX volume avoids memory-mapping DrvFS failures.
native_index=OUT/'tmp/posix_work/SALMON_INDEX';native_index.mkdir(exist_ok=True)
for q in index.iterdir():
 if not q.is_file():continue
 dst=native_index/q.name
 if not dst.exists() or dst.stat().st_size!=q.stat().st_size:shutil.copyfile(q,dst)
assert (native_index/'BENCHMARK_INDEX.json').read_bytes()==(index/'BENCHMARK_INDEX.json').read_bytes()
while True:
 pending=[]
 for run,r in rr.items():
  sid=r['sample_id'];stamp=OUT/'manifests/sra_cloud'/(run+'.json');qdir=OUT/'inputs/salmon_quant'/sid;done=qdir/'BENCHMARK_QUANT.json'
  if done.exists() and json.loads(done.read_text()).get('key',{}).get('SRA_cloud_run')==run:continue
  if not stamp.exists():pending.append(run);continue
  rec=json.loads(stamp.read_text())
  if rec.get('status')!='verified':pending.append(run);continue
  if run in failed:continue
  qdir.mkdir(parents=True,exist_ok=True);work=OUT/'tmp/posix_work'/('RNA_'+run);work.mkdir(parents=True,exist_ok=True);started=now()
  try:
   with open(OUT/'logs'/(run+'_SRA_STREAM_QUANT.log'),'a') as log:
    sra=stage_and_validate(OUT/'inputs/raw_sra'/(run+'.sra'),rec['sha256'],work,log)
    stream=Streamer([sra],r['layout'],work/'streams',log)
    cmd=['salmon','quant','-i',str(native_index),'-l','A','-p','4','--validateMappings','--seqBias','-o',str(qdir)]
    if r['layout']=='PAIRED':cmd+=['--gcBias','-1',str(stream.paths[0]),'-2',str(stream.paths[1])]
    else:cmd+=['-r',str(stream.paths[0]),'--fldMean','200','--fldSD','80']
    print('SALMON_COMMAND',cmd,file=log,flush=True);proc=subprocess.Popen(cmd,stdout=log,stderr=log);stream.start();code=proc.wait()
    if code:
     if stream.process:
      try:os.killpg(stream.process.pid,signal.SIGTERM)
      except ProcessLookupError:pass
     raise RuntimeError('Salmon failed before stream completion; inspect resource/index log, exit '+str(code))
    counts=stream.join()
    assert code==0
   meta=json.loads((qdir/'aux_info/meta_info.json').read_text());assert sum(c['reads_by_mate'][0] for c in counts)==meta['num_processed'];status={'sample_id':sid,'run':run,'key':{'SRA_cloud_run':run,'SRA_sha256':rec['sha256'],'index':json.loads((index/'BENCHMARK_INDEX.json').read_text()),'library_type':'A','seqBias':True,'gcBias':r['layout']=='PAIRED','read_subsampling':False,'SRA_stream_tool':'fastq-dump3.4.1','paired_header_normalization':'.readID to /mate'},'started_UTC':started,'ended_UTC':now(),'command':cmd,'exit_code':code,'status':'completed','complete_biological_read_stream':counts,'meta':meta,'quant_sha256':digest(qdir/'quant.sf')}
   done.write_text(json.dumps(status,indent=2));dump('manifests/sra_cloud/'+run+'_QUANT_PROVENANCE.json',status);sra.unlink();print(sid,run,'real all-read quantification complete mapping',meta.get('percent_mapped'),'stream',counts,flush=True)
   progress('Real complete-read RNA quantification completed: '+sid)
  except Exception as e:
   failed.add(run);dump('runs/preprocessing/'+run+'_STREAM_FAILED.json',{'run':run,'error':repr(e),'UTC':now()});print(run,'failed',repr(e),flush=True)
 if not a.watch or not pending:break
 print('Waiting for SRA archives:',','.join(pending),flush=True);time.sleep(30)
