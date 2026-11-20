from common import *
import subprocess,collections,argparse
parser=argparse.ArgumentParser();parser.add_argument('--watch',action='store_true');args=parser.parse_args()
def linuxpath(s):
 if len(s)>2 and s[1]==':':return Path('/mnt/'+s[0].lower()+s[2:].replace('\\','/'))
 return Path(s)
plan=list(csv.DictReader(open(OUT/'manifests/RNA_DOWNLOAD_PLAN.tsv'),delimiter='\t'))
plan=[r for r in plan if r['group'] in ('source','liver')]
samples=collections.defaultdict(list)
for r in plan:samples[r['sample_id']].append(r)
index=OUT/'inputs/salmon_M24_k31';fasta=OUT/'inputs/annotation/gencode.vM24.transcripts.fa.gz'
salmon='salmon';version=subprocess.check_output([salmon,'--version'],text=True).strip()
ik={'fasta_sha256':digest(fasta),'version':version,'k':31,'gencode':True,'decoys':'none; transcriptome-only index'}
stamp=index/'BENCHMARK_INDEX.json'
if not stamp.exists() or json.loads(stamp.read_text())!=ik:
 with open(OUT/'logs/salmon_index.log','w') as log:
  subprocess.run([salmon,'index','-t',str(fasta),'-i',str(index),'--gencode','-k','31','-p','4'],stdout=log,stderr=log,check=True)
 stamp.write_text(json.dumps(ik,indent=2))
failed=set()
while True:
 pending=[]
 for sid,rr in samples.items():
  qdir=OUT/'inputs/salmon_quant'/sid;qdir.mkdir(parents=True,exist_ok=True);done=qdir/'BENCHMARK_QUANT.json'
  hs=[];ready=True
  for r in rr:
   cp=OUT/'manifests/file_hashes'/(r['filename']+'.json')
   if not cp.exists():ready=False;break
   h=json.loads(cp.read_text())
   if h.get('status')!='verified' or h.get('md5')!=r['expected_md5']:ready=False;break
   hs.append(h['sha256'])
  if not ready:pending.append(sid);continue
  key={'input_sha256':hs,'index':ik,'library_type':'A','seqBias':True,'gcBias':rr[0]['layout']=='PAIRED','validateMappings':True,'threads':4,'single_end_fldMean':200,'single_end_fldSD':80,'read_subsampling':False}
  if done.exists() and json.loads(done.read_text()).get('key')==key:continue
  if sid in failed:continue
  cmd=[salmon,'quant','-i',str(index),'-l','A','-p','4','--validateMappings','--seqBias','-o',str(qdir)]
  if rr[0]['layout']=='SINGLE':cmd+=['-r']+[str(linuxpath(r['local_path'])) for r in rr]+['--fldMean','200','--fldSD','80']
  else:
   mates=sorted(rr,key=lambda r:r.get('paired_end') or r['filename']);assert len(mates)==2
   cmd+=['-1',str(linuxpath(mates[0]['local_path'])),'-2',str(linuxpath(mates[1]['local_path'])),'--gcBias']
  rec={'sample_id':sid,'command':cmd,'key':key,'started_UTC':now()}
  progress('Real Salmon quantification: '+sid)
  with open(OUT/'logs'/('salmon_'+sid+'.log'),'w') as log:code=subprocess.run(cmd,stdout=log,stderr=log).returncode
  rec.update(exit_code=code,ended_UTC=now())
  if code==0 and (qdir/'quant.sf').exists():
   rec.update(status='completed',quant_sha256=digest(qdir/'quant.sf'),meta=json.loads((qdir/'aux_info/meta_info.json').read_text()));done.write_text(json.dumps(rec,indent=2))
   print(sid,'quantified mapping_rate',rec['meta'].get('percent_mapped'),flush=True)
  else:failed.add(sid);dump('runs/preprocessing/'+sid+'_FAILED.json',rec)
 if not args.watch or not pending:break
 print('Waiting for validated FASTQ:',','.join(pending),flush=True);time.sleep(30)
progress('Salmon quantification monitor stopped; inspect per-sample stamps for complete and failed samples.')
