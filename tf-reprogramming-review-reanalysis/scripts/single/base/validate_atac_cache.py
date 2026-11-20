"""Validate real completed ATAC cache contracts without rehashing large BAMs."""
from common import *
import argparse,subprocess,re
p=argparse.ArgumentParser();p.add_argument('--require-all',action='store_true');a=p.parse_args()
plans=list(csv.DictReader(open(OUT/'manifests/ATAC_DOWNLOAD_PLAN.tsv'),delimiter='\t'));ids=list(dict.fromkeys(r['sample_id'] for r in plans));records=[]
index_marker=OUT/'inputs/bowtie2_mm10/BUILD_COMPLETE.json'
if not index_marker.exists():
 assert not a.require_all,'Index is not complete'
 exit()
index_key=json.loads(index_marker.read_text())
for sid in ids:
 marker=OUT/f'inputs/atac/{sid}/PREPROCESSING_COMPLETE.json'
 if not marker.exists():
  assert not a.require_all,'ATAC pending: '+sid
  continue
 r=json.loads(marker.read_text());rr=[v for v in plans if v['sample_id']==sid];cmd=r['alignment_command'];assert r['all_raw_biological_reads_used'] is True
 assert '--very-sensitive-local' in cmd and '--mm' in cmd and cmd[cmd.index('--seed')+1]=='0' and cmd[cmd.index('-X')+1]=='2000'
 assert r['layout']==rr[0]['layout']
 if r['layout']=='PAIRED':assert '--no-mixed' in cmd and '--no-discordant' in cmd
 raw={}
 for v in rr:
  source=v['group'] in ('source','liver');stamp=OUT/('manifests/file_hashes' if source else 'manifests/atac_sra_cloud')/(v['filename']+'.json' if source else v['run']+'.json');s=json.loads(stamp.read_text());assert s['status']=='verified'
  raw[v['filename'] if source else v['run']]=s['sha256']
 assert raw==r['input_hashes'],'ATAC raw input identity changed: '+sid
 outputs={}
 for k in ['bam','peaks']:
  q=Path(r[k]);st=q.stat();assert st.st_mtime_ns<=marker.stat().st_mtime_ns
  outputs[str(q)]={'sha256':r[k+'_sha256'],'bytes':st.st_size,'mtime_ns':st.st_mtime_ns}
 bam=Path(r['bam']);assert Path(str(bam)+'.bai').exists();subprocess.run(['samtools','quickcheck','-v',str(bam)],check=True)
 header=subprocess.check_output(['samtools','view','-H',str(bam)],text=True);sizes={}
 for line in header.splitlines():
  if line.startswith('@SQ\t'):
   v=dict(k.split(':',1) for k in line.split('\t')[1:]);sizes[v['SN']]=int(v['LN'])
 assert sizes['chr1']==195471971 and sizes['chrX']==171031299,'Not expected mm10 coordinate dictionary'
 assert 'SO:coordinate' in header
 n=0
 with open(r['peaks']) as f:
  for line in f:
   v=line.rstrip().split('\t');assert len(v)==10 and v[0] in sizes
   assert 0<=int(v[1])<int(v[2])<=sizes[v[0]] and v[0] not in ('chrM','M')
   n+=1
 assert n>0,'No real peaks: '+sid
 flag=(bam.parent/'FLAGSTAT.txt').read_text()
 for label in ['secondary','supplementary','duplicates']:
  line=next(x for x in flag.splitlines() if x.endswith(' '+label));assert line.startswith('0 + 0 '),line
 streams=r['complete_SRA_stream_counts']
 for s in streams:
  reads=s['reads_by_mate'];assert reads[0]>0
  if r['layout']=='PAIRED':assert len(reads)==2 and reads[0]==reads[1]
 key={'raw_input_hashes':raw,'genome_index_build_key':index_key,'alignment':{k:(cmd[cmd.index(k)+1] if k in ('--seed','-X') else True) for k in ['--very-sensitive-local','--mm','--seed','-X']},'paired_rescue':'disabled for PE; only proper pairs kept','filter':'q30 F3852; PE f2; markdup-r; no chrM','peaks':'MACS3 q.01 keep-dup all; PE BAMPE; SE BAM nomodel shift-100 extsize200','layout':r['layout']}
 sha=hashlib.sha256(json.dumps(key,sort_keys=True).encode()).hexdigest()
 records.append({'sample_id':sid,'scientific_cache_key_sha256':sha,'scientific_cache_key':key,'completed_marker_sha256':digest(marker),'executed_code_sha256':r['uniform_pipeline_code_sha256'],'outputs':outputs,'n_peaks':n,'actual_alignment_threads':int(cmd[cmd.index('-p')+1]),'genome_coordinate_header_checked':True,'peak_coordinates_checked':True,'complete_SRA_stream_checked':bool(streams),'large_hash_source':'verified completion marker; output size/mtime unchanged; no redundant multi-GB rescan'})
if (OUT/'PRE_EVAL_FREEZE_C.json').exists():
 old=json.loads((OUT/'manifests/ATAC_CACHE_CONTRACTS.json').read_text());assert old['records']==records,'Frozen ATAC cache contract changed'
else:dump('manifests/ATAC_CACHE_CONTRACTS.json',{'UTC':now(),'n_completed':len(records),'expected':len(ids),'records':records})
print('Validated real ATAC cache/input/coordinate contracts:',len(records),'/',len(ids))
