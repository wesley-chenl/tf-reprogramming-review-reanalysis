"""Prepare the native shared Taiji genome index from real source BAM metadata.
This invokes only the official index node; it produces no TF predictions.
"""
from common import *
import subprocess,argparse
p=argparse.ArgumentParser();p.add_argument('--max-heap',default='14G');p.add_argument('--compact-gc',action='store_true');p.add_argument('--watch',action='store_true');a=p.parse_args()
work=OUT/'runs/preprocessing/TAIJI_GENOME_INDEX';work.mkdir(parents=True,exist_ok=True)
genome=OUT/'inputs/annotation/GRCm38.p6.genome.fa';index=OUT/'tmp/posix_work/TAIJI_mm10.genome_index';rows=[]
record=work/'RUN_RECORD.json'
if record.exists():
 old=json.loads(record.read_text())
 if old.get('exit_code')==0 and index.exists() and index.stat().st_size==old.get('native_index_bytes') and index.stat().st_mtime_ns==old.get('native_index_mtime_ns'):
  assert old['native_binary_sha256']==digest(OUT/'environments/bin/taiji') and old['genome_plain_bytes']==genome.stat().st_size and old['genome_plain_mtime_ns']==genome.stat().st_mtime_ns;print('Verified native index checkpoint retained',flush=True);exit()
while not all((OUT/f'inputs/atac/{sid}/filtered.bam.bai').exists() for sid in ['GSM1867737','GSM1867738']):
 assert a.watch,'Real source BAM prerequisites not ready'
 time.sleep(30)
for i,sid in enumerate(['GSM1867737','GSM1867738'],1):
 bam=OUT/f'inputs/atac/{sid}/filtered.bam';assert bam.exists() and Path(str(bam)+'.bai').exists()
 rows.append({'type':'ATAC-seq','id':'ATAC_source','group':'source','rep':i,'path':str(bam),'tags':'Filtered,PairedEnd','format':'Bam'})
tsv('runs/preprocessing/TAIJI_GENOME_INDEX/native_input.tsv',rows)
cfg={'input':str(work/'native_input.tsv'),'output_dir':str(OUT/'tmp/posix_work/TAIJI_RESOURCE_INDEX'),'assembly':'mm10','genome':str(genome),'genome_index':str(index),'annotation':str(OUT/'inputs/annotation/gencode.vM24.annotation.gtf'),'motif_file':str(OUT/'sources/Taiji_documentation/documentation-master/data/motifs/cisBP_mouse.meme'),'tmp_dir':str(OUT/'tmp/posix_work')}
dump('runs/preprocessing/TAIJI_GENOME_INDEX/native_config.yml',cfg)
env=os.environ.copy();env['LD_LIBRARY_PATH']=str(OUT/'environments/taiji_legacy_lib/lib/x86_64-linux-gnu');cmd=[str(OUT/'environments/bin/taiji'),'run','--config',str(work/'native_config.yml'),'--select','ATAC_Make_Index','-n','2','+RTS','-N2','-M'+a.max_heap,'-s']+(['-c'] if a.compact_gc else [])+['-RTS'];started=now()
dump('runs/preprocessing/TAIJI_GENOME_INDEX/RUN_RECORD.json',{'status':'RESOURCE_BUILD_RUNNING','command':cmd,'start_UTC':started,'exit_code':None,'no_TF_scores_generated':True})
with open(OUT/'logs/TAIJI_GENOME_INDEX.log','a') as log:
 print('\nINDEX_RESOURCE_ATTEMPT',started,'max_heap',a.max_heap,file=log,flush=True);r=subprocess.run(cmd,stdout=log,stderr=log,env=env,cwd=work)
rec={'command':cmd,'start_UTC':started,'end_UTC':now(),'exit_code':r.returncode,'genome':str(genome),'native_index':str(index),'native_binary_sha256':digest(OUT/'environments/bin/taiji'),'genome_gzip_sha256':json.loads((OUT/'inputs/bowtie2_mm10/BUILD_COMPLETE.json').read_text())['genome_sha256'],'genome_plain_bytes':genome.stat().st_size,'genome_plain_mtime_ns':genome.stat().st_mtime_ns,'no_TF_scores_generated':True,'no_reference_matching':True,'scope':'official genome index prerequisite only, not a benchmark task'}
if r.returncode==0:
 assert index.exists() and index.stat().st_size>0,'Native index node returned without a real index'
 rec.update(native_index_bytes=index.stat().st_size,native_index_mtime_ns=index.stat().st_mtime_ns,native_index_sha256=digest(index))
dump('runs/preprocessing/TAIJI_GENOME_INDEX/RUN_RECORD.json',rec)
print('Taiji native index prerequisite exit',r.returncode,flush=True);exit(r.returncode)
