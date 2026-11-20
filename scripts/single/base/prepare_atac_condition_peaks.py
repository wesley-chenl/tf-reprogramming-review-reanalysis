"""One source-condition peak object from both complete biological libraries.
Native Taiji calls peaks after merging replicate beds; provided per-replicate
NarrowPeak objects bypass that merge and must not duplicate a network group.
"""
from common import *
import subprocess,argparse
p=argparse.ArgumentParser();p.add_argument('--watch',action='store_true');a=p.parse_args()
sids=['GSM1867737','GSM1867738'];dest=OUT/'inputs/atac_conditions/source';dest.mkdir(parents=True,exist_ok=True);marker=dest/'PEAK_POOL_COMPLETE.json';peaks=dest/'source_pool_peaks.narrowPeak'
while not all((OUT/f'inputs/atac/{s}/PREPROCESSING_COMPLETE.json').exists() for s in sids):
 assert a.watch,'Completed source biological BAM prerequisites required'
 time.sleep(30)
records=[json.loads((OUT/f'inputs/atac/{s}/PREPROCESSING_COMPLETE.json').read_text()) for s in sids]
for r in records:assert r['layout']=='PAIRED' and r['all_raw_biological_reads_used'] is True
cmd=['macs3','callpeak','-t']+[r['bam'] for r in records]+['-f','BAMPE','-g','mm','-q','.01','--keep-dup','all','-n','source_pool','--outdir',str(dest)]
key={'source_library_BAM_hashes':{r['sample_id']:r['bam_sha256'] for r in records},'parameters':['BAMPE','mm','q.01','keep-dup all'],'MACS3_version':subprocess.check_output(['macs3','--version'],text=True).strip(),'code_sha256':digest(Path(__file__))}
if marker.exists():
 old=json.loads(marker.read_text())
 assert old['scientific_cache_key']==key and peaks.stat().st_size==old['peaks_bytes'] and peaks.stat().st_mtime_ns==old['peaks_mtime_ns'];print('Verified source condition peak cache retained',flush=True);exit()
started=now()
with open(OUT/'logs/ATAC_SOURCE_CONDITION_POOL.log','w') as log:
 print('COMMAND',cmd,file=log,flush=True);r=subprocess.run(cmd,stdout=log,stderr=log)
assert r.returncode==0 and peaks.exists(),'Condition peak call failed; inspect real log'
n=sum(1 for _ in open(peaks));assert n>0
dump('inputs/atac_conditions/source/PEAK_POOL_COMPLETE.json',{'scientific_cache_key':key,'source_biological_samples':sids,'source_biological_n':2,'native_condition_group':'source','native_peak_replicate_tag':0,'individual_library_QC_peaks_preserved':True,'BAMs_already_deduplicated_within_biological_library':True,'duplicates_between_independent_libraries_not_removed':True,'command':cmd,'start_UTC':started,'end_UTC':now(),'exit_code':r.returncode,'peaks':str(peaks),'peaks_sha256':digest(peaks),'peaks_bytes':peaks.stat().st_size,'peaks_mtime_ns':peaks.stat().st_mtime_ns,'n_peaks':n,'reference_matching_used':False,'reason':'native workflow pools replicate reads before one condition peak call; prevents repeated group network path overwrite from supplied library peak objects'})
progress('Native C input contract qualified: both complete source ATAC biological libraries pooled only for one condition-level peak call; individual BAM/QC products and biological n2 retained. No reference matching.');print('Source condition peak call complete, peaks',n,flush=True)
