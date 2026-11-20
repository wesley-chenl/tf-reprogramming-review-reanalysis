from method_io import *
import subprocess,argparse,os,shutil
p=argparse.ArgumentParser();p.add_argument('task');p.add_argument('--preflight',action='store_true');a=p.parse_args();task=a.task
if not a.preflight:assert (OUT/'PRE_EVAL_FREEZE_C.json').exists(),'Freeze actual ATAC/C inputs before C scoring'
index_record=OUT/'runs/preprocessing/TAIJI_GENOME_INDEX/RUN_RECORD.json';index_info=json.loads(index_record.read_text());assert index_info.get('exit_code')==0 and index_info.get('native_index_sha256'),'Complete and verify native index prerequisite first'
condition_marker=OUT/'inputs/atac_conditions/source/PEAK_POOL_COMPLETE.json';condition=json.loads(condition_marker.read_text());assert condition['exit_code']==0 and condition['source_biological_n']==2
access=guard_access(OUT/f'logs/access_{task}_C{"_preflight" if a.preflight else ""}.jsonl')
interface_policy_path=OUT/'config/C_NATIVE_INTERFACE_POLICY.json';interface_policy=json.loads(interface_policy_path.read_text());assert interface_policy['status']=='qualified_global_native_interface'
work=OUT/'runs'/task/'C';work.mkdir(parents=True,exist_ok=True)
native_tmp=OUT/'tmp/posix_work'/('TAIJI_'+task);native_tmp.mkdir(parents=True,exist_ok=True)
native=native_tmp/'output';native.mkdir(exist_ok=True)
def write_same_or_new(path,text):
 if not path.exists() or path.read_text()!=text:path.write_text(text)
started=now();tpm=OUT/f'inputs/rna/{task}/TPM.tsv';meta=pd.read_csv(OUT/f'inputs/rna/{task}/SAMPLE_METADATA.tsv',sep='\t');x=pd.read_csv(tpm,sep='\t',index_col=0);genes=pd.read_csv(OUT/'inputs/annotation/GENES_M24.tsv',sep='\t');x['gene_symbol']=genes.set_index('gene_id').gene_symbol.reindex(x.index)
# Taiji's official GeneQuant interface uses annotation gene_name. Pool exact
# duplicate names for its native expression table; ambiguous reverse TF IDs
# remain unresolved in the gene-level output, never assigned to multiple genes.
x=x.dropna(subset=['gene_symbol']).groupby('gene_symbol',sort=True).sum(numeric_only=True)
rows=[]
for role in ['source','target']:
 mm=meta[meta.group.eq('source' if role=='source' else task)]
 for i,sid in enumerate(mm.sample_id,1):
  path=work/f'RNA_{role}_{i}.tsv';write_same_or_new(path,x[[sid]].to_csv(sep='\t',header=False,index=True))
  rows.append({'type':'RNA-seq','id':'RNA_'+role,'group':role,'rep':i,'path':str(path),'tags':'GeneQuant','format':'Tsv'})
ast=pd.read_csv(OUT/'manifests/ATAC_SAMPLES_PREPROCESSING.tsv',sep='\t');ast=ast[ast.group.isin(['source',task])];input_paths=[tpm,OUT/f'inputs/rna/{task}/SAMPLE_METADATA.tsv',OUT/'inputs/annotation/GENES_M24.tsv',interface_policy_path];verified_outputs={}
shared_index=Path(index_info['native_index']);st=shared_index.stat();assert st.st_size==index_info['native_index_bytes'] and st.st_mtime_ns==index_info['native_index_mtime_ns'];verified_outputs[str(shared_index)]={'sha256':index_info['native_index_sha256'],'bytes':st.st_size,'mtime_ns':st.st_mtime_ns};input_paths +=[shared_index,index_record]
for r in ast.itertuples():
 marker=OUT/f'inputs/atac/{r.sample_id}/PREPROCESSING_COMPLETE.json';rec=json.loads(marker.read_text());role='source' if r.group=='source' else 'target';eid='ATAC_'+role
 # BAMs remain full-read preprocessing/QC lineage; the qualified native
 # GeneQuant+NarrowPeak interface does not consume them as native rows.
 input_paths +=[marker,Path(rec['bam'])]
 if role=='target':
  rows +=[{'type':'ATAC-seq','id':eid,'group':role,'rep':r.biological_replicate,'path':rec['peaks'],'tags':'','format':'NarrowPeak'}];input_paths +=[Path(rec['peaks'])]
 for typ in ['bam','peaks']:
  pp=Path(rec[typ]);assert pp.stat().st_mtime_ns<=marker.stat().st_mtime_ns;verified_outputs[str(pp)]={'sha256':rec[typ+'_sha256'],'bytes':pp.stat().st_size,'mtime_ns':pp.stat().st_mtime_ns}
pool_peak=Path(condition['peaks']);st=pool_peak.stat();assert st.st_size==condition['peaks_bytes'] and st.st_mtime_ns==condition['peaks_mtime_ns'];verified_outputs[str(pool_peak)]={'sha256':condition['peaks_sha256'],'bytes':st.st_size,'mtime_ns':st.st_mtime_ns};input_paths +=[condition_marker,pool_peak]
rows +=[{'type':'ATAC-seq','id':'ATAC_source','group':'source','rep':0,'path':str(pool_peak),'tags':'','format':'NarrowPeak'}]
write_same_or_new(work/'native_input.tsv',pd.DataFrame(rows).to_csv(sep='\t',index=False))
genome=OUT/'inputs/annotation/GRCm38.p6.genome.fa';annotation_stage=json.loads((OUT/'manifests/C_NATIVE_ANNOTATION_STAGE.json').read_text());assert annotation_stage['status']=='verified';gtf=Path(annotation_stage['native_logical_path']);input_paths +=[OUT/'manifests/C_NATIVE_ANNOTATION_STAGE.json'];motifs=OUT/'sources/Taiji_documentation/documentation-master/data/motifs/cisBP_mouse.meme'
# Cache keys contain actual RNA/ATAC/prior/binary hashes, not just filenames.
# Reuse already verified BAM hashes and immutable resource fingerprints rather
# than repeatedly scanning the same multi-GB objects.
cp=OUT/'manifests/C_RESOURCE_HASH_CACHE.json';hc=json.loads(cp.read_text()) if cp.exists() else {}
def fingerprint(pp):
 pp=Path(pp);st=pp.stat();key=str(pp);r=verified_outputs.get(key,hc.get(key))
 if r:
  if st.st_size!=r['bytes'] or st.st_mtime_ns!=r['mtime_ns']:
   assert not (OUT/'PRE_EVAL_FREEZE_C.json').exists(),'Frozen C cached resource changed: '+key
   h=hashlib.sha256()
   with open(pp,'rb') as f:
    while b:=f.read(1024*1024):h.update(b)
   r={'bytes':st.st_size,'mtime_ns':st.st_mtime_ns,'sha256':h.hexdigest()}
 else:
  h=hashlib.sha256()
  with open(pp,'rb') as f:
   while b:=f.read(1024*1024):h.update(b)
  r={'bytes':st.st_size,'mtime_ns':st.st_mtime_ns,'sha256':h.hexdigest()}
 hc[key]=r;return r['sha256']
native_dependencies={str(pp):fingerprint(pp) for pp in input_paths+[genome,gtf,motifs,OUT/'environments/bin/taiji',Path(__file__)]}
cp.write_text(json.dumps(hc,indent=2));cache_key=hashlib.sha256(json.dumps({'files':native_dependencies,'native':'Taiji1.3.0','inputs':'GeneQuant TPM; one condition NarrowPeak per group; BAMs preprocessing/QC lineage only','ranking':'max(PRt/PRs,PRs/PRt), no epsilon'},sort_keys=True).encode()).hexdigest()
native=native_tmp/('output_'+cache_key[:16]);native.mkdir(exist_ok=True);execution=work/('native_execution_'+cache_key[:16]);execution.mkdir(exist_ok=True)
dump(f'runs/{task}/C/NATIVE_INPUT_FINGERPRINT.json',{'cache_key':cache_key,'input_hashes':native_dependencies,'reference_or_evaluation_cache_used':False,'source':'same actual task inputs only'})
cfg={'input':str(work/'native_input.tsv'),'output_dir':str(native),'assembly':'mm10','genome':str(genome),'annotation':str(gtf),'motif_file':str(motifs),'genome_index':str(OUT/'tmp/posix_work/TAIJI_mm10.genome_index'),'tmp_dir':str(native_tmp),'callpeak_fdr':.01,'callpeak_genome_size':'mm'}
write_same_or_new(work/'native_config.yml',json.dumps(cfg,indent=2));env=os.environ.copy();env['LD_LIBRARY_PATH']=str(OUT/'environments/taiji_legacy_lib/lib/x86_64-linux-gnu')+(':'+env['LD_LIBRARY_PATH'] if env.get('LD_LIBRARY_PATH') else '')
# Official node selection includes the dependency closure for the native rank
# output. Unrelated BigWig, single-cell and presentation branches are excluded.
cmd=[str(OUT/'environments/bin/taiji'),'run','--config',str(work/'native_config.yml'),'--select','Output_Ranks','-n','9','+RTS','-N9','-M10G','-RTS']
ranks_existed_before_call=(native/'GeneRanks.tsv').exists()
with open(OUT/f'logs/Taiji_{task}{"_preflight" if a.preflight else ""}.log','w') as log:code=subprocess.run(cmd,stdout=log,stderr=log,env=env,cwd=execution).returncode
with open(OUT/f'logs/Taiji_{task}_durable_copy.log','a') as log:
 shutil.copytree(native, work/'native_Taiji', dirs_exist_ok=True)
ranks=native/'GeneRanks.tsv'
if code or not ranks.exists():
 if not a.preflight:complete_run(task,'C',task_status='FAILED_EXECUTION',reason=f'official_Taiji_native_pipeline_exit_{code}; GeneRanks_exists={ranks.exists()}',extra={'command':cmd,'start_UTC':started,'exit_code':code})
 else:dump(f'runs/preflight/{task}/C/CHECK.json',{'exit_code':code,'GeneRanks_exists':ranks.exists(),'command':cmd,'score_code_sha256':digest(Path(__file__)),'native_input_cache_key':cache_key,'reference_comparison_inspected':False})
 exit(1)
pr=pd.read_csv(ranks,sep='\t',index_col=0);pr.columns=[str(c) for c in pr.columns]
assert len(pr.columns)==2 and set(pr.columns)=={'source','target'},'Exactly one native network/PR column per condition is required: '+str(pr.columns.tolist())
lookup=genes.assign(native_key=genes.gene_symbol.astype(str).str.casefold()).groupby('native_key').gene_id.agg(list).to_dict();data=[]
for symbol,row in pr.iterrows():
 ids=lookup.get(str(symbol).casefold(),[])
 if len(ids)!=1:continue
 s=float(row['source']);t=float(row['target']);valid_PR=np.isfinite(s) and np.isfinite(t) and s>0 and t>0;score=max(t/s,s/t) if valid_PR else np.nan;finite=bool(valid_PR and np.isfinite(score))
 data.append({'gene_id':ids[0],'raw_score':score,'ranking_score':score,'PR_source':s,'PR_target':t,'signed_log_ratio':np.log(t)-np.log(s) if valid_PR else np.nan,'direction':'target_PR_greater' if finite and t>s else 'source_PR_greater' if finite and s>t else 'equal_PR' if finite else '', 'scorable':finite,'eligible_for_recommendation':finite,'candidate_status':'ranked' if finite else 'undefined_score','missing_reason':'' if finite else 'native_ratio_nonfinite' if valid_PR else 'native_PageRank_absent_zero_or_nonfinite'})
scores=pd.DataFrame(data,columns=['gene_id','raw_score','ranking_score','PR_source','PR_target','signed_log_ratio','direction','scorable','eligible_for_recommendation','candidate_status','missing_reason']).set_index('gene_id');scores.to_csv(work/'NATIVE_RATIO_COMPONENTS.tsv',sep='\t')
motif_symbols={line.split()[1].split('+',1)[0].casefold() for line in motifs.read_text().splitlines() if line.startswith('MOTIF ')}
u=pd.read_csv(OUT/'config/CANDIDATE_UNIVERSE.tsv',sep='\t');motif_coverage=[]
for r in u.itertuples():
 key=str(r.gene_symbol).casefold();mapped=lookup.get(key,[]);motif_available=key in motif_symbols
 motif_coverage.append({'gene_id':r.gene_id,'gene_symbol':r.gene_symbol,'native_cisBP_motif_available':motif_available,'n_case_insensitive_native_symbol_IDs':len(mapped),'native_network_rank_exported':r.gene_id in scores.index})
 if r.gene_id not in scores.index and motif_available:
  scores.loc[r.gene_id,['raw_score','ranking_score']]=np.nan;scores.loc[r.gene_id,['scorable','eligible_for_recommendation']]=False;scores.loc[r.gene_id,'candidate_status']='id_unresolved' if len(mapped)!=1 else 'input_absent';scores.loc[r.gene_id,'missing_reason']='native_case_insensitive_symbol_reverse_mapping_ambiguous_or_absent' if len(mapped)!=1 else 'no_native_network_rank_after_ATAC_motif_expression_linkage'
pd.DataFrame(motif_coverage).to_csv(work/'CISBP_MOUSE_MOTIF_COVERAGE.tsv',sep='\t',index=False)
if a.preflight:
 n_u=int(scores.scorable.reindex(u.gene_id).fillna(False).astype(bool).sum())
 dump(f'runs/preflight/{task}/C/CHECK.json',{'exit_code':code,'n_finite_ratio_all_native_genes':int(scores.scorable.sum()),'n_finite_ratio_in_frozen_U':n_u,'command':cmd,'score_code_sha256':digest(Path(__file__)),'native_input_cache_key':cache_key,'reference_comparison_inspected':False});print('C real preflight finite in U',n_u,'all native genes',scores.scorable.sum());exit()
input_paths +=[genome,gtf,motifs,OUT/'environments/bin/taiji',work/'native_input.tsv',work/'native_config.yml']
complete_run(task,'C',scores,inputs={str(p):digest(p) for p in input_paths},extra={'command':cmd,'start_UTC':started,'exit_code':code,'source_rna_consumed':True,'target_rna_consumed':True,'source_atac_consumed':True,'target_atac_consumed':True,'native_score_direction':'descending max(PR_target/PR_source,PR_source/PR_target)','signed_log_ratio_definition':'natural log(PR_target/PR_source)','external_network_used':False,'cocktail_enumeration_used':False,'native_input_interfaces':'GeneQuant TPM and condition NarrowPeak; full-read BAMs preprocessing/QC lineage only','native_BAM_input_rows':0,'native_interface_equivalence_evidence':str(OUT/'reports/C_NATIVE_PEAK_ONLY_EQUIVALENCE.json'),'source_ATAC_biological_replicates':2,'target_ATAC_biological_replicates':1,'native_input_code_prior_config_cache_key':cache_key,'native_outputs_existed_before_this_call':ranks_existed_before_call,'scope':'Taiji native RNA+ATAC network and personalized PageRank; ratio front-end only'})
