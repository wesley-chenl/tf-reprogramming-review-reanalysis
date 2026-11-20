from method_io import *
import subprocess,argparse,platform
p=argparse.ArgumentParser();p.add_argument('task');p.add_argument('slot');p.add_argument('--preflight',action='store_true');a=p.parse_args();task=a.task;slot=a.slot
access=guard_access(OUT/f'logs/access_{task}_{slot}_python.jsonl')
if not a.preflight:assert (OUT/'PRE_EVAL_FREEZE.json').exists(),'Freeze real inputs before production scoring'
if not a.preflight:assert (OUT/f'config/TASK_INPUT_FREEZE_{task}.json').exists()
run=OUT/('runs/preflight' if a.preflight else 'runs')/task/slot;run.mkdir(parents=True,exist_ok=True)
cmd=['Rscript',str(OUT/'code/score_R.R'),str(OUT),task,slot]+(['preflight'] if a.preflight else [])
started=now()
if slot=='Integrative-1' and task in ('betapanc','dopaminergic'):
 if not a.preflight:complete_run(task,slot,task_status='UNSUPPORTED_TARGET',reason='official_mouse_target_GRN_category_absent; no substitution with neuron/pancreas',extra={'start_UTC':started,'n_source_samples':0,'n_target_samples':0,'source_rna_consumed':False,'target_rna_consumed':False,'external_target_reference':False})
 exit()
with open(OUT/f'logs/scoring_{task}_{slot}{"_preflight" if a.preflight else ""}.log','w') as log:code=subprocess.run(cmd,stdout=log,stderr=log).returncode
if a.preflight:
 raw=run/'RAW_METHOD_SCORES.tsv';n=0
 if raw.exists():n=int(np.isfinite(pd.to_numeric(pd.read_csv(raw,sep='\t').raw_score,errors='coerce')).sum())
 dump(f'runs/preflight/{task}/{slot}/CHECK.json',{'task':task,'slot':slot,'start_UTC':started,'end_UTC':now(),'command':cmd,'exit_code':code,'finite_score_count':n,'ranking_or_reference_match_inspected':False});print('Preflight',task,slot,'exit',code,'finite',n);exit(code)
if code:
 complete_run(task,slot,task_status='FAILED_EXECUTION',reason=f'native_R_exit_{code}; see logs/scoring_{task}_{slot}.log',extra={'command':cmd,'start_UTC':started,'exit_code':code});exit(code)
raw=pd.read_csv(run/'RAW_METHOD_SCORES.tsv',sep='\t').set_index('gene_id');raw=raw[~raw.index.isna()];assert not raw.index.duplicated().any()
for col in ['raw_score','ranking_score']:raw[col]=pd.to_numeric(raw[col],errors='coerce')
raw['scorable']=np.isfinite(raw.raw_score);raw['eligible_for_recommendation']=raw.scorable;raw['candidate_status']=np.where(raw.scorable,'ranked','undefined_score');raw['missing_reason']=np.where(raw.scorable,'','native_nonfinite_score')
inputs=[OUT/'config/CANDIDATE_UNIVERSE.tsv',OUT/'code/score_R.R']
if slot=='B':
 coverage=pd.read_csv(run/'REGULON_COVERAGE.tsv',sep='\t').set_index('gene_id');raw=raw.reindex(coverage.index).join(coverage)
 raw['scorable']=np.isfinite(raw.raw_score);raw['eligible_for_recommendation']=raw.scorable
 raw['candidate_status']=np.where(raw.scorable,'ranked','undefined_score');raw['missing_reason']=np.where(raw.scorable,'',np.where(raw.n_input_targets<25,'regulon_input_targets_below_minsize25','native_nonfinite_or_pleiotropy_filtered'))
 inputs +=[OUT/f'inputs/rna/{task}/DESEQ2_TARGET_V_SOURCE.tsv',OUT/'inputs/priors/COLLECTRI_MOUSE_M24.tsv'];consumed={'source_rna_consumed':True,'target_rna_consumed':True,'external_target_reference':False}
elif slot=='Integrative-2':
 pre=pd.read_csv(run/'pre_filter.tsv',sep='\t').set_index('gene');post=set(raw.index);raw=pre.copy();raw.index.name='gene_id';raw['raw_score']=raw.sum_rank;raw['ranking_score']=-raw.sum_rank;raw['scorable']=np.isfinite(raw.raw_score);raw['eligible_for_recommendation']=raw.index.isin(post)&raw.scorable
 raw['candidate_status']=np.where(raw.eligible_for_recommendation,'ranked',np.where(raw.mean_source_expr>20,'native_filtered','undefined_score'));raw['missing_reason']=np.where(raw.eligible_for_recommendation,'',np.where(raw.mean_source_expr>20,'mean_source_estimated_counts_gt20','native_nonfinite_score'));raw['direction']='unsigned_DE_and_network_priority'
 inputs +=[OUT/f'inputs/rna/{task}/ROUNDED_COUNTS.tsv',OUT/f'inputs/rna/{task}/SAMPLE_METADATA.tsv',OUT/'inputs/priors/COMPFORCE_MOUSE_NETWORK.tsv',OUT/'code/upstream/pm_mogrify_export_adapter.R'];consumed={'source_rna_consumed':True,'target_rna_consumed':True,'external_target_reference':False,'joint_redundancy_pruning_used':False}
 net=pd.read_csv(OUT/'inputs/priors/COMPFORCE_MOUSE_NETWORK.tsv',sep='\t');nodes=set(net.TF)|set(net.TG);cnt=pd.read_csv(OUT/f'inputs/rna/{task}/ROUNDED_COUNTS.tsv',sep='\t',index_col=0);nonzero=set(cnt.index[cnt.sum(axis=1)>0]);uu=pd.read_csv(OUT/'config/CANDIDATE_UNIVERSE.tsv',sep='\t')
 for g in uu.gene_id:
  if g in nodes and g not in raw.index:
   raw.loc[g,['raw_score','ranking_score']]=np.nan;raw.loc[g,['scorable','eligible_for_recommendation']]=False;raw.loc[g,'candidate_status']='input_absent';raw.loc[g,'missing_reason']='all_counts_zero_or_gene_absent' if g not in nonzero else 'network_disconnected_after_native_DE_gene_intersection'
else:
 inputs +=[OUT/'inputs/rna/source/TPM.tsv',OUT/'inputs/annotation/SYMBOL_TO_STABLE_ID.json',OUT/'code/upstream/pacnet_utils.R',OUT/'code/upstream/ccn_grn_utils.R']+[AUDIT/'sources/cellnet_resources'/f for f in ['Mm_grnAll_curatedTFs_Apr-22-2020.rda','Mm_trainingNormalization_Apr-22-2020.rda','Mm_expTrain_Oct-24-2016.rda']];consumed={'source_rna_consumed':True,'target_rna_consumed':False,'external_target_reference':True}
 priorTFs=pd.read_csv(run/'ORIGINAL_TARGET_GRN_TFS.tsv',sep='\t').gene_id.dropna()
 for g in priorTFs:
  if g not in raw.index:
   raw.loc[g,['raw_score','ranking_score']]=np.nan;raw.loc[g,['scorable','eligible_for_recommendation']]=False;raw.loc[g,'candidate_status']='input_absent';raw.loc[g,'missing_reason']='native_target_GRN_TF_unavailable_after_query_training_gene_alignment'
params=json.loads((run/'NATIVE_PARAMETERS.json').read_text());raw['scorable']=raw.scorable.fillna(False);raw['eligible_for_recommendation']=raw.eligible_for_recommendation.fillna(False)
complete_run(task,slot,raw,inputs={str(p):digest(p) for p in inputs},extra={'command':cmd,'start_UTC':started,'exit_code':code,'native_parameters':params,**consumed})
