from method_io import *
import argparse
p=argparse.ArgumentParser();p.add_argument('task');p.add_argument('--isolation-recheck',action='store_true');a=p.parse_args();task=a.task
assert (OUT/'PRE_EVAL_FREEZE.json').exists(),'Freeze real inputs before production scoring'
assert (OUT/f'config/TASK_INPUT_FREEZE_{task}.json').exists()
access=guard_access(OUT/f'logs/access_{task}_A{"_recheck" if a.isolation_recheck else ""}.jsonl')
tpm=OUT/f'inputs/rna/{task}/TARGET_TPM.tsv';st=OUT/f'inputs/rna/{task}/SAMPLE_METADATA.tsv';rp=OUT/'inputs/priors/TRIAGE_RTS_M24.tsv'
x=pd.read_csv(tpm,sep='\t',index_col=0);meta=pd.read_csv(st,sep='\t');cols=meta.loc[meta.group==task,'sample_id'].tolist();x=x[cols]
assert np.isfinite(x.to_numpy()).all() and (x.to_numpy()>=0).all()
rts=pd.read_csv(rp,sep='\t',index_col=0).RTS;ids=x.index.intersection(rts.index)
ds=np.log1p(x.loc[ids]).mul(rts.loc[ids],axis=0);s=ds.mean(axis=1)
(OUT/f'runs/{task}').mkdir(parents=True,exist_ok=True)
if not a.isolation_recheck:ds.to_csv(OUT/f'runs/{task}/A_REPLICATE_DISCORDANCE.tsv',sep='\t')
scores=pd.DataFrame({'raw_score':s,'ranking_score':s,'scorable':True,'eligible_for_recommendation':True,'candidate_status':'ranked','missing_reason':'','direction':'target-associated discordance; not an OE effect'})
scores=scores.reindex(rts.index);absent=~scores.index.isin(x.index);scores.loc[absent,['scorable','eligible_for_recommendation']]=False;scores.loc[absent,['candidate_status','missing_reason']]=['input_absent','gene_absent_from_uniform_target_TPM']
if a.isolation_recheck:
 old=pd.read_csv(OUT/f'runs/{task}/A/FULL_U_SCORES.tsv',sep='\t').set_index('gene_id')
 v=old
 ii=v.index.intersection(scores.index);eq=np.allclose(v.loc[ii,'raw_score'],scores.loc[ii,'raw_score'],rtol=1e-14,atol=1e-14,equal_nan=True)
 dump('reports/REAL_RUNNER_ISOLATION_RECHECK.json',{'task':task,'slot':'A','real_input_recomputed':True,'same_scores_with_reference_access_prohibited':bool(eq),'n_compared':len(ii),'access_log':str(OUT/f'logs/access_{task}_A_recheck.jsonl')});assert eq
else:complete_run(task,'A',scores,inputs={str(p):digest(p) for p in [tpm,st,rp]},extra={'target_samples':cols,'n_source_samples':0,'source_rna_consumed':False,'target_rna_consumed':True,'formula':'mean biological replicate ln(TPM+1)*RTS','native_score_direction':'descending'})
