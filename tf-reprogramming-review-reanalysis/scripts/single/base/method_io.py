from common import *
import numpy as np,pandas as pd
TASKS=['liver','cardiomyocyte','betapanc','dopaminergic','mESC']
RUN_BEGIN=now()
SLOTS=['A','B','C','D','Integrative-1','Integrative-2']
METHODS={'A':('TRIAGE','original_discordance_score'),'B':('VIPER','VIPER_CollecTRI_mouse_adapter'),'C':('Taiji PageRank-ratio module','Taiji_1.3.0'),'D':('TFActProfiler local response','target_relative_paired_sham_HCOP_adapter'),'Integrative-1':('PACNet NIS','official_successor_NIS_relaWeight0'),'Integrative-2':('CompForce Mogrify-style single TF','pm_mogrify_pre_redundancy_export_adapter')}
def guard_access(path):
 """Trace real Python file opens and prohibit evaluation dependencies."""
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);log=open(path,'w',buffering=1);seen=set()
 def audit(event,args):
  if event!='open' or not isinstance(args[0],(str,bytes,os.PathLike)):return
  p=os.path.abspath(os.fsdecode(args[0]));s=p.replace('\\','/').lower()
  denied=any(x in s for x in ['/references/','/evaluation/','/config/tasks.tsv','/taiji-reprogram-main/benchmark/']) or ('gemini_independent_audit_' in s and '/sources/' not in s)
  if p not in seen:seen.add(p);log.write(json.dumps({'path':p,'evaluation_access_denied':denied})+'\n')
  if denied:raise PermissionError('Evaluation/reference input prohibited in scoring runner: '+p)
 sys.addaudithook(audit)
 return log
def complete_run(task,slot,scores=None,task_status=None,reason='',inputs=None,components=None,extra=None):
 d=OUT/'runs'/task/slot;d.mkdir(parents=True,exist_ok=True)
 u=pd.read_csv(OUT/'config/CANDIDATE_UNIVERSE.tsv',sep='\t',keep_default_na=False)
 mid,impl=METHODS[slot]
 c_config=OUT/'config/FROZEN_CONFIG_C.json'
 if slot=='C' and c_config.exists():
  cc=json.loads(c_config.read_text());mid=cc.get('method_id',mid);impl=cc.get('implementation_id',impl)
 out=u[['gene_id','gene_symbol']].copy();out['raw_score']=np.nan;out['ranking_score']=np.nan;out['scorable']=False;out['eligible_for_recommendation']=False;out['candidate_status']='prior_absent';out['missing_reason']='not_in_method_prior';out['direction']=''
 if scores is not None and len(scores):
  assert not scores.index.duplicated().any()
  ix=out.gene_id.isin(scores.index);cols=list(scores.columns)
  for col in cols:
   if col not in out:out[col]=pd.Series([None]*len(out),index=out.index,dtype=object) # mixed/nullable native components
   out.loc[ix,col]=out.loc[ix,'gene_id'].map(scores[col]).to_numpy()
 for c in ['scorable','eligible_for_recommendation']:out[c]=out[c].fillna(False).astype(bool)
 if task_status is not None:
  out['candidate_status']='not_run_due_to_task_failure';out['missing_reason']=reason
 unm=u.mapping_status.eq('id_unresolved');out.loc[unm,['candidate_status','missing_reason']]=['id_unresolved','official_TF_ID_absent_GENCODE_M24']
 out.loc[unm,['scorable','eligible_for_recommendation']]=False
 out.loc[unm,['raw_score','ranking_score']]=np.nan
 out['eligible_for_recommendation']=out['eligible_for_recommendation']&np.isfinite(pd.to_numeric(out.ranking_score,errors='coerce'))
 eligible=out[out.eligible_for_recommendation].sort_values(['ranking_score','gene_id'],ascending=[False,True],kind='stable')
 for c in ['rank','rank_min','rank_max','tie_group']:out[c]=np.nan
 out['tie_crosses_10']=False;out['in_top10']=False
 if len(eligible):
  eligible['rank']=np.arange(1,len(eligible)+1)
  eligible['rank_min']=eligible.groupby('ranking_score')['rank'].transform('min');eligible['rank_max']=eligible.groupby('ranking_score')['rank'].transform('max')
  eligible['tie_group']=pd.factorize(eligible.ranking_score,sort=False)[0]+1
  eligible['tie_crosses_10']=(eligible.rank_min<=10)&(eligible.rank_max>10)
  for c in ['rank','rank_min','rank_max','tie_group','tie_crosses_10']:out.loc[eligible.index,c]=eligible[c]
  out.loc[eligible.head(10).index,'in_top10']=True
 if task_status is None:
  if len(eligible)<10:task_status='UNDERFILLED';reason='fewer_than_10_native_eligible_finite_candidates'
  elif eligible.ranking_score.nunique()==1:task_status='DEGENERATE_OUTPUT';reason='all_native_eligible_scores_identical';out['in_top10']=False
  else:task_status='REAL_SCORING_COMPLETE'
 cfg=OUT/'config/FROZEN_CONFIG.yaml';ch=digest(cfg) if cfg.exists() else 'preflight_not_frozen'
 if slot=='C' and c_config.exists():ch=hashlib.sha256((ch+digest(c_config)).encode()).hexdigest()
 ih=hashlib.sha256(json.dumps(inputs or {},sort_keys=True).encode()).hexdigest()
 for col,value in [('task_id',task),('slot_id',slot),('method_id',mid),('implementation_id',impl),('config_hash',ch),('input_hash',ih)]:out[col]=value
 order=['task_id','slot_id','method_id','implementation_id','gene_id','gene_symbol','raw_score','ranking_score','rank','rank_min','rank_max','scorable','eligible_for_recommendation','in_top10','candidate_status','missing_reason','direction','input_hash','config_hash']
 out=out[order+[c for c in out if c not in order]]
 out.to_csv(d/'FULL_U_SCORES.tsv',sep='\t',index=False,na_rep='NA',float_format='%.17g')
 code_paths=[Path(__file__),OUT/'code/common.py',Path(sys.modules['__main__'].__file__)]
 smp=OUT/f'inputs/rna/{task}/SAMPLE_METADATA.tsv';meta=pd.read_csv(smp,sep='\t') if smp.exists() else pd.DataFrame({'group':[]})
 rec={'task_id':task,'slot_id':slot,'method_id':mid,'implementation_id':impl,'config_hash':ch,'input_hashes':inputs or {},'code_hashes':{str(p):digest(p) for p in code_paths},'command':[sys.executable]+sys.argv,'start_UTC':RUN_BEGIN,'end_UTC':now(),'exit_code':0 if scores is not None else None,'status':task_status,'reason':reason,'n_candidates':len(u),'n_source_samples':3 if slot=='Integrative-1' and task_status=='REAL_SCORING_COMPLETE' else int(meta.group.eq('source').sum()),'n_target_samples':0 if slot=='Integrative-1' else int(meta.group.eq(task).sum()),'n_scorable':int(out.scorable.sum()),'n_native_eligible':len(eligible),'shortlist_size':int(out.in_top10.sum()),'output_path':str(d/'FULL_U_SCORES.tsv'),'output_sha256':digest(d/'FULL_U_SCORES.tsv'),'real_input':scores is not None,'evaluation_used':False,'environment_manifest':str(OUT/'environments/VERSIONS.json')}
 if extra:rec.update(extra)
 dump(f'runs/{task}/{slot}/RUN_RECORD.json',rec)
 print(task,slot,task_status,'n_finite',rec['n_scorable'],'n_eligible',len(eligible),flush=True)
 return rec
