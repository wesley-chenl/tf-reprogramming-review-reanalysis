args<-commandArgs(TRUE);root<-normalizePath(args[1]);task<-args[2];slot<-args[3]
set.seed(1729)
library(jsonlite);library(igraph)
preflight<-length(args)>=4&&args[4]=='preflight';out<-if(preflight)file.path(root,'runs/preflight',task,slot) else file.path(root,'runs',task,slot);dir.create(out,recursive=TRUE,showWarnings=FALSE)
access<-file(file.path(root,'logs',paste0('access_',task,'_',slot,'_R.jsonl')),open='wt')
audit_file<-function(p){if(!is.character(p)||!length(p))return(invisible(NULL));s<-gsub('\\\\','/',p[1]);deny<-grepl('/(references|evaluation)/|/config/TASKS.tsv|Taiji-reprogram-main/benchmark/',s,ignore.case=TRUE)||grepl('gemini_independent_audit_',s)&&!grepl('/sources/',s);writeLines(toJSON(list(path=p[1],evaluation_access_denied=deny),auto_unbox=TRUE),access);if(deny)stop('Evaluation input forbidden in real R scorer')}
trace('file',where=baseenv(),tracer=quote(audit_file(description)),print=FALSE)
trace('gzfile',where=baseenv(),tracer=quote(audit_file(description)),print=FALSE)
on.exit(close(access))
u<-read.delim(file.path(root,'config/CANDIDATE_UNIVERSE.tsv'),stringsAsFactors=FALSE)
readmat<-function(p)as.matrix(read.delim(p,row.names=1,check.names=FALSE))
write_scores<-function(x){for(nm in colnames(x))if(is.numeric(x[[nm]]))x[[nm]]<-ifelse(is.na(x[[nm]]),NA_character_,sprintf('%.17g',x[[nm]]));write.table(x,file.path(out,'RAW_METHOD_SCORES.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')}
if(slot=='B'){
 library(viper)
 source(file.path(root,'code/upstream/VIPER_dimension_safe.R'))
 de<-read.delim(file.path(root,'inputs/rna',task,'DESEQ2_TARGET_V_SOURCE.tsv'),stringsAsFactors=FALSE)
 prior<-read.delim(file.path(root,'inputs/priors/COLLECTRI_MOUSE_M24.tsv'),stringsAsFactors=FALSE)
 signature<-matrix(de$stat[is.finite(de$stat)],ncol=1,dimnames=list(de$gene_id[is.finite(de$stat)],'target_vs_source_Wald'))
 regulon<-lapply(split(prior,prior$source),function(e)list(tfmode=setNames(e$weight,e$target),likelihood=setNames(rep(1,nrow(e)),e$target)))
 tf<-sort(unique(prior$source));coverage<-data.frame(gene_id=tf,n_prior_targets=sapply(regulon,function(x)length(x$tfmode)),n_input_targets=sapply(regulon,function(x)sum(names(x$tfmode)%in%rownames(signature))))
 coverage$minsize_pass<-coverage$n_input_targets>=25
 activity<-viper_dimension_safe(signature,regulon,method='none',nes=TRUE,pleiotropy=TRUE,minsize=25,bootstraps=0,verbose=TRUE)
 if(preflight){
  original_two_columns<-viper(cbind(signature,signature),regulon,method='none',nes=TRUE,pleiotropy=TRUE,minsize=25,bootstraps=0,verbose=FALSE)
  common<-intersect(rownames(activity),rownames(original_two_columns));difference<-max(abs(activity[common,1]-original_two_columns[common,1]));stopifnot(is.finite(difference),difference<1e-8)
  write_json(list(interface_patch='drop=FALSE in default eset.filter TRUE selection; preserves full native statistical background',single_column_vs_unmodified_identical_two_column_max_abs_difference=difference,n_TFs_compared=length(common),duplicate_column_used_only_for_shape_validation_not_biological_replication=TRUE),file.path(out,'SINGLE_COLUMN_NATIVE_EQUIVALENCE.json'),pretty=TRUE,auto_unbox=TRUE)
 }
 saveRDS(activity,file.path(out,'NATIVE_VIPER_OUTPUT.rds'))
 result<-data.frame(gene_id=rownames(activity),raw_score=activity[,1],ranking_score=activity[,1],direction='positive_target_associated_activity')
 coverage$finite_score<-coverage$gene_id%in%result$gene_id[is.finite(result$raw_score)]
 write.table(coverage,file.path(out,'REGULON_COVERAGE.tsv'),sep='\t',row.names=FALSE,quote=FALSE);write_scores(result)
 write_json(list(viper_version=as.character(packageVersion('viper')),parameters=list(method='none',nes=TRUE,pleiotropy=TRUE,minsize=25,bootstraps=0,likelihood=1),input_signature='DESeq2 signed Wald target/source',n_finite_signature=nrow(signature)),file.path(out,'NATIVE_PARAMETERS.json'),pretty=TRUE,auto_unbox=TRUE)
}
