args <- commandArgs(trailingOnly=TRUE)
root <- args[[1L]]
target <- args[[2L]]
suppressPackageStartupMessages(library(compforce))
suppressPackageStartupMessages(library(jsonlite))
options(digits=17)
set.seed(20261118L)
frozen <- file.path(root,"FROZEN_CONFIG_ROUND2.yaml")
file_sha <- function(path) strsplit(system2("sha256sum",shQuote(path),stdout=TRUE)," +")[[1L]][[1L]]
cfg_hash <- file_sha(frozen)
config <- fromJSON(frozen,simplifyVector=FALSE)
stopifnot(target %in% unlist(config$targets,use.names=FALSE))
runner <- file.path(root,"scripts/compforce_run.R")
runner_sha <- file_sha(runner)
out <- file.path(root,"results",target,"CompForce")
checkpoint <- file.path(root,"checkpoints","CompForce",target,cfg_hash,runner_sha)
dir.create(checkpoint,recursive=TRUE,showWarnings=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
record <- function(path,data) {
  temporary <- paste0(path,".",Sys.getpid(),".partial")
  write_json(data,temporary,auto_unbox=TRUE,pretty=TRUE,digits=17,na="null")
  stopifnot(file.rename(temporary,path))
}
table <- function(x,filename) {
  path <- file.path(out,filename)
  tmp <- paste0(path,".partial")
  write.table(x,tmp,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
  stopifnot(file.rename(tmp,path))
  path
}

base <- config$paths$previous_frozen_inputs
count_path <- file.path(base,"inputs","rna",target,"ROUNDED_COUNTS.tsv")
network_path <- file.path(root,"inputs_manifest/COMPFORCE_MOUSE_CHIPX_M24.tsv")
U_path <- file.path(root,"inputs_manifest/CANDIDATE_UNIVERSE.tsv")
raw <- read.delim(count_path,check.names=FALSE)
ids <- raw[[1L]]
counts <- as.matrix(raw[,-1L,drop=FALSE])
rownames(counts) <- ids
condition <- config$samples[[target]]
source_ids <- unlist(condition$rna_source,use.names=FALSE)
target_ids <- unlist(condition$rna_target,use.names=FALSE)
stopifnot(identical(colnames(counts),c(source_ids,target_ids)),nrow(counts)==54331L,
          all(is.finite(counts)),all(counts==round(counts)), min(counts)>=0,
          length(source_ids)==3L,length(target_ids)==2L,!anyDuplicated(ids))
meta <- data.frame(condition=rep(c("source",target),c(length(source_ids),length(target_ids))),row.names=colnames(counts))
U <- read.delim(U_path)$gene_id
network <- read.delim(network_path,stringsAsFactors=FALSE)[,c("TF","TG")]
stopifnot(length(U)==1374L,!anyDuplicated(U),nrow(network)==444007L,
          sum(network$TF %in% ids & network$TG %in% ids)>0L)
obj <- compforce::cf_initialize(expCounts=counts,sampTab=meta,start="source",target=target,
                                species="Mm",description_column="condition",TFs=U,network=network)
trans <- paste("source",target,sep="..")
network_tf <- tabulate(match(network$TF,U), nbins=length(U))
network_tg <- tabulate(match(network$TG,U), nbins=length(U))
table(data.frame(TF=U,count_present=U %in% ids,source_mean_raw_counts=rowMeans(counts[,source_ids,drop=FALSE])[U],
                 target_mean_raw_counts=rowMeans(counts[,target_ids,drop=FALSE])[U],
                 expressed_row=rowSums(counts)[U]>0,chipX_prior_out_edges=network_tf,chipX_prior_in_edges=network_tg),
      "COMPFORCE_CANDIDATE_INPUT.tsv")
table(data.frame(sample_id=colnames(counts),condition=meta$condition,library_total=colSums(counts),
                 is_background=meta$condition=="source",scale="rounded_tximport_lengthScaledTPM_estimated_counts"),
      "COMPFORCE_OBSERVATIONS.tsv")

# Insert side-effect-only observers around the archived component statements. They do not alter their arguments or results.
observe <- function(fun, before, expression) {
  code <- as.list(body(fun))
  lines <- vapply(code,function(line)paste(deparse(line),collapse=" "),character(1L))
  location <- grep(before,lines,fixed=TRUE)
  if (length(location)>1L) location <- which(startsWith(lines,before))
  stopifnot(length(location)==1L)
  body(fun) <- as.call(append(code,list(expression),after=location-1L))
  fun
}
diffex <- observe(compforce::pm_diffexp,"topTFs <- na.omit",quote(get("observe_de",envir=.GlobalEnv)(topTFs,method)))
cellnet <- observe(compforce::pm_cellnet,"nis_res <- adapted_NIS",quote(get("observe_nis_inputs",envir=.GlobalEnv)(trained_nets,expNorm,sampTab,target)))
cellnet <- observe(cellnet,"nis_summary <- summarize_NIS",quote(get("observe_nis_scores",envir=.GlobalEnv)(nis_res)))
mogrify <- observe(compforce::pm_mogrify,"TFscores <- TFscores[!(TFscores$mean_source_expr",quote(get("observe_mogrify_input",envir=.GlobalEnv)(TFscores)))
mogrify <- observe(mogrify,"TFscores <- TFscores[!(TFscores$gene %in% redundant)",quote(get("observe_mogrify_redundancy",envir=.GlobalEnv)(redundant)))
prefix <- "main"
observe_de <- function(x,method) invisible(table(x,paste0(prefix,"_DiffExp_",method,"_BEFORE_PADJ.tsv")))
observe_nis_inputs <- function(nets,expression,samples,task) {
  grn <- nets$ctGRNs[[task]]
  g <- union(grn$TF,grn$TG)
  sd <- apply(expression[g,samples$cf_annotation==task,drop=FALSE],1,sd)
  tfs <- sort(unique(grn$TF))
  degrees <- lengths(lapply(split(grn$TG,grn$TF),unique))
  d <- data.frame(TF=tfs,target_n=as.integer(degrees[tfs]),target_variance_0=as.logical(sd[tfs]==0))
  invisible(table(d,"CellNet_TARGET_SUBNETWORK_TFS.tsv"))
  invisible(table(data.frame(stable_id=g,target_sd=as.numeric(sd)),"CellNet_TARGET_SUBNETWORK_GENES.tsv"))
  record(file.path(out,"CellNet_SUBNETWORK_QC.json"),list(target=task,overall_n=nrow(nets$overallGRN),target_edges=nrow(grn),
         target_nodes=length(g),target_regulators=length(tfs),target_genes_zero_variance=sum(sd==0,na.rm=TRUE),
         target_genes_nonfinite_variance=sum(!is.finite(sd)),target_specific_genes=length(nets$specGenes[[task]])))
}
observe_nis_scores <- function(x) invisible(table(data.frame(TF=rownames(x),x,check.names=FALSE),"CellNet_NIS_PER_SOURCE_SAMPLE.tsv"))
observe_mogrify_input <- function(x) invisible(table(x,"Mogrify_BEFORE_NATIVE_FILTERS.tsv"))
observe_mogrify_redundancy <- function(x) invisible(table(data.frame(TF=as.character(x)),"Mogrify_REDUNDANT_TFS.tsv"))

main_calls <- list(DiffExp_ttest=function()diffex(obj,method="ttest",transition_row=1,pThresh=0.05),
                  DiffExp_deseq=function()diffex(obj,method="deseq",transition_row=1,pThresh=0.05),
                  CellNet=function()cellnet(obj,use_cfobj_network=TRUE,transition_row=1),
                  Mogrify=function()mogrify(obj,use_cfobj_network=TRUE,background_method="source",transition_row=1,diff_thresh=20),
                  dAlessio=function()compforce::pm_dalessio(obj,all_background=FALSE,transition_row=1))
native_metric <- c(DiffExp_ttest="padj",DiffExp_deseq="padj",CellNet="medVal",Mogrify="sum_rank",dAlessio="specificity_score")

run_component <- function(label,key,call) {
  cache <- file.path(checkpoint,paste0(label,".rds"))
  successful <- c("completed","completed_empty","completed_partial_statistic")
  if (file.exists(cache)) {
    saved <- readRDS(cache)
    if (saved$state %in% successful) {
      stopifnot(saved$runner_sha256==runner_sha,saved$config_sha256==cfg_hash)
      if (all(vapply(saved$file_sha256,function(entry) file.exists(entry$path)&&file_sha(entry$path)==entry$sha256,logical(1L)))) {
        cat("Verified resumed component",label,"\n")
        return(saved)
      }
    }
  }
  cat("Starting",target,label,format(Sys.time(),tz="Asia/Singapore",usetz=TRUE),"\n")
  start <- proc.time()[["elapsed"]]
  warnings <- character()
  value <- tryCatch(withCallingHandlers(call(),warning=function(w) warnings <<- unique(c(warnings,conditionMessage(w)))),error=identity)
  result <- list(label=label,key=key,started=TRUE,seconds=proc.time()[["elapsed"]]-start,
                 config_sha256=cfg_hash,runner_sha256=runner_sha,warnings=as.list(warnings))
  if (inherits(value,"error")) {
    result$state <- "failed_execution"
    result$reason <- conditionMessage(value)
    result$error_class <- as.list(class(value))
    cat("COMPONENT_ERROR",label,result$reason,"\n")
  } else {
    result$raw <- value@results_raw[[trans]][[key]]
    result$clean <- value@results[[trans]][[key]]
    stopifnot(is.data.frame(result$clean),is.data.frame(result$raw),
              all(c("TF","rank") %in% names(result$clean)),all(result$clean$TF %in% U),!anyDuplicated(result$clean$TF),
              identical(as.character(result$raw$TF),as.character(result$clean$TF)),
              identical(as.integer(result$clean$rank),seq_len(nrow(result$clean))))
    result$number <- nrow(result$clean)
    result$number_finite <- sum(is.finite(result$raw[[native_metric[[key]]]]))
    result$state <- if(!result$number) "completed_empty" else if(!result$number_finite) "undefined_statistic" else
                        if(result$number_finite<result$number) "completed_partial_statistic" else "completed"
    if (result$state=="undefined_statistic") result$reason <- "No finite metric in native component list"
    a <- table(result$raw,paste0(label,"_NATIVE_RAW.tsv"))
    b <- table(result$clean,paste0(label,"_NATIVE_LIST.tsv"))
    result$file_sha256 <- lapply(c(a,b),function(p)list(path=p,sha256=file_sha(p)))
  }
  temp <- paste0(cache,".partial")
  saveRDS(result,temp,compress=TRUE,version=3L)
  stopifnot(file.rename(temp,cache))
  record(file.path(out,paste0(label,"_RUN.json")),result[setdiff(names(result),c("clean","raw"))])
  cat("Finished",target,label,result$state,"\n")
  result
}

results <- lapply(names(main_calls),function(key)run_component(key,key,main_calls[[key]]))
names(results) <- names(main_calls)
prefix <- "sensitivity"
convert <- getFromNamespace("downsample_and_transform","compforce")
variant <- run_component("DiffExp_ttest_logCP100k_sensitivity","DiffExp_ttest",function(){
  log_counts <- convert(counts,xFact=1e5,observations_by_features=FALSE)
  stopifnot(identical(dimnames(log_counts),dimnames(counts)),all(is.finite(log_counts)))
  variant_input <- obj
  variant_input@expression_matrix <- log_counts
  diffex(variant_input,method="ttest",transition_row=1,pThresh=0.05)
})

aggregate_complete <- function(lists,label) {
  usable <- vapply(lists,function(l)l$state %in% c("completed","completed_empty","completed_partial_statistic"),logical(1L))
  nonempty <- sum(vapply(lists,function(l)!is.null(l$number)&&l$number>0L,logical(1L)))
  report <- list(name=label,number_of_components=length(lists),native_nonempty_list_n=nonempty,
                 inputs=lapply(lists,function(x)x[intersect(names(x),c("label","state","reason","number","number_finite","file_sha256"))]))
  if (!all(usable)) {
    report$execution_status <- "blocked_component"
    report$reason <- paste(paste(names(lists)[!usable],vapply(lists[!usable],function(l)l$state,character(1L)),sep="="),collapse="; ")
    return(report)
  }
  if (!nonempty) {
    report$execution_status <- "completed_empty"
    report$reason <- "Five native lists empty; candidate union empty"
    path <- table(data.frame(TF=character(),RRA_Score=double(),rank=integer()),paste0(label,"_NATIVE_AGGREGATE.tsv"))
    report$native_candidate_n <- 0L
    report$output <- list(path=path,sha256=file_sha(path))
    return(report)
  }
  x <- obj
  x@results[[trans]] <- lapply(lists,`[[`,"clean")
  x@results_raw[[trans]] <- lapply(lists,`[[`,"raw")
  native <- tryCatch(compforce::cf_RRA(x),error=identity)
  if (inherits(native,"error")) {
    report$execution_status <- "failed_execution"
    report$reason <- conditionMessage(native)
    return(report)
  }
  score <- native@results[[trans]]$aggregate
  tf_union <- unique(unlist(lapply(x@results[[trans]],`[[`,"TF"),use.names=FALSE))
  stopifnot(is.data.frame(score),identical(names(score),c("TF","RRA_Score","rank")),
            setequal(score$TF,tf_union),!anyDuplicated(score$TF),
            all(is.finite(score$RRA_Score)),all(score$RRA_Score>=0),all(score$RRA_Score<=1),
            all(diff(score$RRA_Score)>=0))
  path <- table(score,paste0(label,"_NATIVE_AGGREGATE.tsv"))
  report$execution_status <- if(length(unique(score$RRA_Score))==1L)"completed_uninformative" else "completed"
  report$native_candidate_n <- nrow(score)
  report$output <- list(path=path,sha256=file_sha(path))
  report$empty_names_retained <- as.list(names(lists)[vapply(lists,function(l)l$number==0L,logical(1L))])
  report
}

primary <- aggregate_complete(results,"CompForce_RRA5_archived")
alt <- results
alt$DiffExp_ttest <- variant
other <- setdiff(names(results),"DiffExp_ttest")
same_bytes <- setNames(lapply(other,function(k){
  x <- results[[k]]
  a <- serialize(x$clean,NULL,version=3L); b <- serialize(alt[[k]]$clean,NULL,version=3L)
  f <- x$file_sha256
  available <- length(f)>0L
  list(in_memory_identical=identical(a,b),component_files_available=available,
       export_sha256_verified_if_available=if(available)all(vapply(f,function(z)file_sha(z$path)==z$sha256,logical(1L))) else NA,
       preserved_same_export_paths=identical(f,alt[[k]]$file_sha256))
}),other)
stopifnot(all(vapply(same_bytes,function(x)x$in_memory_identical&&x$preserved_same_export_paths,logical(1L))))
sensitivity <- aggregate_complete(alt,"CompForce_RRA5_ttest_logCP100k_sensitivity")
record(file.path(out,"METHOD_RUN.json"),list(target=target,config_sha256=cfg_hash,runner_sha256=runner_sha,
       implementation="CompForce_RRA5_archived",execution_status=primary$execution_status,
       reason=if(is.null(primary$reason)) "" else primary$reason,
       primary=primary,sensitivity=sensitivity,identical_other_four_components=same_bytes,
       count_scale="rounded tximport lengthScaledTPM estimated counts", full_counts_n=nrow(counts),
       network_overlap_counts=list(archive_network_rows=nrow(network),
                                   both_endpoints_in_counts=sum(network$TF %in% ids & network$TG %in% ids),
                                   regulator_missing_count_row=sum(!network$TF %in% ids),
                                   target_missing_count_row=sum(!network$TG %in% ids)),
       software=list(CompForce=as.character(packageVersion("compforce")),DESeq2=as.character(packageVersion("DESeq2")),
                     CellNet=as.character(packageVersion("CellNet")),RobustRankAggreg=as.character(packageVersion("RobustRankAggreg"))),
       mogrify_observed_previous_only_R_expression=as.list(5L-1:1)))
cat(target,primary$execution_status,sensitivity$execution_status,"\n")
