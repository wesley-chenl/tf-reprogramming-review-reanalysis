args <- commandArgs(trailingOnly=TRUE)
root <- args[[1]]
target <- args[[2]]
workers <- as.integer(args[[3]])
preflight <- identical(target,"preflight")
suppressPackageStartupMessages(library(IMCell))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(jsonlite))
if (!preflight) stopifnot(file.exists(file.path(root,"FROZEN_CONFIG_ROUND2.yaml")))
RNGkind("L'Ecuyer-CMRG")
constant_seed_kind <- c("L'Ecuyer-CMRG","Inversion","Rejection")
stopifnot(identical(RNGkind(),constant_seed_kind))
prior <- read.delim(file.path(root,"inputs_manifest/IMCELL_COLLECTRI_MOUSE_SIGNED.tsv"))
stopifnot(all(prior$weight>0), all(is.finite(prior$weight)), setequal(prior$type,c(1,-1)),
          !anyDuplicated(prior[,c("TF","TG")]))
graph <- igraph::graph_from_data_frame(prior, directed=TRUE)
program <- read.delim(file.path(root,"inputs_manifest/IMCELL_TARGET_PROGRAM.tsv"))
program <- program[tolower(program$mapped_to_graph)=="true",]
jobs <- read.delim(file.path(root,"inputs_manifest/IMCELL_RNG_STREAMS.tsv"))
coverage <- read.delim(file.path(root,"inputs_manifest/IMCELL_SEED_COVERAGE.tsv"))
record <- function(path,data) {
  temporary <- paste0(path,".",Sys.getpid(),".partial")
  write_json(data,temporary,auto_unbox=TRUE,pretty=TRUE,digits=17,na="null")
  stopifnot(file.rename(temporary,path))
}
file_sha <- function(path) {
  strsplit(system2("sha256sum",shQuote(path),stdout=TRUE)," +")[[1]][[1]]
}

if (!preflight) {
  adapt <- fromJSON(file.path(root,"checkpoints/IMCELL_RUNTIME_READY.json"))
  code <- file.path(root,"scripts/imcell_indexed_cpu.R")
  stopifnot(adapt$config_sha256==file_sha(file.path(root,"FROZEN_CONFIG_ROUND2.yaml")),
            adapt$indexed_code_sha256==file_sha(code))
  source(code)
  r2_install_indexed_mc()
  cat("IMCell archive rank_by_spread with exact ordered integer lookups:",adapt$indexed_code_sha256,"\n")
}

score_job <- function(job,niter=1000L) {
  assign(".Random.seed",as.integer(job[1L,paste0("rng_",1:7)]),envir=.GlobalEnv)
  targets <- program[program$target==job$target,,drop=FALSE]
  activate <- as.character(targets$stable_id[targets$direction=="activate"])
  repress <- as.character(targets$stable_id[targets$direction=="repress"])
  stopifnot(length(activate)+length(repress)>0L)
  call_time <- proc.time()[["elapsed"]]
  ans <- IMCell::rank_by_spread(
    ig=graph,repressor_wins=FALSE,tfs=as.character(job$stable_id),
    targets_activate=activate,targets_repress=repress,niter=niter,
    edge_probability_method="in_degree",edge_probability_multiplier=1,
    return_spread=TRUE,num_cores=1L)
  stopifnot(is.data.frame(ans),nrow(ans)==1L,identical(as.character(ans$TF),as.character(job$stable_id)),
            is.finite(ans$spread),ans$spread>=0,ans$spread<=length(activate)+length(repress))
  list(score=as.numeric(ans$spread), wall=proc.time()[["elapsed"]]-call_time)
}
if (preflight) {
  first <- jobs[jobs$target=="liver" & jobs$replicate==1L,]
  chosen <- head(first,2L)
  probe <- lapply(seq_len(nrow(chosen)),function(i) {
    job <- chosen[i,,drop=FALSE]
    # Test the exact native singleton call and explicit job state at full niter.
    a <- score_job(job); b <- score_job(job)
    stopifnot(identical(a$score,b$score))
    list(stable_id=job$stable_id,niter=1000L,spread=a$score,wall_s=a$wall,
         repeated_identical=TRUE,rng=as.integer(job[1,paste0("rng_",1:7)]))
  })
  forked <- mclapply(seq_len(nrow(chosen)),function(i) score_job(chosen[i,,drop=FALSE]),
                     mc.cores=min(2L,workers),mc.set.seed=FALSE)
  stopifnot(identical(lapply(probe,function(x)x$spread),lapply(forked,function(x)x$score)))
  record(file.path(root,"validation/IMCELL_PREFLIGHT.json"),
         list(chosen=probe,fork_equals_serial=TRUE,official_function="IMCell::rank_by_spread",
              policy="Lexicographically first two eligible stable IDs in fixed first task; reference independent",
              igraph=as.character(packageVersion("igraph")), R=as.character(getRversion())))
  quit(save="no",status=0L)
}

task_jobs <- jobs[jobs$target==target,,drop=FALSE]
out <- file.path(root,"results",target,"IMCell")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cfg_hash <- file_sha(file.path(root,"FROZEN_CONFIG_ROUND2.yaml"))
checkpoint <- file.path(root,"checkpoints","IMCell",target,cfg_hash)
dir.create(checkpoint,recursive=TRUE,showWarnings=FALSE)
fields <- c("job_id","target","stable_id","replicate","spread","seconds",
            "job_status","reason","seed_in_activation_program")

collect <- function() {
  filenames <- list.files(checkpoint,pattern="^worker_[0-9]+[.]tsv$",full.names=TRUE)
  rows <- list()
  for (name in filenames) {
    lines <- readLines(name,warn=FALSE)
    valid <- lengths(strsplit(lines,"\t",fixed=TRUE))==length(fields)
    if (any(valid)) rows[[length(rows)+1L]] <- read.delim(text=paste(lines[valid],collapse="\n"),header=FALSE,
                                                         col.names=fields,stringsAsFactors=FALSE,quote="",fill=FALSE)
  }
  if (length(rows)==0L) return(data.frame(job_id=integer(),job_status=character()))
  combined <- do.call(rbind,rows)
  completed <- combined[combined$job_status=="completed",,drop=FALSE]
  if (anyDuplicated(completed$job_id)) {
    spread_by_job <- split(completed$spread,completed$job_id)
    stopifnot(all(vapply(spread_by_job,function(s)length(unique(s))==1L,logical(1))))
  }
  combined <- combined[order(combined$job_status!="completed"),,drop=FALSE]
  combined[!duplicated(combined$job_id),,drop=FALSE]
}

if (!nrow(task_jobs)) {
  record(file.path(out,"METHOD_RUN.json"),list(execution_status="blocked_input",reason="both mapped programs empty"))
  quit(save="no",status=0L)
}
previous <- collect()
done <- previous$job_id[previous$job_status=="completed"]
todo <- task_jobs[!task_jobs$job_id %in% done,,drop=FALSE]
if (nrow(todo)) {
  # Dispatch larger seed degrees first for workload balance. Stored job RNG is independent of order.
  degree <- coverage$out_degree[match(todo$stable_id,coverage$stable_id)]
  todo <- todo[order(-degree,todo$job_id),,drop=FALSE]
  worker_count <- max(1L,min(workers,nrow(todo)))
  shards <- split(seq_len(nrow(todo)),(seq_len(nrow(todo))-1L) %% worker_count)
  result <- mclapply(shards,function(indexes) {
    dest <- file(file.path(checkpoint,sprintf("worker_%d.tsv",Sys.getpid())),open="at")
    on.exit(close(dest),add=TRUE)
    completed <- 0L
    for (index in indexes) {
      job <- todo[index,,drop=FALSE]
      result <- tryCatch(score_job(job),error=function(e)list(score=NA_real_,wall=NA_real_,error=conditionMessage(e)))
      success <- is.null(result$error)
      reason <- if(success) "none" else gsub("[\t\r\n]"," ",result$error)
      seed_counted <- any(program$target==target & program$stable_id==job$stable_id & program$direction=="activate")
      line <- c(job$job_id,job$target,job$stable_id,job$replicate,
                sprintf("%.17g",result$score),sprintf("%.6f",result$wall),
                if(success)"completed" else "failed_execution",reason,seed_counted)
      writeLines(paste(line,collapse="\t"),dest);flush(dest)
      completed <- completed+success
    }
    completed
  },mc.cores=worker_count,mc.set.seed=FALSE,mc.preschedule=TRUE)
  cat("Worker completion counts:",paste(result,collapse=" "), "\n")
}
complete <- collect()
stopifnot(all(complete$job_id %in% task_jobs$job_id))
joined <- merge(task_jobs,complete,by=c("job_id","target","stable_id","replicate"),sort=TRUE,all.x=TRUE)
write.table(joined,file.path(out,"NATIVE_SINGLETON_SPREAD.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
errors <- is.na(joined$job_status)|joined$job_status!="completed"
number <- sum(!errors)
record(file.path(out,"METHOD_RUN.json"),
       list(execution_status=if(any(errors))"failed_execution" else "completed",
            expected_jobs=nrow(task_jobs), completed_jobs=number,failed_or_missing_jobs=sum(errors),
            requested_workers=workers,native_inner_cores=1L,niter_per_job=1000L,independent_replicates=3L,
            config_sha256=cfg_hash,graph_nodes=vcount(graph),graph_edges=ecount(graph),
            native_function="IMCell::rank_by_spread",actual_rng=constant_seed_kind))
cat(target,number,"/",nrow(task_jobs),"native singleton jobs completed\n")
