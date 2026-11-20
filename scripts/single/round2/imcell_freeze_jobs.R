args <- commandArgs(trailingOnly=TRUE)
root <- args[[1]]
tasks <- c("liver", "cardiomyocyte", "betapanc")
program <- read.delim(file.path(root,"inputs_manifest/IMCELL_TARGET_PROGRAM.tsv"))
coverage <- read.delim(file.path(root,"inputs_manifest/IMCELL_SEED_COVERAGE.tsv"))
seeds <- sort(coverage$stable_id[coverage$out_degree>0], method="radix")
eligible_tasks <- tasks[tasks %in% program$target[tolower(program$mapped_to_graph)=="true"]]
jobs <- do.call(rbind,lapply(eligible_tasks,function(target) {
  data.frame(target=target,stable_id=rep(seeds,each=3L),replicate=rep(1:3,length(seeds)))
}))
jobs$job_id <- seq_len(nrow(jobs))
RNGkind("L'Ecuyer-CMRG")
set.seed(20261118L)
streams <- matrix(0L,nrow=nrow(jobs),ncol=7L)
for(index in seq_len(nrow(jobs))) {
  streams[index,] <- .Random.seed
  .Random.seed <- parallel::nextRNGStream(.Random.seed)
}
colnames(streams) <- paste0("rng_",seq_len(7L))
stopifnot(!anyDuplicated(as.data.frame(streams)))
write.table(cbind(jobs,streams),file.path(root,"inputs_manifest/IMCELL_RNG_STREAMS.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
cat("Bound RNG streams to",nrow(jobs),"named target/TF/replicate jobs\n")
