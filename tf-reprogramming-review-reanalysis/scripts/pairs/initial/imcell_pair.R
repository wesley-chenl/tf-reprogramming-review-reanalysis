args <- commandArgs(trailingOnly=TRUE)
input <- args[[1L]]
output <- args[[2L]]
source_dir <- args[[3L]]
workers <- as.integer(args[[4L]])
mode <- args[[5L]]
stopifnot(workers >= 1L, workers <= 12L, mode %in% c("score", "replay", "preflight"))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(parallel))
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
options(digits=17)
archive_functions <- new.env(parent=globalenv())
sys.source(file.path(source_dir,"IMCell.R"), envir=archive_functions)
sys.source(file.path(source_dir,"IMCell_setup.R"), envir=archive_functions)
sys.source(file.path(source_dir,"imcell_indexed_cpu.R"), envir=globalenv())
read_tsv <- function(p) read.delim(p, colClasses="character", check.names=FALSE, quote="", comment.char="")
write_tsv <- function(x,p) write.table(x,file=p,quote=FALSE,sep="\t",row.names=FALSE,na="")
network <- read_tsv(file.path(input,"TRRUST_MOUSE_SIGNED.tsv"))
jobs <- read_tsv(file.path(output,"imcell_jobs.tsv"))
program <- read_tsv(file.path(output,"imcell_target_program.tsv"))
audit_jobs <- head(jobs[order(jobs$audit_key),],12L)
stopifnot(nrow(audit_jobs)==12L)
probabilities <- c(imcell_p010=0.10,imcell_p005=0.05,imcell_p020=0.20)
endpoints <- unique(program$endpoint)
block_grid <- expand.grid(endpoint=endpoints,method=names(probabilities),stringsAsFactors=FALSE)
vertices <- sort(unique(c(network$tf_gene_id, network$target_gene_id)))
graph <- graph_from_data_frame(data.frame(from=network$tf_gene_id,to=network$target_gene_id,
  type=as.integer(network$sign),weight=1),directed=TRUE,vertices=data.frame(name=vertices))
stopifnot(ecount(graph)==3964L, sum(which_loop(graph))==16L, all(abs(E(graph)$type)==1L))

draw <- function(g, index, seed, goals, trials) {
  totals <- numeric(3L)
  native <- archive_functions$mc_ICP_parallel
  bindings <- new.env(parent=environment(native))
  bindings$random_ICP_parallel <- function(ig, seed, repressor_wins, targets_activate,
                                           targets_repress, return_spread) {
    ans <- r2_icp_indexed(index,seed,repressor_wins,targets_activate,targets_repress,return_spread)
    a <- ans$activated_spread
    b <- ans$repressed_spread
    totals <<- totals+c(a,b,(a+b)^2)
    ans
  }
  environment(native) <- bindings
  # The archived MC entry point summarizes calls to the verified ordered edge index.
  means <- native(g,sort(seed),repressor_wins=FALSE,targets_activate=goals$up,
                  targets_repress=goals$down,niter=trials,return_spread=FALSE,num_cores=1L)
  stopifnot(identical(totals[1:2]/trials,unname(unlist(means))))
  totals
}

record_parity <- function(graph, index, goals, endpoint, method) {
  results <- lapply(seq_len(nrow(audit_jobs)),function(i) {
    seeds <- sort(c(audit_jobs$left[[i]],audit_jobs$right[[i]]))
    trial_n <- 64L
    token <- 975100L+i
    set.seed(token)
    start <- .Random.seed
    native <- archive_functions$mc_ICP_parallel(graph, seeds, repressor_wins=FALSE,
       targets_activate=goals$up, targets_repress=goals$down, niter=trial_n, return_spread=FALSE, num_cores=1L)
    native_rng <- .Random.seed
    assign(".Random.seed",start,envir=.GlobalEnv)
    accelerated <- draw(graph,index,seeds,goals,trial_n)
    fast_rng <- .Random.seed
    assign(".Random.seed",start,envir=.GlobalEnv)
    reverse <- draw(graph,index,rev(seeds),goals,trial_n)
    reverse_rng <- .Random.seed
    assign(".Random.seed",start,envir=.GlobalEnv)
    raw_reverse <- archive_functions$mc_ICP_parallel(graph, rev(seeds), repressor_wins=FALSE,
       targets_activate=goals$up, targets_repress=goals$down, niter=trial_n, return_spread=FALSE, num_cores=1L)
    matched <- identical(unname(accelerated[1:2])/trial_n,unname(unlist(native))) &&
               identical(native_rng,fast_rng)
    canonical <- identical(reverse,accelerated) && identical(fast_rng,reverse_rng)
    stopifnot(matched, canonical)
    data.frame(endpoint,method,work_id=audit_jobs$work_id[[i]],hash=audit_jobs$audit_key[[i]],
      n=trial_n,native_total=sum(unlist(native)),indexed_total=sum(accelerated[1:2])/trial_n,
      native_raw_reverse_total=sum(unlist(raw_reverse)),same_after_canonical_order=canonical,
      same_native_means_and_rng=matched)
  })
  do.call(rbind,results)
}

record_greedy <- function(graph, goals, endpoint, method) {
  raw <- graph
  edge_probability <- probabilities[[method]]
  edge_attr(raw,"weight") <- rep(1,ecount(raw))
  candidates <- sort(unique(c(head(audit_jobs$left,2L),head(audit_jobs$right,2L))))
  set.seed(975312L)
  capture.output(result <- archive_functions$IMCell(raw,kmax=2L,repressor_wins=FALSE,tfs=candidates,
     targets_activate=goals$up,targets_repress=goals$down,niter=64L,
     edge_probability_method="coin_flip",edge_probability_multiplier=2*edge_probability,
     min_marginal_spread=0,return_spread=FALSE,num_cores=1L))
  enumerated <- combn(candidates,2,simplify=FALSE)
  indexed <- r2_icp_index(archive_functions$weights_to_probability(raw,method="coin_flip",p=edge_probability))
  pair_names <- vapply(enumerated,function(seed) paste(sort(seed),collapse="|"),character(1))
  spreads <- vapply(seq_along(enumerated),function(i) {
    set.seed(975412L+i)
    sum(draw(raw,indexed,enumerated[[i]],goals,64L)[1:2])/64
  },numeric(1))
  selections <- paste(sort(result$solution_set),collapse="|")
  data.frame(endpoint,method,allowed_seed_ids=paste(candidates,collapse="|"),
    greedy=selections,greedy_in_this_benchmark=selections %in% jobs$work_id,
    exhaustive_winner=pair_names[[which.max(spreads)]],pairs_compared=length(enumerated),
    exhaustive_spreads=paste(paste(pair_names,spreads,sep="="),collapse=";"),n=64L)
}

synthetic <- function() {
  g <- graph_from_data_frame(data.frame(from=c("a","b","x","b"),to=c("x","x","y","b"),
                                        type=c(1L,1L,-1L,-1L),weight=1),directed=TRUE)
  goals <- list(up="x",down="y")
  evidence <- lapply(c(0,1),function(p) {
    current <- archive_functions$weights_to_probability(g,method="coin_flip",p=p,pmax=1)
    ix <- r2_icp_index(current)
    results <- lapply(list(c("a","b"),c("b","a"),"a","b"),function(seeds) {
      set.seed(972051L)
      draw(current,ix,seeds,goals,16L)[1:2]/16
    })
    stopifnot(identical(results[[1]],results[[2]]),identical(results[[1]],c(p,p)))
    data.frame(p,joint_activation=results[[1]][1],joint_repression=results[[1]][2],
      singleton_sum=sum(results[[3]])+sum(results[[4]]),swap_exact=TRUE,
      test="two seeds converge on same target, followed by repressive continuation and untouched self-loop")
  })
  do.call(rbind,evidence)
}

block <- function(index) {
  endpoint <- block_grid$endpoint[[index]]
  method <- block_grid$method[[index]]
  prob <- probabilities[[method]]
  g <- archive_functions$weights_to_probability(graph,method="coin_flip",p=prob)
  stopifnot(all(E(g)$weight==prob))
  indexed <- r2_icp_index(g)
  target <- program[program$endpoint==endpoint,]
  goals <- list(up=intersect(target$stable_id[target$direction=="1"],vertices),
                down=intersect(target$stable_id[target$direction=="-1"],vertices))
  parity <- NULL
  greedy <- NULL
  if (mode=="score" || mode=="preflight") {
    parity <- record_parity(g,indexed,goals,endpoint,method)
    if (method=="imcell_p010") greedy <- record_greedy(g,goals,endpoint,method)
  }
  if (mode=="preflight") return(list(parity=parity,greedy=greedy))
  output_path <- file.path(output,paste0("mc_block_",index,".tsv"))
  if (file.exists(output_path)) {
    receipt <- read_tsv(output_path)
    stopifnot(nrow(receipt)==nrow(jobs)+length(unique(c(jobs$left,jobs$right))),
              all(receipt$endpoint==endpoint),all(receipt$method==method))
    message("Verified cached block ",endpoint," ",method)
    return(list(parity=parity,greedy=greedy))
  }
  nodes <- sort(unique(c(jobs$left,jobs$right)))
  # Streams are allocated before forking. Ordering and worker scheduling do not alter them.
  task_seeds <- c(paste0("node:",nodes),paste0("pair:",jobs$work_id))
  set.seed(481991L+index)
  states <- vector("list", length(task_seeds))
  current <- .Random.seed
  for (i in seq_along(states)) {
    states[[i]] <- vector("list",3L)
    for (rep in 1:3) {
      current <- nextRNGStream(current)
      states[[i]][[rep]] <- current
    }
  }
  count_trials <- 1000L
  one <- function(i) {
    if (i<=length(nodes)) {
      input_seed <- nodes[[i]]
    } else {
      j <- i-length(nodes)
      input_seed <- c(jobs$left[[j]],jobs$right[[j]])
    }
    tmp <- vapply(1:3,function(rep) {
      assign(".Random.seed",states[[i]][[rep]],envir=.GlobalEnv)
      draw(g,indexed,input_seed,goals,count_trials)
    },numeric(3L))
    s <- rowSums(tmp)
    n <- count_trials*3
    spread <- s[[1]]+s[[2]]
    mcse <- sqrt(max(0,(s[[3]]-spread^2/n)/(n-1))/n)
    data.frame(endpoint,method,job=task_seeds[[i]],activation_sum=s[[1]],repression_sum=s[[2]],
      score_num=spread,score_den=n,sum_square_total=s[[3]],mc_se=mcse,rep1=sum(tmp[1:2,1])/count_trials,
      rep2=sum(tmp[1:2,2])/count_trials,rep3=sum(tmp[1:2,3])/count_trials,
      n_active_goals=length(goals$up),n_repress_goals=length(goals$down))
  }
  # Blocks are sequential to avoid nested process grids.
  results <- mclapply(seq_along(task_seeds),one,mc.cores=workers,mc.preschedule=TRUE,mc.set.seed=FALSE)
  stopifnot(all(vapply(results,is.data.frame,logical(1))))
  result <- do.call(rbind,results)
  write_tsv(result,paste0(output_path,".partial"))
  stopifnot(file.rename(paste0(output_path,".partial"),output_path))
  message(endpoint," ",method," completed ",nrow(result)," seed jobs at ",3L*count_trials," draws each")
  list(parity=parity,greedy=greedy)
}

writeLines(c(R.version.string,paste("igraph",as.character(packageVersion("igraph"))),
  paste("RNG",paste(RNGkind(),collapse=" / ")),paste("worker_count",workers),paste("mode",mode),
  paste("non-root UID",Sys.info()[["user"]])), file.path(output,"imcell_R_ENV.txt"))
checks <- lapply(seq_len(nrow(block_grid)),block)
if (mode!="replay") {
  write_tsv(synthetic(),file.path(output,"IMCELL_SYNTHETIC_CHECK.tsv"))
  write_tsv(do.call(rbind,lapply(checks,"[[","parity")),file.path(output,"IMCELL_NATIVE_PARITY.tsv"))
  write_tsv(do.call(rbind,lapply(checks,"[[","greedy")),file.path(output,"IMCELL_GREEDY_CHECK.tsv"))
}
