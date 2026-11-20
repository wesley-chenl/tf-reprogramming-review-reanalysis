# IMCell weight control, using the same graph and target programs as the original frozen run.
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1L]]
output <- args[[2L]]
source_dir <- args[[3L]]
workers <- as.integer(args[[4L]])
mode <- args[[5L]]
stopifnot(workers >= 1L, workers <= 12L, mode %in% c("preflight", "initial", "extend", "replay"))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(parallel))
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
options(digits=17)
archive_functions <- new.env(parent=globalenv())
sys.source(file.path(source_dir, "IMCell.R"), envir=archive_functions)
sys.source(file.path(source_dir, "IMCell_setup.R"), envir=archive_functions)
sys.source(file.path(source_dir, "imcell_indexed_cpu.R"), envir=globalenv())

read_tsv <- function(p) {
  read.delim(p, colClasses="character", check.names=FALSE, quote="", comment.char="")
}
write_tsv <- function(x, p) {
  write.table(x, file=paste0(p, ".partial"), quote=FALSE, sep="\t", row.names=FALSE, na="")
  stopifnot(file.rename(paste0(p, ".partial"), p))
}
network <- read_tsv(file.path(input,"TRRUST_MOUSE_SIGNED.tsv"))
jobs <- read_tsv(file.path(output,"imcell_jobs.tsv"))
program <- read_tsv(file.path(output,"imcell_target_program.tsv"))
vertices <- sort(unique(c(network$tf_gene_id, network$target_gene_id)))
graph <- graph_from_data_frame(
  data.frame(from=network$tf_gene_id, to=network$target_gene_id, type=as.integer(network$sign), weight=1),
  directed=TRUE, vertices=data.frame(name=vertices)
)
stopifnot(ecount(graph)==3964L, vcount(graph)==1848L, sum(which_loop(graph))==16L,
          all(abs(E(graph)$type)==1L), nrow(jobs)==260L)
g_in <- archive_functions$weights_to_probability(graph, method="in_degree", multiplier=1, pmax=0.8)
g_mean <- archive_functions$weights_to_probability(graph, method="coin_flip", p=mean(E(g_in)$weight))
g_base <- archive_functions$weights_to_probability(graph, method="coin_flip", p=0.1)
in_strength <- unname(strength(graph, mode="in", loops=TRUE)[network$target_gene_id])
prob <- E(g_in)$weight
stopifnot(identical(prob, pmin(0.8, 1/in_strength)), all(E(g_mean)$weight == mean(prob)))
write_tsv(data.frame(from=network$tf_gene_id, to=network$target_gene_id,
  type=network$sign, target_in_strength=in_strength, loop=network$tf_gene_id==network$target_gene_id,
  probability_in_degree=prob, at_cap=prob==0.8, probability_matched_mean=E(g_mean)$weight,
  probability_p010=E(g_base)$weight), file.path(output, "IMCELL_EDGE_PROBABILITIES.tsv"))

graphs <- list(p010=g_base, in_degree=g_in, mean=g_mean)
endpoints <- unique(program$endpoint)
stopifnot(identical(endpoints, c("GSE53951","GSE114809")))
goals_for <- function(endpoint) {
  target <- program[program$endpoint==endpoint, ]
  list(up=intersect(target$stable_id[target$direction=="1"], vertices),
       down=intersect(target$stable_id[target$direction=="-1"], vertices))
}

# Preserve the archived MC aggregator and each original R draw. The indexed kernel
# replaces per-draw graph lookups, while tracking integer totals and second moments.
draw <- function(g, index, seed, goals, trials) {
  totals <- numeric(3L)
  native <- archive_functions$mc_ICP_parallel
  bindings <- new.env(parent=environment(native))
  bindings$random_ICP_parallel <- function(ig, seed, repressor_wins, targets_activate,
                                           targets_repress, return_spread) {
    ans <- r2_icp_indexed(index, seed, repressor_wins, targets_activate, targets_repress, return_spread)
    a <- ans$activated_spread
    b <- ans$repressed_spread
    totals <<- totals + c(a,b,(a+b)^2)
    ans
  }
  environment(native) <- bindings
  values <- native(g, sort(seed), repressor_wins=FALSE,
      targets_activate=goals$up, targets_repress=goals$down, niter=trials,
      return_spread=FALSE, num_cores=1L)
  stopifnot(identical(totals[1:2]/trials, unname(unlist(values))))
  totals
}
parity <- function(g, index, seeds, goals, token, n=64L) {
  set.seed(token)
  start <- .Random.seed
  means <- archive_functions$mc_ICP_parallel(g, sort(seeds), repressor_wins=FALSE,
      targets_activate=goals$up, targets_repress=goals$down,
      niter=n, return_spread=FALSE, num_cores=1L)
  native_rng <- .Random.seed
  assign(".Random.seed", start, envir=.GlobalEnv)
  fast <- draw(g,index,seeds,goals,n)
  same_rng <- identical(.Random.seed, native_rng)
  assign(".Random.seed", start, envir=.GlobalEnv)
  swapped <- draw(g,index,rev(seeds),goals,n)
  same_canonical_order <- identical(fast, swapped) && identical(.Random.seed, native_rng)
  stopifnot(identical(unname(unlist(means)), fast[1:2]/n), same_rng, same_canonical_order)
  data.frame(seeds=paste(seeds,collapse="|"), token, n,
    native_activation=means[[1]], native_repression=means[[2]],
    indexed_activation=fast[[1]]/n, indexed_repression=fast[[2]]/n,
    same_canonical_order, same_rng)
}

if (mode=="preflight" || mode=="replay") {
  # Two unequally probable edges reach one node with opposing signs. A third
  # edge continues with a negative sign. Keep an explicit self-loop.
  synth <- graph_from_data_frame(data.frame(from=c("a","b","x","x"),
      to=c("x","x","y","x"), type=c(1L,-1L,-1L,1L), weight=c(0.2,0.75,0.6,0.42)))
  synth_goals <- list(up="x",down="y")
  synth_prob <- c(0.2,0.75,0.6,0.42)
  s_index <- r2_icp_index(synth)
  stopifnot(identical(sort(s_index$weight), sort(synth_prob)), sum(which_loop(synth))==1L)
  all_off <- archive_functions$weights_to_probability(synth, method="coin_flip", p=0)
  synthetic <- rbind(
    cbind(scenario="conflicting_unequal_edges_and_loop",
          parity(synth,s_index,c("a","b"),synth_goals,975051L,128L)),
    cbind(scenario="all_edges_disabled",
          parity(all_off,r2_icp_index(all_off),c("a","b"),synth_goals,975052L,128L)),
    cbind(scenario="no_valid_targets",
          parity(synth,s_index,c("a","b"),list(up=character(),down=character()),975053L,128L))
  )
  stopifnot(all(synthetic$indexed_activation[-1]==0), all(synthetic$indexed_repression[-1]==0))
  write_tsv(synthetic, file.path(output,"IMCELL_SYNTHETIC_CHECK.tsv"))

  pair_subset <- head(jobs[order(jobs$audit_key),],12L)
  singleton_subset <- head(sort(unique(c(pair_subset$left,pair_subset$right))),4L)
  seed_lists <- c(lapply(seq_len(nrow(pair_subset)), function(i) c(pair_subset$left[[i]],pair_subset$right[[i]])),
                  as.list(singleton_subset))
  audit_keys <- c(pair_subset$audit_key, paste0("singleton:",singleton_subset))
  grid <- expand.grid(endpoint=endpoints, method=names(graphs), stringsAsFactors=FALSE)
  checks <- lapply(seq_len(nrow(grid)), function(block) {
    endpoint <- grid$endpoint[[block]]
    method <- grid$method[[block]]
    g <- graphs[[method]]
    index <- r2_icp_index(g)
    goals <- goals_for(endpoint)
    do.call(rbind,lapply(seq_along(seed_lists), function(i) {
      cbind(endpoint, method, audit_key=audit_keys[[i]],
            parity(g,index,seed_lists[[i]],goals,975100L+i))
    }))
  })
  write_tsv(do.call(rbind,checks),file.path(output,"IMCELL_NATIVE_PARITY.tsv"))
  message("Native parity passed for ",length(checks)*length(seed_lists)," checks")
}

seed_nodes <- sort(unique(c(jobs$left,jobs$right)))
seed_lists <- c(as.list(seed_nodes), lapply(seq_len(nrow(jobs)),function(i) c(jobs$left[[i]],jobs$right[[i]])))
seed_keys <- c(paste0("node:",seed_nodes),paste0("pair:",jobs$work_id))
stopifnot(length(seed_nodes)==36L, length(seed_keys)==296L)

build_rng <- function(block_seed) {
  set.seed(block_seed)
  state <- .Random.seed
  lapply(seq_along(seed_keys),function(i) {
    lapply(1:3,function(replicate) {
      state <<- nextRNGStream(state)
      state
    })
  })
}
# Preserve the old p=.10 substreams. The in-degree and matched-average graphs
# use disjoint, predetermined substreams; each pair's two singletons have their own.
block_seed <- function(endpoint, method) {
  endpoint_offset <- match(endpoint,endpoints)
  if (method=="p010") return(481991L + endpoint_offset)
  581991L + endpoint_offset + 2L*(match(method,c("in_degree","mean"))-1L)
}

checkpoint_name <- function(endpoint, method, precision) {
  file.path(output,paste0("mc_",endpoint,"_",method,"_",precision,".tsv"))
}
run_block <- function(endpoint, method, precision, subset=NULL) {
  path <- checkpoint_name(endpoint,method,precision)
  replaying <- !is.null(subset)
  if (file.exists(path) && !replaying) {
    cached <- read_tsv(path)
    stopifnot(nrow(cached)==3L*length(seed_keys), all(cached$endpoint==endpoint),
              all(cached$method==method), all(as.integer(cached$trial_count)==precision),
              identical(cached$job,rep(seed_keys,each=3L)),
              all(grepl("^[-0-9]+(,[-0-9]+){6}$",cached$rng_end)))
    message("Reusing ",basename(path))
    return(invisible(NULL))
  }
  g <- graphs[[method]]
  index <- r2_icp_index(g)
  goals <- goals_for(endpoint)
  seed <- block_seed(endpoint,method)
  streams <- build_rng(seed)
  if (precision==5000L && method!="p010" && !replaying) {
    prior <- read_tsv(checkpoint_name(endpoint,method,1000L))
    stopifnot(identical(prior$job, rep(seed_keys,each=3L)))
  } else {
    prior <- NULL
  }
  if (method=="p010") {
    legacy <- read_tsv(file.path(input,paste0("mc_block_",match(endpoint,endpoints),".tsv")))
    stopifnot(identical(legacy$job, seed_keys), all(legacy$endpoint==endpoint),
              all(legacy$method=="imcell_p010"))
  }

  run_job <- function(i) {
    rows <- lapply(1:3,function(replicate) {
      start <- streams[[i]][[replicate]]
      initial <- numeric(3)
      if (!is.null(prior)) {
        previous <- prior[(i-1L)*3L+replicate,]
        initial <- as.numeric(previous[c("activation_sum","repression_sum","sum_square_total")])
        continuation <- as.integer(strsplit(previous$rng_end,",",fixed=TRUE)[[1L]])
        assign(".Random.seed", continuation, envir=.GlobalEnv)
        current <- draw(g,index,seed_lists[[i]],goals,precision-1000L)
      } else if (method=="p010") {
        assign(".Random.seed", start, envir=.GlobalEnv)
        initial <- draw(g,index,seed_lists[[i]],goals,1000L)
        old_rep <- round(1000*as.numeric(legacy[[paste0("rep",replicate)]][[i]]))
        stopifnot(sum(initial[1:2])==old_rep)
        current <- if (precision > 1000L) draw(g,index,seed_lists[[i]],goals,precision-1000L) else numeric(3)
      } else {
        assign(".Random.seed", start, envir=.GlobalEnv)
        current <- draw(g,index,seed_lists[[i]],goals,precision)
      }
      combined <- initial+current
      data.frame(endpoint, method, job=seed_keys[[i]], replicate, trial_count=precision,
         activation_sum=combined[[1]], repression_sum=combined[[2]], sum_square_total=combined[[3]],
         initial_activation=initial[[1]], initial_repression=initial[[2]],
         initial_sum_square=initial[[3]], block_seed=seed,
         rng_start=paste(start,collapse=","), rng_end=paste(.Random.seed,collapse=","),
         n_active_goals=length(goals$up), n_repress_goals=length(goals$down))
    })
    result <- do.call(rbind,rows)
    if (method=="p010") {
      stopifnot(sum(result$initial_activation)==as.numeric(legacy$activation_sum[[i]]),
                sum(result$initial_repression)==as.numeric(legacy$repression_sum[[i]]),
                sum(result$initial_sum_square)==as.numeric(legacy$sum_square_total[[i]]))
    }
    result
  }
  ids <- if (replaying) subset else seq_along(seed_keys)
  started <- proc.time()[["elapsed"]]
  results <- mclapply(ids,run_job,mc.cores=workers,mc.preschedule=TRUE,mc.set.seed=FALSE)
  stopifnot(all(vapply(results,is.data.frame,logical(1L))))
  frame <- do.call(rbind,results)
  if (replaying) {
    cached <- read_tsv(path)
    cached <- cached[cached$job %in% frame$job,]
    cmp <- c("rng_start","rng_end","activation_sum","repression_sum","sum_square_total")
    stopifnot(identical(cached$job,frame$job), nrow(cached)==nrow(frame),
              all(vapply(cmp,function(col) identical(cached[[col]],as.character(frame[[col]])),logical(1L))))
    message("Full-precision deterministic substream replay passed for ", basename(path),
            " with ",length(ids)," seed jobs")
  } else {
    write_tsv(frame,path)
    message("Finished ",basename(path)," with ",nrow(frame)," replicate rows in ",
            round(proc.time()[["elapsed"]]-started,2)," seconds")
  }
}

if (mode=="initial") {
  for (endpoint in endpoints) {
    for (method in c("in_degree","mean")) run_block(endpoint, method, 1000L)
  }
} else if (mode %in% c("extend","replay")) {
  triggered <- read_tsv(file.path(output,"IMCELL_PRECISION_TRIGGER.tsv"))
  stopifnot(identical(triggered$endpoint,endpoints), all(triggered$trigger %in% c("0","1")))
  for (endpoint in endpoints) {
    active <- triggered$trigger[triggered$endpoint==endpoint]=="1"
    if (mode=="extend" && active) {
      for (method in names(graphs)) run_block(endpoint,method,5000L)
    }
    if (mode=="replay") {
      prec <- if (active) 5000L else 1000L
      audit_pair <- head(order(jobs$audit_key),12L)
      ids <- sort(c(seq_len(4L), length(seed_nodes)+audit_pair))
      for (method in if (active) names(graphs) else c("in_degree","mean")) {
        run_block(endpoint,method,prec,subset=ids)
      }
    }
  }
}
writeLines(c(R.version.string, paste("igraph",packageVersion("igraph")), paste("RNG",paste(RNGkind(),collapse=" / ")),
             paste("workers",workers), paste("user",Sys.info()[["user"]]), paste("mode",mode)),
           file.path(output,paste0("IMCELL_R_ENV_",mode,".txt")))

