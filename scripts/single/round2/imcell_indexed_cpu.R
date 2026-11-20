# Execution compatibility for the archived singleton call on a full signed graph.
# The same ordered edges, R RNG calls, polarity, and returns are used. Integer
# endpoints replace per-layer data-frame creation and repeated igraph lookups.

r2_icp_index <- function(ig) {
  nv <- igraph::vcount(ig)
  vertices <- igraph::V(ig)
  node_names <- vertices$name
  adjacent <- lapply(igraph::adjacent_vertices(ig,vertices,mode="out"),as.integer)
  counts <- lengths(adjacent)
  origin <- rep.int(seq_len(nv),counts)
  target <- unlist(adjacent,use.names=FALSE)
  pairs <- as.vector(rbind(origin,target))
  native_edges <- igraph::E(ig,P=pairs)
  list(adjacency=split(seq_along(target),factor(origin,levels=seq_len(nv))), origin=origin,target=target,
       type=native_edges$type,weight=native_edges$weight,
       vertex_names=node_names,node_lookup=setNames(seq_len(nv),node_names))
}

r2_icp_indexed <- function(index,seed,repressor_wins,targets_activate,targets_repress,return_spread) {
  positive <- unname(index$node_lookup[seed])
  negative <- integer()
  new <- positive
  state <- integer(length(index$vertex_names))
  state[new] <- 1L
  while (length(new)>0L) {
    # Preserve native origin order, each origin's igraph adjacency order, and the number of RNG draws.
    proposal <- unlist(index$adjacency[new],use.names=FALSE)
    proposal <- proposal[state[index$target[proposal]]==0L]
    if (!length(proposal)) break
    flips <- stats::runif(length(proposal),min=0,max=1)
    live <- proposal[flips < index$weight[proposal]]
    if (!length(live)) break
    polarity <- state[index$origin[live]] * index$type[live]
    dest <- index$target[live]
    if (!repressor_wins) {
      native_random_order <- sample(seq_len(length(live)))
      dest <- dest[native_random_order]
      polarity <- polarity[native_random_order]
      first_live <- !duplicated(dest)
      next_positive <- dest[first_live & polarity==1L]
      next_negative <- dest[first_live & polarity==(-1L)]
    } else {
      next_negative <- dest[polarity==(-1L)]
      next_positive <- setdiff(dest[polarity==1L],next_negative)
    }
    positive <- c(positive,next_positive)
    negative <- c(negative,next_negative)
    state[next_positive] <- 1L
    state[next_negative] <- -1L
    new <- union(next_positive,next_negative)
  }
  positive <- index$vertex_names[unique(positive)]
  negative <- if (length(negative)) index$vertex_names[unique(negative)] else NULL
  if (!is.null(targets_activate) && is.null(targets_repress)) {
    positive <- intersect(positive,targets_activate)
    negative <- NULL
  } else if (is.null(targets_activate) && !is.null(targets_repress)) {
    positive <- NULL
    negative <- intersect(negative,targets_repress)
  } else if (!is.null(targets_activate) && !is.null(targets_repress)) {
    positive <- intersect(positive,targets_activate)
    negative <- intersect(negative,targets_repress)
  }
  result <- list(activated_spread=length(positive),repressed_spread=length(negative))
  if (return_spread) result <- c(result,list(activated=positive,repressed=negative))
  result
}

r2_install_indexed_mc <- function() {
  namespace <- asNamespace("IMCell")
  native_mc <- get("mc_ICP_parallel",envir=namespace)
  func <- native_mc
  statements <- as.list(body(func))
  printed <- vapply(statements,function(statement)paste(deparse(statement),collapse=" "),character(1L))
  location <- which(startsWith(printed,"rcs <- mclapply(") & grepl("random_ICP_parallel",printed,fixed=TRUE))
  stopifnot(length(location)==1L)
  # Keep the archive's full-index graph preparation, MC summary, target lists, and public rank_by_spread unchanged.
  statements[[location]] <- quote(rcs <- {
    r2_ordered_index <- r2_icp_index(ig)
    mclapply(1:niter,function(x)r2_icp_indexed(r2_ordered_index,seed,repressor_wins,
                                             targets_activate,targets_repress,return_spread),mc.cores=num_cores)
  })
  body(func) <- as.call(statements)
  bindings <- new.env(parent=environment(native_mc))
  bindings$r2_icp_index <- r2_icp_index
  bindings$r2_icp_indexed <- r2_icp_indexed
  environment(func) <- bindings
  unlockBinding("mc_ICP_parallel",namespace)
  assign("mc_ICP_parallel",func,envir=namespace)
  lockBinding("mc_ICP_parallel",namespace)
  invisible(native_mc)
}

r2_restore_native_mc <- function(native) {
  namespace <- asNamespace("IMCell")
  unlockBinding("mc_ICP_parallel",namespace)
  assign("mc_ICP_parallel",native,envir=namespace)
  lockBinding("mc_ICP_parallel",namespace)
}
