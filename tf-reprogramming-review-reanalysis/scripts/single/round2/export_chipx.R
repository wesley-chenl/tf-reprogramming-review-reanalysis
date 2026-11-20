args <- commandArgs(trailingOnly=TRUE)
data_env <- new.env(parent=emptyenv())
objects <- load(args[[1]], envir=data_env)
stopifnot(identical(objects, "chipX_gs"))
edges <- data_env$chipX_gs
stopifnot(is.data.frame(edges), identical(names(edges), c("TF", "TG")),
          all(vapply(edges, is.character, logical(1))))
write.table(edges, args[[2]], sep="\t", row.names=FALSE, quote=FALSE, na="NA")
cat("archive chipX_gs:", nrow(edges), "rows,", length(unique(c(edges$TF, edges$TG))), "symbols; only TF and TG fields\n")
