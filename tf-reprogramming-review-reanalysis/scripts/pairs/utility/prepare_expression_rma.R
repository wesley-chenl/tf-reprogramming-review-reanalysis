args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop("Usage: prepare_expression_rma.R WORKSPACE_ROOT")
root <- normalizePath(args[1L], mustWork=TRUE)
raw <- file.path(root, "expression_raw")
out <- file.path(root, "validation")
accessions <- c("GSM1303947", "GSM1303949", "GSM1303943", "GSM1303945")
files <- file.path(raw, "GSE53951_cels", paste0(accessions, ".CEL"))
stopifnot(all(file.exists(files)), length(unique(unname(tools::md5sum(files)))) == 4L)

required <- c("affy", "mouse4302cdf", "mouse4302.db", "AnnotationDbi", "org.Mm.eg.db", "affyio", "preprocessCore")
for (pkg in required) stopifnot(requireNamespace(pkg, quietly=TRUE))
arrays <- affy::ReadAffy(filenames=files, sampleNames=accessions)
stopifnot(identical(Biobase::sampleNames(arrays), accessions))
stopifnot(identical(tolower(gsub("_", "", affy::cdfName(arrays))), "mouse4302"))
needs <- environmentName(affy::getCdfInfo(arrays))
stopifnot(length(needs) == 1L)
# One single joint RMA run for these four arrays. No treatment controls enter normalization.
fit <- affy::rma(arrays, background=TRUE, normalize=TRUE, verbose=TRUE)
probes <- Biobase::exprs(fit)
stopifnot(nrow(probes) == 45101L, ncol(probes) == 4L, identical(colnames(probes), accessions),
          all(is.finite(probes)), !anyDuplicated(rownames(probes)))

db <- mouse4302.db::mouse4302.db
mapping <- AnnotationDbi::select(db, keys=rownames(probes), keytype="PROBEID", columns="ENTREZID")
ok <- !is.na(mapping$ENTREZID) & nzchar(mapping$ENTREZID)
distinct <- split(as.character(mapping$ENTREZID[ok]), mapping$PROBEID[ok])
by_probe <- vapply(distinct, function(v) paste(sort(unique(v)), collapse="|"), character(1))
selected <- unname(by_probe[rownames(probes)])
selected[is.na(selected)] <- ""
write.table(
  data.frame(probeset_id=rownames(probes), bioc_entrez_ids=selected, probes, check.names=FALSE),
  file=file.path(raw, "GPL1261_RMA_PROBES.tsv"),
  quote=FALSE, row.names=FALSE, sep="\t", na=""
)
write.table(
  data.frame(package=required, version=vapply(required, function(p) as.character(packageVersion(p)), character(1)),
             db_file=vapply(required, function(p) {
               if (p %in% c("mouse4302.db", "org.Mm.eg.db")) {
                 package <- getExportedValue(p, p)
                 basename(AnnotationDbi::dbfile(package))
               } else ""
             }, character(1))),
  file=file.path(out, "EXPRESSION_R_PACKAGES.tsv"), quote=FALSE, row.names=FALSE, sep="\t", na=""
)
write.table(AnnotationDbi::metadata(db), file=file.path(out, "GPL1261_ANNOTATION_METADATA.tsv"), quote=FALSE, row.names=FALSE, sep="\t")
write.table(
  data.frame(gsm=accessions, raw_CELs=basename(files), platform=rep("GPL1261",4L),
             rma_min=apply(probes,2,min), rma_median=apply(probes,2,median), rma_max=apply(probes,2,max),
             rma_units=rep("log2-scale joint RMA",4L)),
  file=file.path(out, "GPL1261_RMA_ARRAY_QC.tsv"), quote=FALSE, row.names=FALSE, sep="\t", na=""
)
capture.output(sessionInfo(), file=file.path(out,"EXPRESSION_R_SESSION.txt"))
cat("Joint GPL1261 RMA complete:", nrow(probes), "probe sets x", ncol(probes), "arrays\n")
