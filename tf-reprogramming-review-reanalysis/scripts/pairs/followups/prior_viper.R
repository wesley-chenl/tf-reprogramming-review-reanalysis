args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args)==3L)
endpoint <- args[1L]; inp <- args[2L]; out <- args[3L]
stopifnot(as.character(utils::packageVersion("viper"))=="1.40.0")
options(digits=17)
tsv <- function(f) read.delim(file.path(inp,f), check.names=FALSE, colClasses="character", stringsAsFactors=FALSE, na.strings="")
pairs <- tsv("CLEAN_PAIRS.tsv")
expr <- tsv(paste0(endpoint,"_EFFECT.tsv"))
edges <- tsv("GRAPH.tsv")
stopifnot(nrow(edges)>0L, !anyDuplicated(expr$stable_id), !anyDuplicated(edges[,c("tf_gene_id","target_gene_id")]))
ids <- sort(unique(na.omit(sub("^GeneID:","",c(pairs$bhlh_id,pairs$partner_id)))),method="radix")
stopifnot(all(grepl("^[0-9]+$",ids)))
expr <- expr[is.finite(suppressWarnings(as.numeric(expr$d_g))),]
sig <- setNames(as.numeric(expr$d_g),expr$stable_id)
full_edges <- split(edges$target_gene_id,edges$tf_gene_id)
edges_measured <- edges[edges$tf_gene_id %in% ids & edges$target_gene_id %in% names(sig),]
reg <- lapply(split(edges_measured, edges_measured$tf_gene_id), function(e) {
    stopifnot(!anyDuplicated(e$target_gene_id), all(as.numeric(e$sign) %in% c(-1,1)))
    e <- e[order(e$target_gene_id,method="radix"),]
    list(tfmode=setNames(as.numeric(e$sign), e$target_gene_id), likelihood=setNames(rep(1,nrow(e)),e$target_gene_id))
})
pair_id <- function(a,b) paste(sort(c(a,b),method="radix"),collapse="--")
cov <- function(a,b) {
    ag <- names(reg[[a]]$tfmode); bg <- names(reg[[b]]$tfmode)
    common <- intersect(ag,bg)
    c(a=length(ag),b=length(bg),common=length(common),
      raw_common=length(intersect(full_edges[[a]],full_edges[[b]])),
      signed_agree=if(length(common))sum(reg[[a]]$tfmode[common]==reg[[b]]$tfmode[common]) else 0)
}
tab <- function(x, f) write.table(x,file.path(out,f),sep="\t",row.names=FALSE,quote=FALSE,na="")
mapping <- unique(pairs[,c("pair_id","bhlh_id","partner_id")])
cat(sprintf("endpoint=%s expression_genes=%d candidate_member_regs=%d native_viper=%s\n",endpoint,length(sig),length(reg),packageVersion("viper")))
run_case <- function(minsize, tag) {
    fitted <- viper::msviper(sig,reg,nullmodel=NULL,pleiotropy=FALSE,minsize=minsize,adaptive.size=FALSE,ges.filter=FALSE,synergy=0,cores=1,verbose=FALSE)
    regulators <- sort(names(fitted$es$nes),method="radix")
    stopifnot(length(regulators)>2L)
    ranked <- viper::msviperCombinatorial(fitted,regulators=regulators,nullmodel=NULL,minsize=minsize,adaptive.size=FALSE,level=2,cores=1,processAll=TRUE,verbose=FALSE)
    single <- fitted$es$nes
    es <- ranked$es$nes[grepl("--",names(ranked$es$nes),fixed=TRUE)]
    sizes <- ranked$es$size[names(es)]
    # Verify native pair names are exactly canonical IDs supplied to msviper.
    if(length(es)) stopifnot(all(names(es)==vapply(strsplit(names(es),"--",fixed=TRUE), function(x) pair_id(x[1],x[2]), "")))
    # Blind swap test: explicitly reverse the regulator list and repeat the official caller.
    swapped <- viper::msviperCombinatorial(fitted,regulators=rev(regulators),nullmodel=NULL,minsize=minsize,adaptive.size=FALSE,level=2,cores=1,processAll=TRUE,verbose=FALSE)
    if(length(es)) {
        inverse <- vapply(strsplit(names(es),"--",fixed=TRUE),function(x) paste(rev(x),collapse="--"),"")
        stopifnot(isTRUE(all.equal(unname(swapped$es$nes[inverse]),unname(es),check.attributes=FALSE,tolerance=1e-13)))
    }
    raw <- vector("list",nrow(pairs))
    positive <- vector("list",nrow(pairs))
    for(i in seq_len(nrow(pairs))) {
        p <- pairs[i,]
        a <- sub("^GeneID:","",p$bhlh_id); b <- sub("^GeneID:","",p$partner_id)
        valid <- !anyNA(c(a,b)); cv <- if(valid) cov(a,b) else rep(NA_integer_,5L)
        sz <- if(valid) cv[["common"]] else NA_integer_
        key <- if(valid) pair_id(a,b) else NA_character_
        v <- if(valid && key %in% names(es)) es[[key]] else NULL
        av <- if(valid) unname(single[a]) else NA_real_
        bv <- if(valid) unname(single[b]) else NA_real_
        reason <- if(!valid) "unresolved_identity" else if(min(cv[1:2])<minsize) "member_regulon_below_minsize" else if(sz<minsize) "common_regulon_below_minsize" else if(!length(v)) "native_pair_not_returned" else if(!is.finite(v)) "native_nonfinite_pair_nes" else ""
        if(!length(v) || !is.finite(v)) v <- NA_real_
        if(!is.na(v)) stopifnot(sz==unname(sizes[[key]]),sz>=minsize)
        details <- sprintf("minsize=%d;measured_a=%s;measured_b=%s;measured_common=%s;raw_common=%s;shared_sign_agree=%s;a_NES=%.17g;b_NES=%.17g",minsize,cv[1],cv[2],cv[3],cv[4],cv[5],av,bv)
        raw[[i]] <- data.frame(method=tag,endpoint=endpoint,variant=p$variant,pair_id=p$pair_id,score=v,reason=reason,details=details)
        is_positive <- is.finite(av) && is.finite(bv) && av>0 && bv>0
        keep <- is_positive && is.finite(v)
        positive[[i]] <- data.frame(method=if(tag=="msviper") "msviper_both_positive" else "msviper_minsize2_both_positive_diag",endpoint=endpoint,variant=p$variant,pair_id=p$pair_id,score=if(keep) v else NA_real_,reason=if(!is.finite(v))reason else if(!is_positive)"member_NES_not_both_positive" else "",details=details)
    }
    result <- rbind(do.call(rbind,raw),do.call(rbind,positive))
    cat(sprintf("minsize=%d member_regulons=%d native_all_pairs=%d variant_pair_scores=%d both_positive=%d\n",
                minsize,length(regulators),length(es),sum(is.finite(result$score)&result$method==tag),sum(is.finite(result$score)&result$method!=tag)))
    if(length(es)) {
        native <- data.frame(endpoint=endpoint,minsize=minsize,id_pair=names(es),native_pair_NES=unname(es),shared_size=unname(sizes))
        tab(native,paste0("viper_native_pairs_",endpoint,"_",minsize,".tsv"))
    }
    result
}
result <- rbind(run_case(5L,"msviper"),run_case(2L,"msviper_minsize2_diag"))
tab(result,paste0("viper_",endpoint,".tsv"))
writeLines(capture.output(sessionInfo()),file.path(out,paste0("viper_session_",endpoint,".txt")))
