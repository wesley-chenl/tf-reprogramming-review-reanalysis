# Interface-only repair of a one-column dimension drop in installed VIPER.
# Statistical background, activity algorithm and all native defaults retained.
viper_dimension_safe<-viper::viper
txt<-deparse(body(viper_dimension_safe),width.cutoff=500)
needle<-'tt <- tt[rownames(tt) %in% unique(tmp), ]'
stopifnot(sum(grepl(needle,txt,fixed=TRUE))==1)
txt<-gsub(needle,'tt <- tt[rownames(tt) %in% unique(tmp), , drop = FALSE]',txt,fixed=TRUE)
body(viper_dimension_safe)<-parse(text=txt)[[1]]
environment(viper_dimension_safe)<-asNamespace('viper')
