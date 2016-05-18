cnv.data.chr <-
function(chr, position, normal.ref.rd, normal.alt.rd, tumor.ref.rd, tumor.alt.rd,
                         normal.total.coverage, tumor.total.coverage)
{
	## require(DNAcopy)
	## all variables should have the same length
	## position has to be sorted
	## recommend this step done chr by chr
	normal.coverage <- normal.ref.rd + normal.alt.rd
	tumor.coverage  <- tumor.ref.rd + tumor.alt.rd
	
	## Log2 Ratio
    log2ratio <- log2(tumor.coverage/normal.coverage) - log2(tumor.total.coverage/normal.total.coverage)
    ## smoothing to remove outlier using a toolkit in DNAcopy
    n <- length(log2ratio)
    log2ratio <- smooth.CNA( CNA(genomdat=log2ratio, chrom=rep(1,n), maploc=1:n) )[[3]]
    log2ratio[abs(log2ratio)==Inf] <- NA

    ## normal BAF and mBAF
    normal.baf <- normal.alt.rd/normal.coverage
    normal.mbaf <- abs(normal.baf - 0.5) + 0.5
    normal.mbaf[normal.mbaf > 0.99] <- NA
    n <- length(normal.mbaf)
    normal.mbaf <- smooth.CNA( CNA(genomdat=normal.mbaf, chrom=rep(1,n), maploc=1:n) )[[3]]
    normal.mbaf[abs(normal.mbaf)==Inf] <- NA

    ## tumor BAF and mBAF
    tumor.baf <- tumor.alt.rd/tumor.coverage
    tumor.mbaf <- abs(tumor.baf - 0.5) + 0.5
    tumor.mbaf[tumor.mbaf > 0.99] <- NA
    n <- length(tumor.mbaf)
    tumor.mbaf <- smooth.CNA( CNA(genomdat=tumor.mbaf, chrom=rep(1,n), maploc=1:n) )[[3]]
    tumor.mbaf[abs(tumor.mbaf)==Inf] <- NA

    ## tumor vs normal
    pair.mbaf <- log2(tumor.mbaf/normal.mbaf)
    n <- length(pair.mbaf)
    pair.mbaf <- smooth.CNA( CNA(genomdat=pair.mbaf, chrom=rep(1,n), maploc=1:n) )[[3]]
    pair.mbaf[abs(pair.mbaf)==Inf] <- NA
    
    ## remove probes with missing value
    idx.na <- which(is.na(log2ratio) | is.na(normal.mbaf) | is.na(tumor.mbaf) | is.na(pair.mbaf))
    if (length(idx.na)>=1) {
        chr         <- chr[-idx.na]
        position    <- position[-idx.na]
        log2ratio   <- log2ratio[-idx.na]
        normal.baf  <- normal.baf[-idx.na]
        normal.mbaf <- normal.mbaf[-idx.na]
        tumor.baf   <- tumor.baf[-idx.na]
        tumor.mbaf  <- tumor.mbaf[-idx.na]
        pair.mbaf   <- pair.mbaf[-idx.na]
    }
    
    ## impute remaining missing data    
    if (sum(is.na(log2ratio))>0) log2ratio <- impute.missing.data(log2ratio)
    if (sum(is.na(normal.mbaf))>0) normal.mbaf <- impute.missing.data(normal.mbaf)
    if (sum(is.na(tumor.mbaf))>0) tumor.mbaf <- impute.missing.data(tumor.mbaf)    
    if (sum(is.na(pair.mbaf))>0) pair.mbaf <- impute.missing.data(pair.mbaf)

    return(data.frame(chr=chr, position=position, log2ratio=log2ratio, log2mBAF=pair.mbaf,
                      normal.BAF=normal.baf, normal.mBAF=normal.mbaf,
                      tumor.BAF=tumor.baf, tumor.mBAF=tumor.mbaf,
                      stringsAsFactors=FALSE))
}
