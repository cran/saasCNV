snp.cnv.data <-
function(snp, min.chr.probe=100, verbose=FALSE) 
{
	## require(DNAcopy)
	
	## initialize
	n.Normal.GT.na <- sum(is.na(snp$Normal.GT))
	if ( n.Normal.GT.na > 0.2 * nrow(snp) ) {
		normal.genotype <- FALSE
		tmp <- BAF2mBAF(x=snp$Tumor.BAF[snp$CHROM %in% paste0("chr", 1:22)], gt=NULL, het.gt="0/1",
                        mBAF.thd=0.97, win.thd=0.6, w=1, k=2, median.adjust=FALSE, impute=FALSE)$mBAF
        normal.mBAF.base <- median(tmp[tmp < 0.55], na.rm=TRUE)
        flag <- ifelse(sum(!is.na(tmp)) < 0.1 * length(tmp), 1, 0)
	} else {
		normal.genotype <- TRUE
		flag <- 0
	}	
	
    ## count #probes in each chr
    chr.probe.count <- table(snp$CHROM)
    chrs <- names(chr.probe.count)[chr.probe.count >= min.chr.probe]
    chrs <- sub("^chr","",chrs)
    all.chrs <- c(1:22,"X","Y")
    chrs <- all.chrs[all.chrs %in% chrs] ## this is to make the order of chr as what we usually see 1-22, X, Y

    data <- NULL
    for (chr in chrs) {
	    ## data processing
	    snp.chr <- snp[snp$CHROM==paste0("chr",chr), ]
	    snp.chr <- snp.chr[order(snp.chr$POS),]
	    n.chr   <- nrow(snp.chr)
	    
	    ## Log2 Ratio
        log2ratio <- smooth.CNA( CNA(genomdat=snp.chr$Tumor.LRR, chrom=rep(1, n.chr), maploc=1:n.chr, data.type="logratio") )[[3]]
        if (sum(is.na(log2ratio)) > 0) log2ratio <- impute.missing.data(log2ratio)
 
        if (normal.genotype) {
           ## normal mBAF
           tmp <- BAF2mBAF(x=snp.chr$Normal.BAF, gt=snp.chr$Normal.GT, het.gt="0/1",
                           mBAF.thd=0.97, win.thd=0.8, w=1, k=2, median.adjust=FALSE, impute=FALSE)
           normal.mBAF <- smooth.CNA( CNA(genomdat=tmp$mBAF, chrom=rep(1, length(tmp$mBAF)), maploc=1:length(tmp$mBAF)) )[[3]]
           idx.normal.mBAF <- tmp$idx
           
           ## tumor mBAF
           tmp <- BAF2mBAF(x=snp.chr$Tumor.BAF, gt=snp.chr$Normal.GT, het.gt="0/1",
                           mBAF.thd=0.97, win.thd=0.8, w=1, k=2, median.adjust=FALSE, impute=FALSE)
           tumor.mBAF <- smooth.CNA( CNA(genomdat=tmp$mBAF, chrom=rep(1, length(tmp$mBAF)), maploc=1:length(tmp$mBAF)) )[[3]]
           idx.tumor.mBAF <- tmp$idx
           
           ## log2mBAF
           tmp <- log2(tumor.mBAF/normal.mBAF)
           log2mBAF <- smooth.CNA( CNA(genomdat=tmp, chrom=rep(1, length(tmp)), maploc=1:length(tmp)) )[[3]]
           
           use.in.seg <- rep(0, length(log2mBAF))
           use.in.seg[!is.na(log2mBAF)] <- 1
        }
 
        if (!normal.genotype) {           
           ## tumor mBAF
           tmp <- BAF2mBAF(x=snp.chr$Tumor.BAF, gt=NULL, het.gt="0/1",
                           mBAF.thd=0.97, win.thd=0.6, w=1, k=2, median.adjust=FALSE, impute=FALSE)
           tumor.mBAF <- smooth.CNA( CNA(genomdat=tmp$mBAF, chrom=rep(1, length(tmp$mBAF)), maploc=1:length(tmp$mBAF)) )[[3]]

           ## normal mBAF
           snp.chr$Normal.BAF  <- 0.5
           normal.mBAF <- normal.mBAF.base
           
           ## log2mBAF
           tmp <- log2(tumor.mBAF/normal.mBAF.base)
           log2mBAF <- smooth.CNA( CNA(genomdat=tmp, chrom=rep(1, length(tmp)), maploc=1:length(tmp)) )[[3]]
           
           use.in.seg <- rep(0, length(log2mBAF))
           use.in.seg[!is.na(log2mBAF)] <- 1
        }

	    data.chr <- data.frame(chr=snp.chr$CHROM, position=snp.chr$POS,
	                           use.in.seg=use.in.seg, flag=flag,
	                           log2ratio=log2ratio, log2mBAF=log2mBAF,
                               normal.BAF=snp.chr$Normal.BAF, normal.mBAF=normal.mBAF,
                               tumor.BAF=snp.chr$Tumor.BAF, tumor.mBAF=tumor.mBAF,
                               stringsAsFactors=FALSE)
    
        ## collect information
        data <- rbind(data, data.chr)
    
        if (verbose) cat("finish data transformation of chr =",chr,"\n")
    } ## for chr

    return(data) 
}
