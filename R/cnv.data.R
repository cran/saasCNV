cnv.data <-
function(vcf, min.chr.probe=100, verbose=FALSE) 
{
	## initialize
	idx <- which(vcf$Normal.GT=="0/1")
    vcf <- vcf[idx,]
    normal.total.coverage <- sum(vcf$Normal.REF.DP + vcf$Normal.ALT.DP, na.rm=TRUE)
    tumor.total.coverage  <- sum(vcf$Tumor.REF.DP + vcf$Tumor.ALT.DP, na.rm=TRUE)

    ## count #probes in each chr
    chr.probe.count <- table(vcf$CHROM)
    chrs <- names(chr.probe.count)[chr.probe.count >= min.chr.probe]
    chrs <- sub("^chr","",chrs)
    all.chrs <- c(1:22,"X","Y")
    chrs <- all.chrs[all.chrs %in% chrs] ## this is to make the order of chr as what we usually see 1-22, X, Y

    data <- NULL
    for (chr in chrs) {
	    ## data processing
	    vcf.chr <- vcf[vcf$CHROM==paste0("chr",chr),]
	    vcf.chr <- vcf.chr[order(vcf.chr$POS),]
	    data.chr <- cnv.data.chr(chr=vcf.chr$CHROM, position=vcf.chr$POS, 
	                             normal.ref.rd=vcf.chr$Normal.REF.DP, normal.alt.rd=vcf.chr$Normal.ALT.DP, 
	                             tumor.ref.rd=vcf.chr$Tumor.REF.DP, tumor.alt.rd=vcf.chr$Tumor.ALT.DP,
	                             normal.total.coverage=normal.total.coverage, tumor.total.coverage=tumor.total.coverage)
    
        ## collect information
        data <- rbind(data, data.chr)
    
        if (verbose) cat("finish data transformation of chr =",chr,"\n")
    } ## for chr

    return(data)
}
