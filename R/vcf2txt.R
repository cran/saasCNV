vcf2txt <-
function(vcf.file, normal.col=10, tumor.col=11, MQ.cutoff=30)
{
	inputFile <- file(vcf.file, "r")

    ## skip some definition lines
    while(length(lines <- readLines(inputFile, n = 1, warn = FALSE)) > 0) {
        ##message(lines)
        if (length(grep("\\#CHROM",lines))==1) break
    }
    header <- sub("\\#","",lines)
    header <- strsplit(header,"\t")[[1]]
    header[c(normal.col, tumor.col)] <- c("Normal","Tumor")

    ## read in data table
    vcf <- read.table(inputFile,header=FALSE,sep="\t",as.is=TRUE,comment.char="") ##na.strings="."
    names(vcf) <- header
    vcf <- vcf[,c(1:9, normal.col, tumor.col)]
    close(inputFile)

    ## filter1: chr1-22,X,Y and PASS
    chrs <- paste0("chr",c(1:22,"X","Y"))
    idx <- with(vcf, CHROM %in% chrs & FILTER=="PASS")
    vcf <- vcf[idx,]

    ## filter2: only have 1 ALT allele
    idx <- grep(",",vcf$ALT)
    if (length(idx)>=1) vcf <- vcf[-idx,]

    ## find REF and ALT read depth
    ## normal
    normal.RD <- strsplit(vcf$Normal,":")
    normal.RD <- do.call(rbind,normal.RD)
    vcf$Normal.GT <- normal.RD[,1]
    normal.AD <- do.call(rbind,strsplit(normal.RD[,2],","))
    vcf$Normal.REF.DP <- as.integer(normal.AD[,1])
    vcf$Normal.ALT.DP <- as.integer(normal.AD[,2])
    ## tumor
    tumor.RD <- strsplit(vcf$Tumor,":")
    tumor.RD <- do.call(rbind,tumor.RD)
    vcf$Tumor.GT <- tumor.RD[,1]
    tumor.AD <- do.call(rbind,strsplit(tumor.RD[,2],","))
    vcf$Tumor.REF.DP <- as.integer(tumor.AD[,1])
    vcf$Tumor.ALT.DP <- as.integer(tumor.AD[,2])

    ## reduce # columns
    vcf <- vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "INFO", 
                  "Normal.GT", "Normal.REF.DP", "Normal.ALT.DP", "Tumor.GT", "Tumor.REF.DP", "Tumor.ALT.DP")]

    ## filter3: normal or tumor GT=="./."
    idx <- which(vcf$Normal.GT=="./." | vcf$Tumor.GT=="./.")
    if (length(idx)>=1) vcf <- vcf[-idx,]

    ## in INFO find MQ
    info <- strsplit(vcf$INFO,";")
    extract.MQ <- function(x) {
        idx.MQ <- grep("^MQ\\=",x)
        if (length(idx.MQ)==1) return(x[idx.MQ]) else return(NA)
        }
    info1 <- sapply(info,FUN=extract.MQ)
    vcf$MQ <- as.numeric(sub("^MQ\\=","",info1))

    ## filter4: MQ
    vcf <- vcf[vcf$MQ>=MQ.cutoff,]

    ## reduce # columns
    vcf <- vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "MQ", 
                  "Normal.GT", "Normal.REF.DP", "Normal.ALT.DP", "Tumor.GT", "Tumor.REF.DP", "Tumor.ALT.DP")]

    ## remove duplicate
    a <- vcf[,c("CHROM","POS")]
    idx <- which(duplicated(a))
    if (length(idx)>=1) vcf <- vcf[-idx,]
    
    return(vcf)
}
