reannotate.CNV.res <-
function(res, gene, only.CNV=FALSE) {
    ## preprocessing
    ## gene annotation
    gene <- gene[,c("name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","name2","cdsStartStat","cdsEndStat")]
    gene <- gene[gene$chrom %in% paste("chr",c(1:22,"X","Y"),sep=""),]
    gene <- gene[order(gene$chrom,gene$txStart,gene$strand),]
    
    ## GISTIC results
    if (only.CNV==TRUE) {
        res <- res[!is.na(res$CNV) & res$CNV!="normal" & res$CNV!="undecided", 
                   c("Sample_ID", "chr", "posStart", "posEnd", "length", "numProbe",
                     "log2ratio.Mean.adj", "log2ratio.SD", "log2ratio.Median.adj", "log2ratio.MAD", 
                     "log2mBAF.Mean", "log2mBAF.SD", "log2mBAF.Median", "log2mBAF.MAD", 
                     "log2ratio.p.value", "log2mBAF.p.value", "p.value", "CNV")]
    } else {
    	res <- res[,c("Sample_ID", "chr", "posStart", "posEnd", "length", "numProbe",
                      "log2ratio.Mean.adj", "log2ratio.SD", "log2ratio.Median.adj", "log2ratio.MAD", 
                      "log2mBAF.Mean", "log2mBAF.SD", "log2mBAF.Median", "log2mBAF.MAD", 
                      "log2ratio.p.value", "log2mBAF.p.value", "p.value", "CNV")]
    }
    if (length(grep("chr",res$chr))==0) res$chr <- paste("chr",res$chr,sep="")
    res$numGene <- 0
    res$gene <- NA

    ## re-annotate gene list in each region
    for (i in 1:nrow(res)) {
    	chr <- res$chr[i]
        start <- res$posStart[i]
        end <- res$posEnd[i]
        
        gene.chr <- gene[gene$chrom==chr,]
        idx.gene <- check.overlap(start=gene.chr$txStart, end=gene.chr$txEnd, start1=start, end1=end)
        genes <- sort(unique(gene.chr$name2[idx.gene]))
        res$numGene[i] <- length(genes)
        if (res$numGene[i]==1) res$gene[i] <- genes 
        if (res$numGene[i]>1) res$gene[i] <- Reduce(f=function(x,y) paste(x,y,sep=","), x=genes) 
        } ## for i
    
    return(res)
    }
