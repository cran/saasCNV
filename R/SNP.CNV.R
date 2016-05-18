SNP.CNV <-
function (snp, output.dir, sample.id,
          do.GC.adjust = FALSE,
          gc.file = system.file("extdata","GC_1kb_hg19.txt.gz",package="saasCNV"),
          min.chr.probe=100,
          min.snps=10,
          joint.segmentation.pvalue.cutoff=1e-4,
          max.chpts=30,
          do.merge=TRUE, use.null.data=TRUE, num.perm=1000, maxL=NULL, 
          merge.pvalue.cutoff=0.05,
          do.cnvcall.on.merge=TRUE, 
          cnvcall.pvalue.cutoff=0.05,
          do.boundary.refine=FALSE,
          do.plot=TRUE, cex=0.3, ref.num.probe=NULL,
          do.gene.anno=FALSE,
          gene.anno.file=NULL,
          seed=NULL,
          verbose=TRUE)
{
	## create output dir
	if (!file.exists(output.dir)) dir.create(output.dir)
	
	## data process
	snp.data <- snp.cnv.data(snp=snp, min.chr.probe=min.chr.probe, verbose=verbose) 

	if (do.GC.adjust) {
		gc <- read.delim(file = gc.file, as.is=TRUE)
        snp.data <- GC.adjust(data = snp.data, gc = gc, maxNumDataPoints = 10000)
	}
	
	data.dir <- paste0(output.dir,"/mid_data")
    if (!file.exists(data.dir)) dir.create(data.dir)
    save(list="snp.data", file=file.path(data.dir, "snp.data.RData"))
	
	## Step I: joint segmentation
	snp.segs <- joint.segmentation(data=snp.data[snp.data$use.in.seg==1,], 
	                               min.snps=min.snps, 
	                               global.pval.cutoff=joint.segmentation.pvalue.cutoff, 
	                               max.chpts=max.chpts, verbose=verbose)
    write.table(snp.segs, file=file.path(data.dir,"snp.segs.txt"), 
                quote=FALSE, row.names=FALSE, sep="\t")
    
    ## segmentation diagnosis plot
    if (do.plot) {
        plot.dir <- paste0(output.dir,"/mid_seg_plot")
        if (!file.exists(plot.dir)) dir.create(plot.dir)
        chrs <- sub("^chr","",unique(snp.segs$chr))
        
        for (chr in chrs) {
            png(filename=file.path(plot.dir, paste0("seg_chr",chr,".png")), width=240*5, height=240*4)
            diagnosis.seg.plot.chr(data=snp.data[snp.data$use.in.seg==1,], segs=snp.segs, sample.id=sample.id, chr=chr, cex=cex)
            dev.off()
        }
    } ## do.plot
    
    ## other diagnosis plot
    if (do.plot) {
        diagnosis.plot.dir <- file.path(output.dir,"mid_diagnosis_plot")
        if (!file.exists(diagnosis.plot.dir)) dir.create(diagnosis.plot.dir)
        chrs <- sub("^chr","",unique(snp.segs$chr))
        
        png(filename=file.path(diagnosis.plot.dir,"diag_QQ_plot.png"), width=240*3, height=240*2)
        diagnosis.QQ.plot(data=snp.data[snp.data$use.in.seg==1,], chrs=chrs)
        dev.off()    	
    } ##do.plot


    ## Step II (optional): merge segments not far apart in log2ratio and log2mBAF
    if (do.merge) {
    	snp.segs.merge <- merging.segments(data=snp.data[snp.data$use.in.seg==1,], segs.stat=snp.segs, 
    	                                 use.null.data=use.null.data,
                                         N=num.perm, maxL=maxL, 
                                         merge.pvalue.cutoff=merge.pvalue.cutoff, 
                                         verbose=verbose)

        write.table(snp.segs.merge, file=file.path(data.dir,"snp.segs.merge.txt"), 
                    quote=FALSE, row.names=FALSE, sep="\t")
        
        if (do.plot) {
            plot.dir <- paste0(output.dir,"/mid_seg_merge_plot")
            if (!file.exists(plot.dir)) dir.create(plot.dir)
            chrs <- sub("^chr","",unique(snp.segs.merge$chr))
        
            for (chr in chrs) {
                png(filename=file.path(plot.dir, paste0("merge_seg_chr",chr,".png")), width=240*5, height=240*4)
                diagnosis.seg.plot.chr(data=snp.data[snp.data$use.in.seg==1,], 
                                       segs=snp.segs.merge, sample.id=sample.id, chr=chr, cex=cex)
                dev.off()
            }
        } ## do.plot
    } ## do.merge


    ## Step III: call CNV
    cat("CNV calling ...\n")
    if (do.merge==TRUE & do.cnvcall.on.merge==TRUE) {
        snp.cnv <- cnv.call(data=snp.data[snp.data$use.in.seg==1,], sample.id=sample.id, segs.stat=snp.segs.merge, 
                            maxL=maxL, N=num.perm,
                            pvalue.cutoff=cnvcall.pvalue.cutoff)
    } else {
    	snp.cnv <- cnv.call(data=snp.data[snp.data$use.in.seg==1,], sample.id=sample.id, segs.stat=snp.segs, 
    	                    maxL=maxL, N=num.perm,
                            pvalue.cutoff=cnvcall.pvalue.cutoff)

    }
    res.dir <- paste0(output.dir,"/mid_res")
    if (!file.exists(res.dir)) dir.create(res.dir)
    write.table(snp.cnv, file=file.path(res.dir,"snp.cnv.txt"), 
                quote=FALSE, row.names=FALSE, sep="\t")

    ## Step IV (optional): refine boundary
    if (do.boundary.refine) {
    	cat("Refining boundary ...\n")
    	snp.cnv.refine <- snp.refine.boundary(data=snp.data, segs.stat=snp.cnv)
        write.table(snp.cnv.refine, file=file.path(res.dir,"snp.cnv.refine.txt"),
                    quote=FALSE, row.names=FALSE, sep="\t")
    }
    
    ## Step V (optional): gene annotation
    if (do.gene.anno) {
    	cat("gene annotation ...\n")
        if (do.boundary.refine) {
            gene <- read.delim(file=gene.anno.file, as.is=TRUE, comment.char="")
            snp.cnv.refine.anno <- reannotate.CNV.res(res=snp.cnv.refine, gene=gene, only.CNV=TRUE)
            write.table(snp.cnv.refine.anno, file=file.path(res.dir,"snp.cnv.refine.anno.txt"),
                        quote=FALSE, row.names=FALSE, sep="\t")
        } else {
            gene <- read.delim(file=gene.anno.file, as.is=TRUE, comment.char="")
            snp.cnv.anno <- reannotate.CNV.res(res=snp.cnv, gene=gene, only.CNV=TRUE)
            write.table(snp.cnv.anno, file=file.path(res.dir,"snp.cnv.anno.txt"),
                        quote=FALSE, row.names=FALSE, sep="\t")        	
        }
    }

    ## visualization
    if (do.boundary.refine) {
        ## genome-wide plot
        png(filename=file.path(res.dir,"cnv_gw_plot.png"), width=240*6, height=240*4)
        genome.wide.plot(data=snp.data, segs=snp.cnv.refine, 
                         sample.id=sample.id, 
                         chrs=sub("^chr","",unique(snp.cnv$chr)), 
                         cex=cex)
        dev.off()

        ## diagnosis cluster plot
        pdf(file=file.path(res.dir,"cnv_cluster_plot.pdf"), width=8, height=8)
        diagnosis.cluster.plot(segs=snp.cnv.refine, 
                               chrs=sub("^chr","",unique(snp.cnv$chr)), 
                               min.snps=min.snps, max.cex=3, 
                               ref.num.probe=ref.num.probe)
        dev.off()
    } else {
        ## genome-wide plot
        png(filename=file.path(res.dir,"cnv_gw_plot.png"), width=240*6, height=240*4)
        genome.wide.plot(data=snp.data[snp.data$use.in.seg==1,], segs=snp.cnv, 
                         sample.id=sample.id, 
                         chrs=sub("^chr","",unique(snp.cnv$chr)), 
                         cex=cex)
        dev.off()

        ## diagnosis cluster plot
        pdf(file=file.path(res.dir,"cnv_cluster_plot.pdf"), width=8, height=8)
        diagnosis.cluster.plot(segs=snp.cnv, 
                               chrs=sub("^chr","",unique(snp.cnv$chr)), 
                               min.snps=min.snps, max.cex=3, 
                               ref.num.probe=ref.num.probe)
        dev.off()
    }
}
