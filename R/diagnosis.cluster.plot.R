diagnosis.cluster.plot <-
function(segs, chrs, min.snps, max.cex=3, ref.num.probe=NULL)
{
	segs1 <- segs
    segs1$chr <- sub("^chr","", segs1$chr)
    segs1 <- segs1[segs1$chr %in% chrs,]
    segs1$chr[segs1$chr=="X"] <- 23
    segs1$chr[segs1$chr=="Y"] <- 24
    segs1$chr <- as.integer(segs1$chr)
    segs1 <- segs1[order(segs1$chr, segs1$posStart),]
    
    total.num.probe <- sum(segs1$numProbe, na.rm=TRUE)
    ref.num.probe <- ifelse(is.null(ref.num.probe), floor(total.num.probe/100), ref.num.probe)

    ## segment colors
    cols.seg <- c("blue","darkgray","red","darkgreen","cyan")
    names(cols.seg) <- c("loss","normal","gain","LOH","undecided")

    def.par <- par(no.readonly = TRUE) 

    nf <- layout(matrix(c(2,4,1,3),2,2,byrow=TRUE), c(3,1), c(1,3))
    ##layout.show(nf)
    
    par(mar=c(4,4,1,1))
    cols1 <- cols.seg[segs1$CNV]
    cexs1 <- segs1$numProbe/ref.num.probe
    tmp <- cut(cexs1, sort(unique(c(0, 1e-3, 0.01, seq(0.1,1,by=0.1),max(cexs1)))), include.lowest = TRUE)
    rate <- c(0.1,0.1, seq(0.1,1,by=0.1),1)
    names(rate) <- levels(tmp)
    cexs1 <- rate[tmp] * max.cex
    plot(segs1$log2mBAF.Median, segs1$log2ratio.Median, col=cols1, cex=cexs1, lwd=2,
         xlab="log2mBAF", ylab="log2ratio")
    grid()
    abline(h=0, v=0, lty=3, lwd=2, col="black")
    abline(h=unique(segs1$log2ratio.base.Median), v=unique(segs1$log2mBAF.base.Median), lty=2, lwd=2, col="darkgray")
    idx.min <- which(segs$numProbe==min.snps)
    points(segs1$log2mBAF.Median[idx.min], segs1$log2ratio.Median[idx.min], col="black", cex=1)
    ##legend("topleft",legend=names(cols.seg), col=cols.seg, pch=1)

    par(mar=c(0,4,1,1))
    a <- rep(segs1$log2mBAF.Median, segs1$numProbe)
    a <- round( a[!is.na(a)], 2)
    b <- 100 ##seq(min(a),max(a),by=0.01)
    xhist <- hist(a, b, plot=FALSE)
    barplot(xhist$counts, xaxt="n", space=0)
    ##hist(a,b, xaxt="n", main="")
    ##abline(v=c(0, unique(segs1$log2mBAF.base.Median)), lty=c(3,2), col=c("black","darkgray"))
    
    par(mar=c(4,0,1,1))
    a <- rep(segs1$log2ratio.Median, segs1$numProbe)
    a <- round( a[!is.na(a)], 2)
    b <- 100 ##seq(min(a),max(a),by=0.01)
    yhist <- hist(a, b, plot=FALSE)
    barplot(yhist$counts, yaxt="n", space=0, horiz=TRUE)
    
    ## legend
    par(mar=c(0,0,2,0))
    plot(x=0:1,y=0:1, type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
    legend("topleft", legend=names(cols.seg), col=cols.seg, pch=1, cex=1.5, lty=0, lwd=3)   
    
    par(def.par)
}
