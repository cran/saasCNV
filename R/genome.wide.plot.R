genome.wide.plot <-
function(data, segs, sample.id, chrs, cex=0.3)
{
    ## data
    data1 <- data
    data1$chr <- sub("^chr","",data1$chr)
    data1 <- data1[data1$chr %in% chrs,]
    data1$chr[data1$chr=="X"] <- 23
    data1$chr[data1$chr=="Y"] <- 24
    data1$chr <- as.integer(data1$chr)
    data1 <- data1[order(data1$chr, data1$position),]
    
    ## segs
    segs1 <- segs
    segs1$chr <- sub("^chr","", segs1$chr)
    segs1 <- segs1[segs1$chr %in% chrs,]
    segs1$chr[segs1$chr=="X"] <- 23
    segs1$chr[segs1$chr=="Y"] <- 24
    segs1$chr <- as.integer(segs1$chr)
    segs1 <- segs1[order(segs1$chr, segs1$posStart),]
    
    ## new chr position
    position.offsets <- cumsum(c(0, head(tapply(data1$position, data1$chr, max), -1))) ## left ends
    chrs1 <- sort(unique(data1$chr))
    names(position.offsets) <- chrs1
    NewPos1 <- data1$position + rep(position.offsets, tapply(data1$position, data1$chr, length))    
    
    ## ticks position
    xTicks <- chrs1
    xNames <- chrs1
    names( xTicks ) <- xTicks
    xNames.name <- as.character(xTicks)
    xNames.name[xNames.name=="23"] <- "X"
    xNames.name[xNames.name=="24"] <- "Y"
    names( xNames ) <- xNames.name
    for ( j in seq_along( xTicks ) ) {
        xTicks[j] <- max( NewPos1[data1$chr == xTicks[j]] )
        xNames[j] <- ( min( NewPos1[data1$chr == xNames[j]] ) + max( NewPos1[data1$chr == xNames[j]] ) ) / 2
        }
    xTicks <- c( min( NewPos1 ) , xTicks )
    nchr <- length(chrs1)
    
    ## colors
    cols1 <- rep(NA,length(data1$log2ratio))
    if (nchr %% 2 == 0) {
    	cols1[data1$chr %in% seq(1, nchr-1,  by=2)] <- "lightgreen"
        cols1[data1$chr %in% seq(2, nchr,by=2)] <- "lightgray"
    } else {
    	cols1[data1$chr %in% seq(1, nchr,  by=2)] <- "lightgreen"
        cols1[data1$chr %in% seq(2, nchr-1,by=2)] <- "lightgray"
    }

    ## segment colors
    cols.seg <- c("blue","darkgray","red","darkgreen","cyan")
    names(cols.seg) <- c("loss","normal","gain","LOH","undecided")
        

    oldpar <- par(mfrow=c(4,1))
    
    ## log2ratio
    par(mar=c(4,4,4,2))
    plot(x=NewPos1,y=data1$log2ratio, ylim=c(-4, 4), col=cols1, pch=15, cex=cex, xaxt="n", 
         xlab="Chromosome", ylab="log2ratio", main=sample.id)
    abline(h=0, col="black", lty=3,lwd=1)
    axis( 1 , xTicks , labels = FALSE )
    axis( 1 , xNames , names( xNames ) , tick = FALSE )
    ## segments
    for (j in 1:nrow(segs1)) {
    	seg <- segs1[j,]
    	x1 <- position.offsets[seg$chr[1]] + seg$posStart[1]
    	x2 <- position.offsets[seg$chr[1]] + seg$posEnd[1]
    	y <- seg$log2ratio.Median[1]
    	type <- seg$CNV[1]
        segments(x0=x1,y0=y,x1=x2,y1=y, col=cols.seg[type], lty=1,lwd=3)
        }
    
    ## log2mBAF
    par(mar=c(4,4,1,2))
    plot(x=NewPos1,y=data1$log2mBAF, ylim=c(-1, 1), col=cols1, pch=15, cex=cex, xaxt="n", 
         xlab="Chromosome", ylab="log2mBAF", main="")
    abline(h=0, col="black", lty=3,lwd=1)
    axis( 1 , xTicks , labels = FALSE )
    axis( 1 , xNames , names( xNames ) , tick = FALSE )
    ## segments
    for (j in 1:nrow(segs1)) {
    	seg <- segs1[j,]
    	x1 <- position.offsets[seg$chr[1]] + seg$posStart[1]
    	x2 <- position.offsets[seg$chr[1]] + seg$posEnd[1]
    	y <- seg$log2mBAF.Median[1]
    	type <- seg$CNV[1]
        segments(x0=x1,y0=y,x1=x2,y1=y, col=cols.seg[type], lty=1,lwd=3)
        }

    ## tumor BAF
    par(mar=c(4,4,1,2))
    plot(x=NewPos1,y=data1$tumor.BAF, ylim=c(0, 1), col=cols1, pch=15, cex=cex, xaxt="n", 
         xlab="Chromosome", ylab="Tumor BAF", main="")
    abline(h=c(0,0.5,1), col="black", lty=3,lwd=1)
    axis( 1 , xTicks , labels = FALSE )
    axis( 1 , xNames , names( xNames ) , tick = FALSE )
    
    ## normal BAF
    par(mar=c(4,4,1,2))
    plot(x=NewPos1,y=data1$normal.BAF, ylim=c(0, 1), col=cols1, pch=15, cex=cex, xaxt="n", 
         xlab="Chromosome", ylab="Normal BAF", main="")
    abline(h=c(0,0.5,1), col="black", lty=3,lwd=1)
    axis( 1 , xTicks , labels = FALSE )
    axis( 1 , xNames , names( xNames ) , tick = FALSE )
     
    par(oldpar)
}
