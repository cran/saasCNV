diagnosis.QQ.plot <-
function(data, chrs)
{
	data1 <- data
	data1$chr <- sub("^chr","",data1$chr)
    data1 <- data1[data1$chr %in% chrs,]
    data1$chr[data1$chr=="X"] <- 23
    data1$chr[data1$chr=="Y"] <- 24
    data1$chr <- as.integer(data1$chr)
    data1 <- data1[order(data1$chr, data1$position),]
    
    oldpar <- par(mfrow=c(1,2))
    ##log2ratio
    qqnorm(data1$log2ratio, main="log2ratio")
    qqline(data1$log2ratio, col="red")
    ##log2mBAF
    qqnorm(data1$log2mBAF, main="log2mBAF")
    qqline(data1$log2mBAF, col="red")
    par(oldpar)
}
