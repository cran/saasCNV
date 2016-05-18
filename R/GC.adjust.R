GC.adjust <- 
function(data, gc, maxNumDataPoints = 10000)
{
    data$marker <- paste0(data$chr, ":", round(data$position/1000,0)*1000+1)
    gc$marker <- paste0(gc$chr, ":", gc$position)
    data <- merge(x=data, y=gc[,c("marker","GC")], by="marker", all=FALSE)
    data <- data[order(data$chr, data$pos),]
    data <- data[,-which(names(data)=="marker")]

    data$log2ratio.woGCAdj <- data$log2ratio
    ratio <- 2^data$log2ratio
    idx <- which(!is.na(data$GC) & !is.na(ratio) & ratio > 0)
    data <- data[idx,]
    ratio <- ratio[idx]
    
    if (nrow(data) > maxNumDataPoints) {
        idx <- sample(1:nrow(data), maxNumDataPoints)
    }
    gcData <- data.frame(gc = data$GC[idx], ratio = ratio[idx])
    gc.fit <- loess(ratio ~ gc, gcData)
    normCoef <- predict(gc.fit, data.frame(gc = data$GC)) 
    normCoef <- normCoef/median(normCoef, na.rm = TRUE)
    ratio <- ratio/normCoef
    data$log2ratio <- log2(ratio)
    
    idx <- which(!is.na(data$log2ratio) & data$log2ratio < Inf & data$log2ratio > -Inf)
    data <- data[idx,]
    
    return(data)
}