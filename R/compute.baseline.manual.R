compute.baseline.manual <- 
function(segs.stat, data, 
         log2mBAF.left=NULL, log2mBAF.right=NULL, log2ratio.bottom=NULL, log2ratio.up=NULL)
{
    segs.stat <- segs.stat[!(segs.stat$chr %in% c("chrX","chrY")),]
    sigma.log2ratio <- delta.sd(data$log2ratio)
    sigma.log2mBAF  <- delta.sd(data$log2mBAF)

    idx.select <- which( segs.stat$log2mBAF.Median > log2mBAF.left &
                         segs.stat$log2mBAF.Median < log2mBAF.right & 
                         segs.stat$log2ratio.Median > log2ratio.bottom &
                         segs.stat$log2ratio.Median < log2ratio.up )
       
    k <- 1
    while ( length(idx.select)==0 ) {
        k <- k + 0.5
        idx.select <- which( segs.stat$log2mBAF.Median > k*log2mBAF.left &
                             segs.stat$log2mBAF.Median < k*log2mBAF.right & 
                             segs.stat$log2ratio.Median > k*log2ratio.bottom &
                             segs.stat$log2ratio.Median < k*log2ratio.up )
        if (k>3) break
    }
    
    idx.base <- NULL
    if (length(idx.select)>0) {
    	for (i in idx.select) {
    	    idx.base <- c(idx.base, which(data$chr==segs.stat$chr[i] & 
                                          data$pos>=segs.stat$posStart[i] & 
   	                                      data$pos<=segs.stat$posEnd[i]))
        }
    }
    
    return(idx.base)	
}
