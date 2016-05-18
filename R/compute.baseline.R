compute.baseline <- 
function(segs.stat, data)
{
    refine.idx.select <- function(idx.select)
    {
        ## search on log2mBAF dimension
    	segs.stat1 <- segs.stat[idx.select,]
    	mode.tbl1 <- round(rep(segs.stat1$log2mBAF.Median,  segs.stat1$numProbe), 2)
    	base.mode.log2mBAF1  <- Mode( mode.tbl1 )
    	mode.tbl2 <- mode.tbl1[mode.tbl1 != base.mode.log2mBAF1]
    	if (length(mode.tbl2)>0) {
    		base.mode.log2mBAF2  <- Mode( mode.tbl2 )
    		peak1 <- sum(mode.tbl1 == base.mode.log2mBAF1)
    		peak2 <- sum(mode.tbl2 == base.mode.log2mBAF2)
            if (peak2 > 0.5*peak1 &
                abs(peak2 - peak1) >= 0.3*sigma.log2mBAF &
                abs(base.mode.log2mBAF2) < abs(base.mode.log2mBAF1)) {
            	base.mode.log2mBAF <- base.mode.log2mBAF2
            } else {
            	base.mode.log2mBAF <- base.mode.log2mBAF1
            }
    	} else {
    		base.mode.log2mBAF <- base.mode.log2mBAF1
    	}
    	tmp1 <- segs.stat1$log2mBAF.Median  - base.mode.log2mBAF    	
    	idx.select <- idx.select[ which( abs(tmp1) < 0.3*sigma.log2mBAF ) ]

    	## search on log2ratio dimension
    	segs.stat2 <- segs.stat[idx.select,]
    	base.mode.log2ratio <- Mode( round(rep(segs.stat2$log2ratio.Median, segs.stat2$numProbe), 2) )
    	tmp2 <- segs.stat2$log2ratio.Median - base.mode.log2ratio
        idx.select <- idx.select[ which( abs(tmp2) < 0.5*sigma.log2ratio ) ]
        
        return(idx.select)
    	# # seg.sigma.log2mBAF  <- mad(tmp1)
    	# # seg.sigma.log2ratio <- mad(tmp2)
    	# # idx.select <- idx.select[ which( abs(tmp1) < min(2*seg.sigma.log2mBAF,  sigma.log2mBAF) &
    	                                 # # abs(tmp2) < min(3*seg.sigma.log2ratio, sigma.log2ratio) ) ]   	
    }
    
    segs.stat <- segs.stat[!(segs.stat$chr %in% c("chrX","chrY")),]
    sigma.log2ratio <- delta.sd(data$log2ratio)
    sigma.log2mBAF  <- delta.sd(data$log2mBAF)

    idx.select <- which( abs(segs.stat$log2mBAF.Median) < 0.5 * sigma.log2mBAF )
    if ( length(idx.select)>0 ) { ## 0.5 sigma
    	idx.select <- refine.idx.select(idx.select)
    } else { ## 1 sigma
        idx.select <- which( abs(segs.stat$log2mBAF.Median) < sigma.log2mBAF )
    	if ( length(idx.select)>0 ) {
    		idx.select <- refine.idx.select(idx.select)
    	} else { ## 1.5 sigma
            idx.select <- which( abs(segs.stat$log2mBAF.Median) < 1.5 * sigma.log2mBAF )
    	    if ( length(idx.select)>0 ) {
    	    	idx.select <- refine.idx.select(idx.select)
    		} else {
    			idx.select <- refine.idx.select(1:nrow(segs.stat))
    		} 		
    	}
    }
    
    idx.base <- NULL
    for (i in idx.select) {
    	idx.base <- c(idx.base, which(data$chr==segs.stat$chr[i] & 
                                      data$pos>=segs.stat$posStart[i] & 
   	                                  data$pos<=segs.stat$posEnd[i]))
    }
    
    return(idx.base)	
}
