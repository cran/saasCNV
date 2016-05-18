check.overlap <-
function(start,end, start1, end1) {
	## start, end -- vector, length(start)==length(end)
	## start1, end1 -- scalar value
	idx1 <- start>=start1 & start<=end1
	idx2 <- end>=start1 & end<=end1
	idx3 <- start<=start1 & end>=end1
	idx4 <- start>=start1 & end<=end1
    return(which(idx1 | idx2 | idx3 | idx4))
    }
