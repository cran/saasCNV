impute.missing.data <-
function(x)
{
	## require(RANN)
	idx.data <- which(!is.na(x))
    idx.na   <- which(is.na(x))
    tmp <- nn2(data=cbind(idx.data), query=cbind(idx.na), k=6)   ## k-nearest neighbor
    idx.nn <- idx.data[tmp$nn.idx]
    x[idx.na] <- rowMeans(matrix(x[idx.nn],nrow=length(idx.na)))
    return(x)
}
