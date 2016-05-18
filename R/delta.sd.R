delta.sd <-
function(x)
{
	n <- length(x)
	n <- ifelse(n %% 2==0, n, n-1)
	tmp <- x[seq(2,n,by=2)] - x[seq(1,n-1,by=2)]
	return(mad(tmp, na.rm=TRUE)/sqrt(2))
}
