`dchi` <- function(y,N)
{
    fy <- (1-N/2)*log(2) + (N-1)*log(y) - y^2/2 - lgamma(N/2)
    exp(fy)
}


`compute.var` <- function(y)
{
    T <- nrow(y)
    T <- ifelse(T %% 2==0, T, T-1)
    y.diff <- y[seq(1,T-1,by=2),] - y[seq(2,T,by=2),]
    y.var  <- apply(y.diff, 2, mad)^2/2
    return( y.var )
}


`matrix.max` <- function(Z)
{
    ## If Z has NA values treat as -Inf
    na.inds = which(is.na(Z),arr.ind=TRUE)
    Z[na.inds] = -Inf

    row.max.col = max.col(Z, ties.method=c("random", "first", "last"))
    max.row = which.max(Z[cbind(1:nrow(Z),row.max.col)])
    c(max.row, row.max.col[max.row])
}


`computeBeta` <- function(ALPHA)
{
    w<-function(u){
         exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    integrand<-function(u){
        u^2*w(u)^2*(2*u^4*w(u)*(1-w(u))+5*u^2*w(u)-3*u^2-2)*dchi(u,1)
    }
        
    numerator=integrate(integrand,lower=0,upper=10)
    
    g<-function(u){
         u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    g.moments=computeMoments(g,0)
    numerator$value/(2*g.moments$psidotdot)
}


`computeMoments` <- function(g, theta)
{

    INTLIM.THRESH = 0.001
    psidottop.int <- function(u){ g(u)*exp(theta*g(u))*exp(-(u)^2/2) }    
    for( INTLIM in 10:50 ){
        temp=psidottop.int(INTLIM)
        if (is.na(temp)){
            INTLIM=INTLIM-1
            break
        }
        if (temp<INTLIM.THRESH) break
    }
    psidottop = integrate(psidottop.int,lower=-INTLIM,upper=INTLIM)
    psidottop = psidottop$value/sqrt(2*pi)

    psidotbot.int =  function(u){ exp(theta*g(u))*exp(-(u)^2/2)}  
    for( INTLIM in 10:50 ){
        temp=psidotbot.int(INTLIM)
        if( is.na(temp)){
            INTLIM=INTLIM-1
            break
        }
        if( temp<INTLIM.THRESH ) break
    }
    psidotbot = integrate(psidotbot.int, lower=-INTLIM, upper=INTLIM)
    psidotbot = psidotbot$value/sqrt(2*pi)
    psidot = psidottop/psidotbot
    
    EgUsq.int<-function(u){ g(u)^2*exp(theta*g(u))*exp(-(u)^2/2) }    
    for( INTLIM in 10:50 ){
        temp=EgUsq.int(INTLIM)
        if( is.na(temp)){
            INTLIM=INTLIM-1
            break
        }
        if( temp<INTLIM.THRESH) break
    }
    EgUsq =  integrate(EgUsq.int,lower=-INTLIM,upper=INTLIM)
    EgUsq = EgUsq$value/sqrt(2*pi)

    psidotdot = (EgUsq*psidotbot - psidottop^2)/(psidotbot^2);
    psi = psidotbot

    list(psi=psi,psidot=psidot,psidotdot=psidotdot)
}


`vu` <- function(x, maxn=1000, do.approx=(abs(x)<0.2))
{
    x=as.matrix(x)
    vux = matrix(0,nrow=nrow(x),ncol=ncol(x))

    if(is.logical(do.approx)){
        if(sum(do.approx) > 0)   vux[do.approx] = exp(-x[do.approx]*0.583)
    } else if(length(do.approx) > 0){
        vux[do.approx] = exp(-x[do.approx]*0.583)
    }
  

    if (sum(do.approx)<length(x)){
        notdo.approx.ix = which(!do.approx)
        n = matrix(c(1:maxn), nrow=1, ncol=maxn)
        summands = pnorm( -0.5*matrix(x[notdo.approx.ix], nrow=length(notdo.approx.ix), ncol=1) %*% sqrt(n) )/
                   matrix(rep(n,length(notdo.approx.ix)), ncol=length(n), byrow=TRUE)
        expterm = -2*apply(summands,1,"sum")
        vux[notdo.approx.ix] = (2/x[notdo.approx.ix]^2)*exp(expterm)
    }

    vux
}


`computeTiltDirect` <- function(b,g,THRESH,theta0) 
{
    theta = theta0 # if start theta at 0, then dh(theta)=0 at the first step for g(u)=u^2.
    prevtheta = Inf
    prevprevtheta = Inf
    thetarec = theta
    
    while( abs(theta-prevtheta)>THRESH && abs(theta-prevprevtheta)>THRESH ){
        
        g.moments <- computeMoments(g,theta)
        htheta = g.moments$psidot-b
        dhtheta = g.moments$psidotdot

        ## cat("theta=", theta," htheta=",htheta,".\n",sep="")
        thetarec = c(thetarec, theta)
        prevprevtheta = prevtheta
        prevtheta = theta
        theta = prevtheta - htheta/dhtheta
    }

    theta=prevtheta
    psi = g.moments$psi

    list(theta=theta,psi=psi,psidot=g.moments$psidot,psidotdot=g.moments$psidotdot)
}


`ComputeZ.fromS.R.partial` <- function(this.S, this.SST, this.imap, start.inds, end.inds, ALPHA, MIN.SNPs)
{
    T <- nrow(this.S)
    N <- ncol(this.S)
    dfnum <- 1

    log.alpha <- log(ALPHA)
    g <- function(u){
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    
    Z <- matrix(nrow=length(start.inds), ncol=length(end.inds), data=0)

    totalsnps <- this.imap[T,] - this.imap[1,]
    dfden     <- totalsnps-2
    for (i in 1:length(start.inds)) {
        for (j in 1:length(end.inds)) {
            st <- start.inds[i]
            ed <- end.inds[j]
            if (st > 0 && st < ed && ed<T) {
                nsnps <- this.imap[ed,] - this.imap[st,]
                diff1 <- this.S[ed,] - this.S[st,]
                SSb <- nsnps*(diff1/nsnps - this.S[T,]/totalsnps)^2
                SSb <- SSb + (totalsnps-nsnps)*((this.S[T,]-diff1)/(totalsnps-nsnps)-this.S[T,]/totalsnps)^2
                SSw <- this.SST - SSb
                U <- (SSb/dfnum)/(SSw/dfden)
                set.to.zero <- which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                U[set.to.zero] <- 0
                Z[i,j] <- sum(g(U))
            }
        }
    }

    return(Z)
}


`ComputeZ.fromS.R` <- function(this.S, this.SST, this.imap, win, ALPHA, MIN.SNPs)
{
    T <- nrow(this.S)
    N <- ncol(this.S)
    dfnum <- 1

    win = min(win, T-1)
    if (win==0) { 
        return(NULL)
    }

    log.alpha <- log(ALPHA)
    g <- function(u) {
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    
    U <- matrix(nrow=T, ncol=win, data=0)
    Z <- U
    for (i in 1:N) {
        totalsnps <- this.imap[T,i] - this.imap[1,i]
        dfden     <- totalsnps - 2
        for (k in 1:win) {
            nsnps <- this.imap[(k+1):T,i] - this.imap[1:(T-k),i]
            diff1 <- this.S[(k+1):T,i] - this.S[1:(T-k),i]
                
            SSb <- nsnps*(diff1/nsnps - this.S[T,i]/totalsnps)^2
            SSb <- SSb + (totalsnps - nsnps)*((this.S[T,i] - diff1)/(totalsnps - nsnps) - this.S[T,i]/totalsnps)^2
            SSw <- this.SST[i] - SSb
            U[1:(T-k),k] <- (SSb/dfnum)/(SSw/dfden)
                
            set.to.zero <- which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
            U[set.to.zero,k] <- 0
        }
        Z <- Z + g(U)
    }
    return(Z)
}


`computeZ.onechange` <- function(y, t, y.var)
{
    T    <- ncol(y)
    St   <- apply(y[,1:t, drop=FALSE], 1, sum)
    ST   <- apply(y, 1, sum)
    Zsq  <- (St - (t/T)*ST)^2/(y.var*t*(1-t/T))
    sumZ <- sum(Zsq)
    pval <- pchisq(sumZ, nrow(y), lower.tail=FALSE)
  
    list(sumZ=sumZ, pval=pval)
}


`computeZ.onechange.sample` <- function(y, t, y.var)
{
    T    <- ncol(y)
    St   <- apply(y[,1:t, drop=FALSE], 1, sum)
    ST   <- apply(y,1,sum)
    Zsq  <- (St - (t/T)*ST)^2/(y.var*t*(1-t/T))
    pval <- pchisq(Zsq, 1, lower.tail=FALSE)

    list(Zsq=Zsq, pval=pval)
}


`computeZ.squarewave.sample` <- function(this.y, seg, y.var)
{
    T  <- ncol(this.y)
    Ss <- apply(as.matrix(this.y[,1:seg[1]], nrow=nrow(this.y), ncol=seg[1]), 1, sum)
    St <- apply(this.y[,1:seg[2]], 1, sum)
    ST <- apply(this.y,1,sum)
    k  <- seg[2]-seg[1]
    Zsq <- (St-Ss - (k/T)*ST)^2/(y.var*k*(1-k/T))
    pval <- pchisq(Zsq, 1, lower.tail=FALSE)

    list(Zsq=Zsq, pval=pval)
}


`fcompute.max.Z` <- function(this.y, win, y.var, ALPHA, MIN.SNPs, SINGLECHANGE.THRESH=0.0001)
{
    N      <- ncol(this.y)
    this.T <- nrow(this.y)

    if (this.T < 2*MIN.SNPs) {
        bestZ  <- NA
        bestchpt <- c(NA, NA)
    } else {
        this.S   <- apply(this.y,2,cumsum)
        this.SST <- apply(this.y,2,var)*(this.T-1)
        
        ## Find (t1, t2] that maximizes Z using filtered scan
        temp <- fscan.max(this.S, this.SST, this.imap=NULL, MIN.SNPs=MIN.SNPs, ALPHA=ALPHA)
        bestchpt <- temp$seg
        bestZ    <- temp$maxZ
        
        ## Test the left and right change-point individually
        if (bestZ > 0) {
            if (bestchpt[1]==0) { ## change 1 => 0
                pval.L <- 1
            } else {
                pval.L <- computeZ.onechange(t(this.y[1:bestchpt[2],]), bestchpt[1], y.var)$pval
            }
            
            if (bestchpt[2]==this.T) {
                pval.R <- 1
            } else {
                pval.R <- computeZ.onechange(t(this.y[(bestchpt[1]+1):this.T,]), bestchpt[2]-bestchpt[1], y.var)$pval
            }
            
            if ( (pval.L > SINGLECHANGE.THRESH || bestchpt[1] < MIN.SNPs) &&
                 (pval.R <= SINGLECHANGE.THRESH && (this.T - bestchpt[2]) >= MIN.SNPs) ) {
                bestchpt <- c(bestchpt[2], NA)
            } else {
                if ( (pval.R > SINGLECHANGE.THRESH || (this.T - bestchpt[2]) < MIN.SNPs) &&
                     (pval.L <= SINGLECHANGE.THRESH && bestchpt[1] >= MIN.SNPs) ) {
                    bestchpt <- c(bestchpt[1], NA)
                } else {
                    if ( (pval.L > SINGLECHANGE.THRESH || bestchpt[1] < MIN.SNPs) &&
                         (pval.R > SINGLECHANGE.THRESH || (this.T - bestchpt[2]) < MIN.SNPs) ) {
                        bestchpt <- c(NA, NA)
                        bestZ <- NA
                    }
                }
            }     
        
        } else {
            bestchpt <- c(NA, NA)
            bestZ <- NA
        }
    }
    
    list(bestchpt=bestchpt, bestZ=bestZ)
}


`fscan.max` <- function(this.S, this.SST, y.var=NULL, this.imap=NULL, 
                        f=NULL, MIN.SNPs=2, ALPHA=0, beta=NULL, verbose=FALSE)
{
    T <- nrow(this.S)   ## Number of SNPs
    N <- ncol(this.S)   ## Number of samples
    
    if (T < 2*MIN.SNPs) {
        maxZ <- NA
        seg  <- c(NA, NA)
    } else {

        if (is.null(f)) {
            f.power <- seq(-floor(log10(T))+1, 0, 1)
            f <- 10^f.power
        }
      
        if (T<1000) f <- 1  ## even if a value for f is passed in, do not allow filtering for short sequences
      
        f <- sort(f)
        if (f[length(f)] != 1) f <- c(f,1)
        R <- length(f)
        L <- ceiling(f[2:R]/f[1:(R-1)])
        L <- c(ceiling(T*f[1]), L)
      
        chpts   <- matrix(ncol=2, nrow=0)
        chpts.Z <- matrix(nrow=1, ncol=0)
    
        if (is.null(this.imap)) {
            this.imap <- matrix(1:T, nrow=T, ncol=N)
        }
    
        g2 <- function(u) { u^2*exp(u^2/2)/(ALPHA+exp(u^2/2)) }
        g2.moments <- computeMoments(g2,0)
        if(is.null(beta)) beta <- rep(1,N)
      
        for (r in 1:R) {
            stepsize <- floor(1/f[r])
        
            t <- seq(0, T, stepsize)         ## t is the filtered anchor set
            if (t[length(t)]<T) t <- c(t,T)  ## always include the last datapoint in the set
    
            if (min(t)==0) {
            	f.S    <- rbind(c(0,0), this.S[t,])
            	f.imap <- rbind(c(0,0), this.imap[t,])
            } else {
            	f.S    <- this.S[t,]
            	f.imap <- this.imap[t,]
            }
                
            ## Refine previously found change-points using the denser anchor set
            if (nrow(chpts)>0) {
                for (i in 1:nrow(chpts)) {
                    ind.L <- chpts[i,1] %/% stepsize
                    ind.R <- chpts[i,2] %/% stepsize
            
                    check.win  <- f[r]/f[r-1]
                    start.inds <- (ind.L - check.win):(ind.L + check.win) 
                    end.inds   <- (ind.R - check.win):(ind.R + check.win)

                    Z.part <- ComputeZ.fromS.R.partial(f.S, this.SST, f.imap, start.inds, end.inds, ALPHA, MIN.SNPs)
                    Z.part <- (Z.part - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)          
    
                    maxind      <- matrix.max(Z.part)
                    improved.cp <- c(t[start.inds[maxind[1]]], t[end.inds[maxind[2]]])
                    improved.Z  <- Z.part[maxind[1], maxind[2]]
                    if (verbose) {
                        cat("fscan.max: Changepoints (", chpts[i,1], ", ",chpts[i,2], 
                            ") refined to (", improved.cp[1], ", ",improved.cp[2], 
                            "). Z-score improved from ", chpts.Z[i], " to ", improved.Z,"\n", sep="")
                    }
                    chpts[i,]  <- improved.cp
                    chpts.Z[i] <- improved.Z
                } ## for i loop
            } ## if (nrow(chpts)>0)
        
            ## Do scan with min window size MIN.SNPS and max window size L[r]
            Z <- ComputeZ.fromS.R(f.S, this.SST, f.imap, L[r], ALPHA, MIN.SNPs)
            Z <- (Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)      
    
            maxind  <- matrix.max(Z)
            newchpt <- c(t[maxind[1]], t[maxind[1] + maxind[2]])
            newZ    <- Z[maxind[1], maxind[2]]
            if (verbose) {
                cat("fscan.max: Change-point found in round ", r, ":\n")
                print(c(newchpt, newZ), digits=1)
            }
            
            if (nrow(chpts)==0) {
                chpts <- rbind(chpts, newchpt)
                if (length(chpts) == 2) chpts <- matrix(nrow=1, ncol=2, data=chpts)      
                chpts.Z <- c(chpts.Z, newZ)
            }
            if (nrow(chpts)>=1) {
            	tmp <- unique(rbind(chpts, newchpt))
            	if (nrow(tmp) == nrow(chpts) + 1) {
            		chpts <- rbind(chpts, newchpt)
            		chpts.Z <- c(chpts.Z, newZ)
            	}
            	## otherwise the newchpt have been detected in previous round
            }
            
        } ## for r loop
    
        ## Take the maximum over the rounds, return chpt and Z score
        max.ind <- which.max(chpts.Z)
        maxZ    <- chpts.Z[max.ind]
        seg     <- chpts[max.ind,]
    } ## if T<=MIN.SNPS
    
    list(maxZ = maxZ, seg = seg)
}


`getCutoffMultisampleWeightedChisq` <- function(pval, m, delta, win, N, alpha)
{
    
    cat("Computing threshold for weighted chi-square...\n")
    THRES <- 0.1*pval
    currb <- 1
    prevsmallerb <- currb
    prevlargerb  <- 200
    currpval <- 1
    
    while( abs(currpval-pval)>THRES ) {        
        if( currpval > pval){
            ## need to increase b.
            prevsmallerb <- currb
            currb <- currb + (prevlargerb-currb)/2
        } else {
            ## need to decrease b.
            prevlargerb <- currb
            currb <- currb - (currb-prevsmallerb)/2
        }

        currpval <- pvalueMultisampleWeightedChisq(currb, m, delta, win, alpha, N)
    }    
    
    return( currb )
}


`pmarg.sumweightedchisq` <- function(b, ALPHA, N)
{
    g <- function(u) {
         u^2*exp(u^2/2)/(ALPHA + exp(u^2/2))
    }
    g.moments = computeMoments(g,0)
    gnormed <- function(u) {
    	(g(u)-g.moments$psidot)/sqrt(g.moments$psidotdot)
    }
    theta0 = sqrt(g.moments$psidotdot)/2   ## need theta/sqrt(psidotdot) < 1/2 for psi(theta)<infty
    b0 = b/sqrt(N)
    THRESH = 0.0001
    marg = computeTiltDirect(b0, gnormed, THRESH, theta0)
    thetabN = marg$theta*sqrt(N)
    pmarg = exp(-thetabN*b + N*log(marg$psi))/sqrt(2*pi*marg$psidotdot)
    pmarg
}


`pvalueMultisampleWeightedChisq` <- function(b, m, delta, win, ALPHA, N)
{
    delta1 = min(1, win/m)
    beta = computeBeta(ALPHA)

    integrand<-function(u){
        vu(sqrt(2)*b*sqrt(beta)/(sqrt(m)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral <- integrate(integrand, lower=delta, upper=delta1)
    pmarg    <- pmarg.sumweightedchisq(b, ALPHA, N)

    pval <- b^3*beta^2*pmarg*integral$value
    pval
}


`mscbs.classify` <- function(this.y, seg, y.var, CHISQ.PVAL.THRESH=0.001, MIN.REQ.ABSDIFF=NA, MIN.SUFF.ABSDIFF=NA)
{    
    N <- ncol(this.y)
    T <- nrow(this.y)

    if (is.na(MIN.REQ.ABSDIFF)) {
        MIN.REQ.ABSDIFF <- 0.5*sqrt(y.var)
    }
    if (is.na(MIN.SUFF.ABSDIFF)) {
        MIN.SUFF.ABSDIFF <- 5*sqrt(y.var)
    }
    
    ## segment: (seg[1],seg[2]]    
    if (is.na(seg[1]) && is.na(seg[2])) {
        ## Invalid change-point
        sampleswithsegment <- rep(0,N)
        this.yhat <- matrix(0, nrow=nrow(this.y), ncol=ncol(this.y))
    } else {
        if (is.na(seg[2])) {
            ## Single change-point
            sample.pval <- computeZ.onechange.sample(t(this.y), seg[1], y.var)$pval
            pass1 <- sample.pval<CHISQ.PVAL.THRESH
            pass2 <- abs(apply(as.matrix(this.y[1:seg[1],]),2,mean) - apply(as.matrix(this.y[(seg[1]+1):T,]),2,mean)) > MIN.SUFF.ABSDIFF
            cut1  <- abs(apply(as.matrix(this.y[1:seg[1],]),2,mean) - apply(as.matrix(this.y[(seg[1]+1):T,]),2,mean)) < MIN.REQ.ABSDIFF
            sampleswithsegment <- (pass1 | pass2) & (!cut1)
            this.yhat <- matrix(0,nrow=nrow(this.y),ncol=ncol(this.y))
            this.yhat[1:seg[1], sampleswithsegment] <- 
                matrix(rep(apply(as.matrix(this.y[1:seg[1],sampleswithsegment],nrow=seg[1]),2,mean), seg[1]),
                       nrow=seg[1], byrow=TRUE)
            this.yhat[(seg[1]+1):T,sampleswithsegment] <- 
                matrix(rep(apply(as.matrix(this.y[(seg[1]+1):T,sampleswithsegment],nrow=T-seg[1]),2,mean), T-seg[1]), 
                       nrow=T-seg[1], byrow=TRUE)        

        } else {
            ## Two change-points, square wave change, see who has change in middle
            sample.pval <- computeZ.squarewave.sample(t(this.y), seg, y.var)$pval
            pass1 <- sample.pval < CHISQ.PVAL.THRESH
            pass2 <- abs( apply(as.matrix(this.y[(seg[1] + 1):seg[2],]),2,mean) - 
                          apply(as.matrix(this.y[c(1:seg[1],(seg[2]+1):T),]),2,mean) ) > MIN.SUFF.ABSDIFF
            cut1  <- abs( apply(as.matrix(this.y[(seg[1] + 1):seg[2],]),2,mean) - 
                        apply(as.matrix(this.y[c(1:seg[1],(seg[2]+1):T),]),2,mean) ) < MIN.REQ.ABSDIFF
            sampleswithsegment <- (pass1 | pass2) & (!cut1)
            this.yhat <- matrix(0, nrow=nrow(this.y), ncol=ncol(this.y))
            this.yhat[(seg[1]+1):seg[2], sampleswithsegment] <- 
                matrix(rep(apply(as.matrix(this.y[(seg[1] + 1):seg[2], sampleswithsegment], nrow=seg[2]-seg[1]),2,mean), seg[2]-seg[1]), 
                       nrow=seg[2]-seg[1], byrow=TRUE)
            this.yhat[c(1:seg[1], (seg[2] + 1):T),sampleswithsegment] <- 
                matrix(rep(apply(as.matrix(this.y[c(1:seg[1], (seg[2] + 1):T), sampleswithsegment], nrow=T-(seg[2]-seg[1])),2,mean),T-(seg[2]-seg[1])),
                       nrow=T-(seg[2]-seg[1]),byrow=TRUE)
        }
    }
    
    list(carriers=sampleswithsegment, yhat=this.yhat)
}


`fmscbs` <- function(y, 
                     win, 
                     f = c(0.01, 0.1, 1), 
                     MIN.SNPs = 3, 
                     ALPHA = 0, 
                     GLOBAL.PVAL.CUTOFF = 0.0001, 
                     MAX.CHPTS = NA, 
                     WCHISQ.CUTOFF = NA)
{

## -------------------------------------------------------------
##     win:  maximum length of CNV (can be set to nrow(y))
##     MIN.SNPs: minimum length of CNV (can be set to 1).
##     GLOBAL.PVAL.CUTOFF:  Stops segmenting if pval falls below this.
##     MAX.CHPTS:  Stops segmenting after a maximum number of 
##                 change-points is reached.
##     plots:  If true, draws heatmap after every split.  
##             If your data set size not too large, this is
##             useful for tracking progress.  But plotting takes a long time
##             and the program runs a lot faster when plots=FALSE.
##
##   If your data is not properly normalized, or if you want a really sparse
##   summary (i.e. ignore everything except for the most significant changes),
##   it may be a good idea to use MAX.CHPTS to limit the number of recursions.
##
## --------------------------------------------------------------

    N <- ncol(y)
    T <- nrow(y)
    DELTA <- MIN.SNPs/T
  
    if(is.na(MAX.CHPTS)) MAX.CHPTS <- floor(T/MIN.SNPs)
  
    if(is.na(WCHISQ.CUTOFF)) {
        WCHISQ.CUTOFF <- getCutoffMultisampleWeightedChisq(pval = GLOBAL.PVAL.CUTOFF, 
                                                           m = T, 
                                                           delta = DELTA, 
                                                           win = win, 
                                                           N = N, 
                                                           alpha = ALPHA)
        cat("MSCBS: weighted chisquare cutoff = ", WCHISQ.CUTOFF, "\n")
    }

    ## initialize
    y.var <- compute.var(y)
    yhat  <- matrix(rep(apply(y,2,mean),T), nrow=T, byrow=TRUE)
    y.r   <- y - yhat
  
    chpts <- c(0, T)
    bestZ <- fcompute.max.Z(y.r, win, y.var, ALPHA, MIN.SNPs)
    best.subchpt <- matrix(bestZ$bestchpt, ncol=1, nrow=2)
    best.Z       <- bestZ$bestZ

    splitnum  <- 0
    chpt.hist <- vector("list", MAX.CHPTS)
        
    while (TRUE) {
        max.Z      <- max(best.Z, na.rm=TRUE)
        max.region <- which.max(best.Z)
    
        if (max.Z < WCHISQ.CUTOFF) {
            cat("Maximum Z-score is ", max.Z, 
                ", which does not exceed cutoff of ", WCHISQ.CUTOFF, 
                ". Segmentation finished.\n", sep="")
            break
        }

        if (length(chpts) > MAX.CHPTS + 2) {
            cat("Maximum number of change-points reached. Segmentation finished.\n")
            break
        }

        if (length(max.region)==0) {
            cat("Optimal region has no valid change-points. Segmentation finished.\n")      
            break
        }

        if (is.na(best.subchpt[1, max.region])) {
            cat("Optimal region has no valid change-points. Segmentation finished.\n")      
            break
        }
        
        splitnum <- splitnum + 1
        newchpt  <- c(best.subchpt[1, max.region], best.subchpt[2, max.region])

        ## Classify samples and update yhat
        y.r.classify <- mscbs.classify(y.r[(chpts[max.region]+1):chpts[(max.region+1)],], 
                                       newchpt - chpts[max.region], y.var, CHISQ.PVAL.THRESH=0.001)
        y.r.hat <- matrix(0, nrow=T, ncol=N)
        y.r.hat[(chpts[max.region]+1):chpts[(max.region+1)],] <- y.r.classify$yhat  ## add +1 to chpts[max.region]
        yhat <- yhat + y.r.hat
        y.r <- y.r - y.r.hat    

        if (!is.na(newchpt[2])) { ## The added change consistes of two change-points
            cat("Split ", splitnum, ": ", newchpt[1], ", ",newchpt[2], ", Z-score = ", max.Z, ".\n", sep="")
            y.r.L <- y.r[(chpts[max.region]+1):newchpt[1],]
            y.r.M <- y.r[(newchpt[1]+1):newchpt[2],]
            y.r.R <- y.r[(newchpt[2]+1):chpts[max.region+1],]

            bestZ.L <- fcompute.max.Z(y.r.L, win, y.var, ALPHA, MIN.SNPs)
            bestZ.M <- fcompute.max.Z(y.r.M, win, y.var, ALPHA, MIN.SNPs)
            bestZ.R <- fcompute.max.Z(y.r.R, win, y.var, ALPHA, MIN.SNPs)
            best.Z.new <- c(bestZ.L$bestZ, bestZ.M$bestZ, bestZ.R$bestZ)
            best.subchpt.new <- cbind(bestZ.L$bestchpt + chpts[max.region], 
                                      bestZ.M$bestchpt + newchpt[1], 
                                      bestZ.R$bestchpt + newchpt[2])
        } else { ## The added change-point is a singleton
            newchpt <- newchpt[1]
            cat("Split ", splitnum, ": ", newchpt, ", Z-score = ", max.Z, ".\n", sep="")
            y.r.L <- y.r[(chpts[max.region] + 1):newchpt,]
            y.r.R <- y.r[(newchpt + 1):chpts[max.region+1],]
      
            bestZ.L <- fcompute.max.Z(y.r.L, win,y.var, ALPHA, MIN.SNPs)
            bestZ.R <- fcompute.max.Z(y.r.R, win,y.var, ALPHA, MIN.SNPs)
            best.Z.new <- c(bestZ.L$bestZ, bestZ.R$bestZ)
            best.subchpt.new <- cbind(bestZ.L$bestchpt + chpts[max.region], 
                                      bestZ.R$bestchpt + newchpt)
        }
    
        if (max.region > 1) {
            leftpart   <- best.subchpt[,1:(max.region-1)]
            leftpart.Z <- best.Z[1:(max.region-1)]
        } else {
            leftpart   <- matrix(0, ncol=0, nrow=2)
            leftpart.Z <- matrix(0, ncol=0, nrow=0)
        }
    
        if (max.region+1 <= ncol(best.subchpt)) {
            rightpart   <- best.subchpt[,(max.region+1):ncol(best.subchpt)]
            rightpart.Z <- best.Z[(max.region+1):length(best.Z)]
        } else {
            rightpart   <- matrix(0, ncol=0, nrow=2)
            rightpart.Z <- matrix(0, ncol=0, nrow=0)
        }
        
        chpt.hist[[splitnum]] <- list(chpts = chpts, max.region = max.region,
                                      newchpt = newchpt, max.Z = max.Z, 
                                      carriers = y.r.classify$carriers)
        best.Z <- c(leftpart.Z, best.Z.new, rightpart.Z)    
        best.subchpt <- cbind(leftpart, best.subchpt.new, rightpart)
        chpts <- c(chpts[1:max.region], newchpt, chpts[(max.region+1):length(chpts)])
        
    } ## while loop

    chpt.hist <- chpt.hist[1:splitnum]
    if (length(chpts)>2) { 
        for (i in 2:length(chpts)) {
            yhat[(chpts[i-1]+1):chpts[i],] <- matrix(nrow = chpts[i]-chpts[i-1], ncol=ncol(yhat), 
                                                     data = apply(as.matrix(y[(chpts[i-1]+1):chpts[i],]),2,mean), byrow=TRUE)    
        }    
    } else {
        chpts <- c(0, T)
    }

    feature.matrix  <- yhat[chpts[-1], ]
    feature.regions <- data.frame(start = chpts[-length(chpts)] + 1, end = chpts[-1], row.names=NULL)

    list(chpt.hist=chpt.hist, chpts=chpts, yhat=yhat, 
         feature.matrix=feature.matrix, feature.regions=feature.regions)
}
