xsolve <- function(S,lambda,gprob=1,tmax=NULL,qmax,prices=NULL,nstep,type="sip",
                   alpha=NULL,salval=0,maxFac=1000,nverb=0) {
#
# Dispatch the problem to one of xsolve.cont(), xsolve.disc(),
# xsolve.pwl(), depending on the nature of the price sensitivity
# function "S".
#
    soltype <- findSolType(S,prices)

    if(is.numeric(lambda)) {
        if(length(lambda) != 1 || lambda <=0)
            stop("When \"lambda\" is numeric it must be a positive scalar.\n")
        lambda <- with(list(lambda=lambda),function(t){rep(lambda,length(t))})
    }

    switch(EXPR = soltype,
	cont = xsolve.cont(S,lambda,gprob,tmax,qmax,nstep,type,
                           alpha,salval,nverb),
	disc = xsolve.disc(S,lambda,gprob,tmax,qmax,prices,nstep,type,
                           alpha,salval,maxFac,nverb),
	pwl  = xsolve.pwl(S,lambda,gprob,tmax,qmax,nstep,type,
                           alpha,salval,maxFac,nverb)
    )
}
