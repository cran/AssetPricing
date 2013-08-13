xsolve <- function(S,lambda,gprob=1,tmax=NULL,qmax,prices=NULL,nout=300,type="sip",
                   alpha=NULL,salval=0,maxFac=1000,method="lsoda") {
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
	cont = xsolve.cont(S,lambda,gprob,tmax,qmax,nout,type,
                           alpha,salval,method=method),
	disc = xsolve.disc(S,lambda,gprob,tmax,qmax,prices,nout,type,
                           alpha,salval,maxFac,method=method),
	pwl  = xsolve.pwl(S,lambda,gprob,tmax,qmax,nout,type,
                           alpha,salval,maxFac,method=method)
    )
}
