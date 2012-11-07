xsolve.disc <- function(S,lambda,gprob,tmax,qmax,prices,nstep,type,
                        alpha,salval,maxFac,nverb) {
#
# Function xsolve.disc to solve numerically the system of d.e.'s
# for the value v_q(t) of a stock of q items at time t, using the
# method of Runge-Kutta, in the setting in which prices vary over a
# ***discrete*** set.  The optimal prices are then determined by
# maximizing over this discrete set.  (Note that time is thought
# of as ***decreasing*** toward the ``departure time'' of 0.)
# But we solve ``forward in time'', starting from 0.
#
# We are solving the vector system of differential equations
#
#	vdot(t) = {script F}(v,t)
#
# In what follows, ``script F'' is denoted by ``scrF()''.
# In addition to "v" and "t", scrF() has an additional auxilliary
# argument "op".  This is a logical scalar which defaults to FALSE;
# if "op" is TRUE then the value returned by scrF() has an attribute
# "xopt" giving the vector of values of the optimal prices (chosen
# from amongst the finite set of possible prices.
#

# Make sure the group size probabilities are OK.
    gpr <- if(is.function(gprob)) gprob(1:qmax) else gprob
    if(!is.numeric(gpr)) stop("Group size probabilities are not numeric.\n")
    if(any(gpr<0) | sum(gpr) > 1)
    	stop("Group size probabilities are not probabilities!\n")
    jmax <- max(which(gpr > sqrt(.Machine$double.eps)))
    jmax <- min(jmax,qmax)
    gpr  <- gpr[1:jmax]
    if(is.null(alpha)) {
       if(jmax > 1) stop(paste("When the maximum group size is great than 1,\n",
                               "\"alpha\" must be specified.\n"))
       alpha <- 1
    }

# If jmax = 1 we might as well set type equal to "sip" --- since
# indexing according to group size is "degenerate" in this case.
    if(jmax==1) type <- "sip"

# Make sure tmax is specified.
    if(is.null(tmax)) stop("Argument \"tmax\" was not specified.\n")

# If there is only one price sensitivity function, create a new
# function equal to the original function raised to the power "n", with
# "n" a function argument.

    if(is.list(S)) {
    	if(length(S) < jmax)
    		stop(paste("Length of \"S\" as a list must be at least ",
                               "\"jmax\" = ",jmax,".\n",sep=""))
    	if(!all(sapply(S,is.function)))
    		stop("At least one entry of \"S\" is NOT a function.\n")
    	} else if(is.function(S)) {
    	oldS <- S
    	S <- function(x,t,n) {oldS(x,t)^n}
    	environment(S) <- new.env()
    	assign("oldS",oldS,envir=environment(S))
    } else {
    	stop("Argument \"S\" must be either a function or a list of functions.\n")
    }

    environment(scrF)    <- new.env()
    environment(cev)     <- new.env()
#
    assign("dS",S,envir=environment(cev))  # "dS" to make notation compatible with
                                           # the "smooth" case.
    assign("gpr",gpr,envir=environment(cev))
    assign("alpha",alpha,envir=environment(cev))
    assign("maxFac",maxFac,envir=environment(cev))
#
    assign("x",prices,envir=environment(scrF))
    assign("lambda",lambda,envir=environment(scrF))
    assign("type",type,envir=environment(scrF))

# Do the Runge-Kutta thing.
    tvec  <- seq(0,tmax,length=nstep+1)
    tt    <- tvec[1]
    delta <- tmax/nstep
    v     <- (1:qmax)*salval
    tt    <- tvec[1]
    vdot  <- scrF(v,tt)
    x     <- scrF(v,tt,op=TRUE)
    tstor   <- list()
    tstor[[1]] <- list(x=x,v=v,vdot=vdot)

    for(i in 1:nstep) {
            assign("xback",x,envir=environment(cev))
    	m1 <- delta*scrF(v,tt)
    	m2 <- delta*scrF(v+m1/2,tt+delta/2)
    	m3 <- delta*scrF(v+m2/2,tt+delta/2)
    	m4 <- delta*scrF(v+m3,tt+delta)
    	v  <- v + (m1+2*m2+2*m3+m4)/6
    	tt <- tvec[i+1]
    	vdot <- scrF(v,tt)
    	x    <- scrF(v,tt,op=TRUE)
    	tstor[[i+1]] <- list(x=x,v=v,vdot=vdot)
    	if(nverb > 0 & i%%nverb == 0) cat(i,"")
    	if(nverb > 0 & i%%(10*nverb) == 0) cat("\n")
    }
    if(nverb >= 0 & nstep%%10 != 0) cat("\n")
    putAway(tstor,type,jmax,qmax,tmax,discrete=TRUE,prices=prices)
}
