xsolve.pwl <- function(S,lambda,gprob,tmax,qmax,nstep,type,
                           alpha,salval,maxFac,nverb) {
#
# Function xsolve.pwl to solve numerically the system of d.e.'s
# for the value v_q(t) of a stock of q items at time t, using
# the method of Runge-Kutta, in the setting in which prices vary
# continuously but the price sensitivity functions are powers of a
# function (the "group size = 1" function) which is piecewise linear
# in x.  The optimal prices are determined by maximizing over the
# discrete set consisting of:
#
# 	(1) The ``knots'' x_1, ..., x_K of the price sensitivity function.
#	(2) The zeros of the derivative of the conditional (given an arrival
#           at time "t") expected value of the stock.  Note that this
#           conditional expected value is a piecewise polynomial in the
#           price "x" whence so is its derivative.
#
# We are solving the vector system of differential equations
#
#	vdot(t) = {script F}(v,t)
#
# The function "script F" is represented as scrF() in the code.

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
    
# Dig out the upper bound for the values of residual time.
    if(is.null(tmax)) {
        tmax <- attr(S,"tmax")
    } else if(tmax > attr(S,"tmax")) {
        stop(paste("Argument \"tmax\" is greater than the \"tmax\" attribute\n",
                   "of the pwl price sensitivity function specified as\n",
                   "argument \"S\".\n"))
    }
    
# Stow necessary objects in the environment of scrF.
    environment(scrF) <- new.env()
    assign("type",type,envir=environment(scrF))
    assign("lambda",lambda,envir=environment(scrF))
    assign("alpha",get("alpha",envir=environment(S)),envir=environment(scrF))
    assign("beta",get("beta",envir=environment(S)),envir=environment(scrF))
    assign("kn",get("kn",envir=environment(S)),envir=environment(scrF))
#
# There should be only one price sensitivity function. Create a new
# function equal to the original function raised to the power "n", with
# "n" a function argument; call it "dS" to make notation compatible
# with the "smooth" case.
    dS <- with(list(S=S),function(x,t,n) {S(x,t)^n})

# Stow necessary objects in the environment of cev.
    environment(cev) <- new.env()
    assign("dS",dS,envir=environment(cev))
    assign("gpr",gpr,envir=environment(cev))
    assign("alpha",alpha,envir=environment(cev))
    assign("maxFac",maxFac,envir=environment(cev))
    
# Do the Runge-Kutta thing.
    delta <- tmax/nstep
    tstor <- list()
    tvec  <- seq(0,tmax,length=nstep+1)
    tt    <- tvec[1]
    v     <- (1:qmax)*salval
    vdot  <- scrF(v,tt)
    x     <- scrF(v,tt,op=TRUE)
    tstor[[1]] <- list(x=x,v=v,vdot=vdot)
    
    for(i in 1:nstep) {
            assign("xback",x,envir=environment(cev))
    	m1 <- delta*scrF(v,tt)
    	m2 <- delta*scrF(v+m1/2,tt+delta/2)
    	m3 <- delta*scrF(v+m2/2,tt+delta/2)
    	m4 <- delta*scrF(v+m3,tt+delta)
    	v    <- v + (m1+2*m2+2*m3+m4)/6
    	tt   <- tvec[i+1]
    	vdot <- scrF(v,tt)
    	x    <- scrF(v,tt,op=TRUE)
    	tstor[[i+1]] <- list(x=x,v=v,vdot=vdot)
    	if(nverb > 0 & i%%nverb == 0) cat(i,"")
    	if(nverb > 0 & i%%(10*nverb) == 0) cat("\n")
    }
    if(nverb >= 0 & nstep%%10 != 0) cat("\n")
    putAway(tstor,type,jmax,qmax,tmax,discrete=FALSE)
}
