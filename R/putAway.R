putAway <- function(tstor,type,jmax,qmax,tmax,discrete,prices=NULL) {
#
# Equation solved; now do a bunch of housekeeping to store
# the results in a convenient manner.

# Put away the prices.
#
# If putAway() is being called by vsolve() then these prices are
# *not* necessarily optimal and tstor$x is already an object of
# class "flap".  In this case "x" is extracted (nothing needs
# to be done to it) and the other component of tstor (also called
# "tstor") is extracted so that the "v" and "vdot" dot values can
# be processed "in the usual way".  Otherwise the "x" values are
# processed (from tstor) in the usual way.
#
# Note that if type == "dip" then the entry x[[i]] is equal to
# the function x_qj(t) where i = (j-1)*(qmax - j/2) + q.
#
# If type == "sip" then at this stage jmax is effectively 1, whatever
# its `real' value.  At this stage jmax relates only to the indexing
# of the entries of the prices "x".  (Previously the value of jmax had
# an impact on the calculation of expected values of stocks.)  If type
# is equal to "sip" then the index "i" of the prices simply runs from
# 1 to qmax just as it would if jmax were equal to 1.
if(type=="sip") jmax <- 1

if(inherits(tstor$x,"flap")) {
    x     <- tstor$x
    tstor <- tstor$tstor
    nsp1  <- length(tstor)
    tvec  <- seq(0,tmax,length=nsp1)
} else {
    nsp1 <- length(tstor)
    tvec <- seq(0,tmax,length=nsp1)
    nc   <- if(type=="dip") jmax*(qmax - (jmax-1)/2) else qmax
    xx   <- as.data.frame(matrix(unlist(lapply(tstor,function(z){z$x})),
                           ncol=nc,byrow=TRUE))
    if(discrete) {
	x <- lapply(xx,function(x,t){approxfun(x=t,y=x,method="constant",
                                yleft=NA,yright=NA)},t=tvec)
	x <- lapply(x,a2sf,tt=tvec)
    } else x <- lapply(xx,function(x,t){splinefun(x=t,y=x)},t=tvec)

    ylim <- range(xx)
    attr(x,'tlim') <- c(0,tmax)
    attr(x,'ylim') <- ylim
    attr(x,'qmax') <- qmax
    attr(x,'jmax') <- jmax
    if(!is.null(prices)) attr(x,'prices') <- prices
    comment(x) <- "Optimal prices."
    class(x) <- "flap"
    if(discrete) class(x) <- c("flap","pwc.flap")  # pwc <--> piecewise constant.
    if(type=="dip") class(x) <- c(class(x),"di.flap") # di  <--> doubly indexed.
}

# Put away the optimal expected values.
vv   <- as.data.frame(matrix(unlist(lapply(tstor,function(z){z$v})),
                             ncol=qmax,byrow=TRUE))
v    <- lapply(vv,function(x,t){splinefun(x=t,y=x)},t=tvec)
ylim <- range(vv)
attr(v,'tlim') <- c(0,tmax)
attr(v,'ylim') <- ylim
attr(v,'qmax') <- qmax
attr(v,'jmax') <- 1
class(v) <- "flap"

# Put away the derivatives of the optimal expected values.
vd   <- as.data.frame(matrix(unlist(lapply(tstor,function(z){z$vdot})),
                             ncol=qmax,byrow=TRUE))
vdot <- lapply(vd,function(x,t){splinefun(x=t,y=x)},t=tvec)
ylim <- range(vd)
attr(vdot,'tlim') <- c(0,tmax)
attr(vdot,'ylim') <- ylim
attr(vdot,'qmax') <- qmax
attr(vdot,'jmax') <- 1
class(vdot) <- "flap"
rslt <- list(x=x,v=v,vdot=vdot)
rslt
}
