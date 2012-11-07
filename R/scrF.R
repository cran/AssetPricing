scrF  <- function(v,tt,op=FALSE) {
#
# Note:  The object lambda and type are (always) assigned in the
# environment of scrF.  When scrF is called by xsolve.disc(), the
# vector "x" of *possible* prices is assigned in the environment
# of scrF.  When scrF is called by vsolve() the pricing policy "x"
# is assigned in the environment of scrF.  When scrF is called by
# xsolve.pwl() the lists "alpha" and "beta" of coefficient functions
# and knot vector "kn" for the piecewise linear representation of
# S(x,t) are assigned in the environment of scrF.
#
# The value of this function is either vdot = ``script F''(v,t),
# the (unstated :-( ) vectorized version of equation (2) of the
# paper, *or* when op=TRUE the vector of optimal prices.  (The latter is
# applicable only in the case of discrete prices or of a piecewise
# linear price sensitivity function.)

E <- parent.env(environment())
if(is.null(E$x)) {
# Here scrF is being called by xsolve.pwl() and the price elasticity
# functions are piecewise linear (rather than discrete).  We want
# to maximize over the possible price values, but first we need to
# determine what the possible prices are.  Note that "x" is a list
# here, with the q-th entry of the list being the vector of possible
# optimal prices for stock size "q".
	qmax <- length(v)
	jmax <- length(gpr)
	x    <- vector("list",qmax)
	for(q in 1:qmax) {
		jtop <- min(q,jmax)
		Kpa  <- gpr[1:jtop]
		if(jtop < jmax)
			Kpa[jtop] <- Kpa[jtop] + environment(cev)$alpha*(1-sum(Kpa))
		x[[q]] <- getPossPrices(v[1:q],tt,E$alpha,E$beta,E$kn,Kpa)
	}
} else if(inherits(x,"flap")) {
# Here scrF is being called by vsolve() --- pricing policy is given
# so x supplies the actual prices. No optimization to be done,
# so just return the expected values.
	xx <- sapply(x,function(f,t){f(t)},t=tt)
	return(lambda(tt)*(-v + cev(xx,tt,v,type)))
}

# At this point either scrF was called by xsolve.pwl() and the
# vector of possible prices x has been constructed, or scrF is being
# called by xsolve.disc() and the set of possible (discrete) prices
# was provided in the vector x which was assigned in the environment
# of this function.  In either case the vector of possible prices
# is now available and we can maximize over this vector.
R <- cev(x,tt,v,type,maximize=TRUE)

# The optimum prices have been determined (both discrete and
# pwl settings).  If "op" is true we want the optimum prices;
# otherwise we want the vdot values *at* the optimum prices.
if(op) attr(R,"xopt") else lambda(tt)*(R - v)
}
