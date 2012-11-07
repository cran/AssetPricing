cev <- function(x,t,v,type,maximize=FALSE) {
#
# Conditional expected value (given an arrival at
# the time in question).  Note that objects dS, gpr, alpha
# and maxFac are assigned in the environment of cev.
#
if(!(type%in%c("sip","dip")))
	stop(paste("Type",type,"not recognized.\n"))
qmax <- length(v)
jmax <- length(gpr)
R    <- numeric(qmax)
xall <- if(is.numeric(x)) x else unlist(x) 
if(is.list(dS)) {
	S <- lapply(dS,function(f,x,t){f(x,t)},x=xall,t=t)
} else {
	S <- lapply(1:jmax,function(j,dS,x,t){dS(x,t,j)},dS=dS,x=xall,t=t)
}
S <- matrix(unlist(S),nrow=length(xall),ncol=jmax)

if(maximize) {
        E      <- parent.env(environment())
        xback  <- E$xback
        usexb  <- !is.null(xback)
	lopt   <- if(type=="sip") qmax else jmax*(qmax+0.5-jmax/2)
	xopt   <- numeric(lopt)
	istart <- 0
	for(q in 1:qmax) {
		vq   <- v[q]
		jtop <- min(q,jmax)
		rv   <- (c(v[q:1],0)[-1])[1:jtop]
		K    <- gpr[1:jtop]
		if(jtop < jmax) K[jtop] <- K[jtop] + alpha*(1-sum(K))
		xc <- if(is.list(x)) x[[q]] else x
		nx  <- length(xc)
		ind <- istart + (1:nx)
		tmp <- numeric(nx)
		if(type=="sip") {
			for(i in 1:nx) {
				ii <- ind[i]
				tmp[i] <- xc[i]*sum((1:jtop)*S[ii,1:jtop]*K) +
                        		     sum((rv - vq)*S[ii,1:jtop]*K)
			}
# New stuff.  If there are multiple maxima, take the one closest
# to the previous value.  (To prevent wild oscillation of traces.)
# Before, we just took the first of the possibilities.
                        if(usexb) {
                            eCrit <- diff(range(tmp))/maxFac
			    imax  <- rwmax(tmp,eCrit)
                        } else imax <- which.max(tmp)[1]
                        if(length(imax) > 1) {
                               xb   <- xback[q]
                               id   <- which.min(abs(xb-xc[imax]))[1]
                               imax <- imax[id]
                        }
# End new stuff.
			R[q]    <- tmp[imax] + vq
			xopt[q] <- xc[imax]
		} else { # type == "dip"
			Rtmp <- numeric(jtop)
			xtmp <- numeric(jtop)
			for(j in 1:jtop) {
				for(i in 1:nx) {
					ii <- ind[i]
					tmp[i] <- (j*xc[i] + rv[j] - vq)*S[ii,j]*K[j]
				}
# New stuff.  If there are multiple maxima, take the one closest
# to the previous value.  (To prevent wild oscillation of traces.)
# Before, we just took the first of the possibilities.
			        imax <- rwmax(tmp)
                                if(length(imax) > 1 & usexb) {
                                       iqj  <- qj2i(q,j,qmax)
                                       xb   <- xback[iqj]
                                       id   <- which.min(abs(xb-xc[imax]))[1]
                                       imax <- imax[id]
                                } else imax <- imax[1]
# End new stuff.
				Rtmp[j] <- tmp[imax]
				xtmp[j] <- xc[imax]
			}
			R[q] <- sum(Rtmp) + vq
			iv   <- qj2i(q,1:jtop,qmax)
			xopt[iv] <- xtmp
		}
	istart <- if(is.list(x)) istart + nx else 0
	}
	attr(R,"xopt") <- xopt
} else {
	for(q in 1:qmax) {
		vq   <- v[q]
		jtop <- min(q,jmax)
		rv   <- (c(v[q:1],0)[-1])[1:jtop]
		K    <- gpr[1:jtop]
		if(jtop < jmax) K[jtop] <- K[jtop] + alpha*(1-sum(K))
		if(type=="sip") {
			R[q] <- x[q]*sum((1:jtop)*S[q,1:jtop]*K) +
                        	     sum((rv - vq)*S[q,1:jtop]*K) + vq
		} else {
			i    <- qj2i(q,1:jtop,qmax)
			R[q] <- sum(((1:jtop)*x[i] + rv - vq)*
                        	            S[cbind(i,1:jtop)]*K) + vq
		}
	}
}
R
}
