\name{AssetPricing-internal}
\alias{a2sf}
\alias{cev}
\alias{findSolType}
\alias{getPossPrices}
\alias{i2qj}
\alias{initx}
\alias{.newt}
\alias{plot.stepfun}
\alias{putAway}
\alias{qj2i}
\alias{rwmax}
\alias{scrF}
\alias{scrG}
\alias{turnPts}
\alias{vupdate}
\alias{xsolve.cont}
\alias{xsolve.disc}
\alias{xsolve.pwl}
\title{Internal AssetPricing functions}
\description{
       Internal AssetPricing functions.
}
\usage{
a2sf(x,tt,xlim,npts=1000)
cev(x,t,v,type,maximize=FALSE)
findSolType(S,prices)
getPossPrices(v,t,alpha,beta,kn,Kpa)
i2qj(i,qmax,jmax)
initx(v,type)
.newt(fn,start,...,eps.p = 1e-08,eps.v = NULL,
     maxit = 50,verb = FALSE)
plot.stepfun(x, xval, xlim, ylim = range(c(y, Fn.kn)), xlab = "x",
             ylab = "f(x)", main = NULL, add = FALSE, verticals = TRUE,
             do.points = TRUE, pch = par("pch"), col = par("col"),
             col.points = col, cex.points = par("cex"), col.hor = col,
             col.vert = col, lty = par("lty"), lwd = par("lwd"), ...)
putAway(tstor,type,jmax,qmax,tmax,discrete,prices=NULL)
qj2i(q,j,qmax)
rwmax(x,eps)
scrG(x,tt)
scrF(v,tt,op=FALSE)
turnPts(a,b,v,Kpa,xlo,xhi)
vupdate(x,tt,type)
xsolve.cont(S,lambda,gprob,tmax,qmax,nstep,type,
                           alpha,salval,nverb)
xsolve.disc(S,lambda,gprob,tmax,qmax,prices,nstep,type,
                           alpha,salval,maxFac,nverb)
xsolve.pwl(S,lambda,gprob,tmax,qmax,nstep,type,
                           alpha,salval,maxFac,nverb)
}
\details{
  These functions are auxilliary and are not intended to be called by
  the user.
% The function ".newt" is a copy of "newt" from my personal "misc"
% package to which yer average user will not have access.
% The function plot.stepfun() is a slight modification of the
% function of the same name from the package "stats".  It is modified
% so that the "xlim" argument behaves properly (like unto the way
% that argument behaves for plot.function().
}
\keyword{internal}
