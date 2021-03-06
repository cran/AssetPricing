\name{AssetPricing-internal}
\alias{a2sf}
\alias{cev}
\alias{findSolType}
\alias{getPossPrices}
\alias{getS}
\alias{i2qj}
\alias{initx}
\alias{.newt}
\alias{plot.flap}
\alias{plot.stepfun}
\alias{progRep}
\alias{putAway}
\alias{qj2i}
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
getPossPrices(v,t,alpha,beta,kn,Kpa,type)
getS(dS,xc,t,jays)
i2qj(i,qmax,jmax)
initx(v,type)
.newt(fn,start,...,eps.p = 1e-08,eps.v = NULL,
     maxit = 50,verb = FALSE)
\method{plot}{flap}(x,xlim=NULL,ylim=NULL,lty=NULL,cols=NULL,xlab=NULL,
                    ylab=NULL,main=NULL,main.panel= NULL,groups=NULL,
                    add=FALSE,gloss=FALSE,glind=NULL,extend=0.3,col.gloss=1,
                    cex.gloss=0.8,mfrow=NULL,...)
\method{plot}{stepfun}(x, xval, xlim, ylim = range(c(y, Fn.kn)), xlab = "x",
                       ylab = "f(x)", main = NULL, add = FALSE, verticals = TRUE,
                       do.points = TRUE, pch = par("pch"), col = par("col"),
                       col.points = col, cex.points = par("cex"), col.hor = col,
                       col.vert = col, lty = par("lty"), lwd = par("lwd"), ...)
progRep(info,verbInt,tt,tmax)
putAway(odeRslt,type,jmax,qmax,soltype,x=NULL,prices=NULL)
qj2i(q,j,qmax)
scrG(tt,x,parms,verbInt,tmax,info,...)
scrF(tt,v,parms,verbInt,tmax,info,...)
turnPts(a,b,v,Kpa,xlo,xhi,type)
vupdate(x,tt,type)
xsolve.cont(S,lambda,gprob,tmax,qmax,nout,type,
                           alpha,salval,method,verbInt)
xsolve.disc(S,lambda,gprob,tmax,qmax,prices,nout,type,
                           alpha,salval,epsilon,method,verbInt)
xsolve.pwl(S,lambda,gprob,tmax,qmax,nout,type,
                           alpha,salval,epsilon,method,verbInt)
}
\details{
  These functions are auxiliary and are not intended to be called by
  the user.
% The function ".newt" is a copy of "newt" from my personal "misc"
% package to which yer average user will not have access.
% The function plot.stepfun() is a slight modification of the
% function of the same name from the package "stats".  It is modified
% so that the "xlim" argument behaves properly (like unto the way
% that argument behaves for plot.function().
% On the advice of Kurt Hornik (31/03/2018) I created a new *generic*
% plot function, with default method equal to graphics::plot(), so
% as to properly accommodate the modified plot.stepfun() method.
%
}
\keyword{internal}
