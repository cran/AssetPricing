\name{AssetPricing-internal}
\alias{cev}
\alias{i2qj}
\alias{initx}
\alias{.newt}
\alias{qj2i}
\alias{scrF}
\alias{scrG}
\alias{vupdate}
\title{Internal AssetPricing functions}
\description{
       Internal AssetPricing functions.
}
\usage{
cev(x,t,v,type)
i2qj(i,qmax,jmax=qmax)
initx(v,type)
.newt(fn,start,...,eps.p = 1e-08,eps.v = NULL,
     maxit = 50,verb = FALSE)
qj2i(q,j,qmax)
scrG(x,tt)
scrF(v,tt)
vupdate(x,tt,type)
}
\details{
  These functions are auxilliary and are not intended to be called by
  the user.
% The function ".newt" is a copy of "newt" from my personal "misc"
% package to which yer average user will not have access.
}
\keyword{internal}
