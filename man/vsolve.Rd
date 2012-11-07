\name{vsolve}
\alias{vsolve}
\title{
  Solve for expected value of assets.
}
\description{
  Solves a system of coupled differential equations for the expected
  value of a number \code{q} of (\dQuote{perishable}) assets, with
  \code{q} running from \code{1} to \code{qmax}, \emph{given} a pricing
  policy.  Treats the system in a vectorized form and uses the method
  of Runge-Kutta.
}
\usage{
    vsolve(S, lambda, gprob, tmax=NULL, x, nstep,
                   alpha=NULL, salval=0, nverb=0)
}
\arguments{
  \item{S}{
  An expression, or list of expressions, or a function or list
  of functions, specifying the price sensitivity functions
  \code{S_j(x,t)}.  See \bold{Details}.
}
  \item{lambda}{
  A function (of residual time \code{t} --- see \code{tmax}) or
  a positive constant specifying the intensity of the (generally
  inhomogeneous) Poisson process of arrival times of groups of
  potential customers.
}
  \item{gprob}{
  A function (to calculate probabilities) or a numeric vector of 
  probabilities determining the distribution of the size of
  an arriving group of customers.  Must be compatible with
  certain characteristics of \code{s} (see below).  See
  \bold{details}.
}
  \item{tmax}{
  The maximum residual time; think of this as being the initial
  time at which the assets go on sale (with time \emph{decreasing}
  to zero, at which point the value of each asset drops to the
  \dQuote{salvage value} (\code{salval}), usually \code{0}).
  The system of differential equations is solved over the time
  interval [0,\code{tmax}].  See \bold{Details}.
}
  \item{x}{
  An object of class \code{flap} (see \code{\link{xsolve}()})
  specifying the (given) pricing policy.  It has the form of
  a list of functions \code{x_i(t)}, with \code{i} running
  from 1 to \code{qmax} if \code{x} is \dQuote{singly indexed},
  i.e. \emph{not} of class \code{di.flap} or \code{x_ij(t)}, with
  \code{i} running from 1 to \code{qmax} and \code{j} running from 1
  to the maximum group size if \code{x} is \dQuote{doubly indexed},
  i.e.  \emph{is} of class \code{di.flap}.  Note that \code{x}
  has (must have) an \emph{attribute} \code{qmax} specifying the
  maximum number of assets available for sale, i.e. the number of
  assets available at the starting (residual) time \code{tmax}.
}
  \item{nstep}{
  The number of (equal) sub-intervals into which the interval
  [0,\code{tmax}] is divided in the process of solving the system
  of differential equations numerically.
}
  \item{alpha}{
  A numeric scalar between 0 and 1 specifying the probability that an
  arriving group of size \code{j > q} (where \code{q} is the number
  of assets remaining for sale) will \emph{consider} purchasing
  (all of) these remaining assets.  It is irrelevant (and defaults
  to 1 as a \dQuote{place holder}) if customers always arrive singly.
}
  \item{salval}{
  A (non-negative) numeric scalar specifying the \dQuote{salvage
  value} of an asset --- i.e. the quantity to which the value of an
  asset drops at residual time \code{t=0}.  Usually \code{salval}
  is equal to 0.
}
  \item{nverb}{
  Non-negative integer; if \code{nverb > 0} then a \dQuote{progress
  report} (actually consisting solely of the step number) is
  printed out every \code{nverb} steps of the Runge-Kutta procedure.
  (If \code{nverb==0} the procedure runs silently.)
}
}

\details{
  The components of the argument \code{S} may be provided either
  as expressions or functions.  If the former, these expressions
  should be amenble to differentiation with respect to \code{x}
  and \code{t} via the function \code{deriv3()}.  This is
  essentially a matter of convenience; the derivatives are not
  actually used by \code{vsolve}.  The expressions are turned into
  functions by \code{deriv3()} in the same manner as is used by
  \code{\link{xsolve}()}.  See the help for \code{\link{xsolve}()}
  for further information about the required nature of \code{S}.

  The argument \code{tmax} (if specified) must be less than or equal
  to the \code{tmax} attribute of argument \code{S} if \code{S} is a
  piecewise linear price sensitivity function, and must also be less
  than or equal to the \code{tlim} attribute of argument \code{x}.

  If \code{tmax} is not specified it will be set equal to the
  \code{tmax} attribute of argument \code{S} if \code{S} is a
  piecewise linear price sensitivity function, in which case
  this attribute must be less than or equal to the \code{tlim}
  attribute of argument \code{x}.  (If this is not so then \code{S}
  and \code{x} are incompatible.)  Otherwise \code{tmax} will be
  set equal to the \code{tlim} attribute of argument \code{x}.

  The argument \code{gprob} determines the range of possible values
  of the size of an arriving group of customers.  The maximum value
  of this group size is in effect that value of \code{j} for which
  the corresponding probability value is numerically distinguishable
  from zero.  If the argument \code{x} is a \dQuote{doubly indexed}
  list of functions (was created with \code{type="dip"}) then
  the maximum value of group size as determined by \code{gprob}
  must be compatible with the indexing scheme of \code{x}.
  That is to say, it must be less than or equal to the \code{jmax}
  attribute of \code{x}, otherwise an error is given.  Note that
  if single indexing is in effect (i.e. \code{x} was created with
  \code{type="sip"}) then this attribute is equal to 1, but for
  single indexing \code{x} does not depend on group size and so no
  restriction is imposed.
}
\value{
A list with components
\item{x}{The argument \code{x} which was passed to \code{vsolve},
   possibly with its \code{tlim} attribute modified.  It is
   an object of class \code{flap}.}
\item{v}{An object of class \code{flap} whose entries are
   (spline) functions \code{v_q(t)} specifying the expected value of
   \code{q} assets at time \code{t} as determined by numerically
   solving the coupled system of differential equations.
}
\item{vdot}{An object of class \code{flap} whose entries are
   the derivatives (with respect to \code{t}) of the functions \code{v_q(t)}
   described above.  The values of these derivatives are determined as
   the left hand side of the differential equations being solved.
}
}
\references{
Baneree, P. K. and Turner, T. R.
A flexible model for the pricing of perishable assets.
Omega, vol. 40, number 5, 2012, pages 533 -- 540,
doi: 10.1016/j.omega.2011.10.001.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}
  \url{http://www.math.unb.ca/~rolf}
}
\seealso{
\code{\link{xsolve}()}, \code{\link{plot.flap}()}
}
\examples{
#
# In these examples "nstep" has been set equal to 100
# to reduce the time for package checking on CRAN.
# A more "realistic" value is nstep=1000.  Likewise "qmax"
# has been set equal to 5 which is an unrealistically low
# value for the total number of assets.
#
S <- expression(exp(-kappa*x/(1+gamma*exp(-beta*t))))
attr(S,"parvec") <- c(kappa=10/1.5,gamma=9,beta=1)
lambda1 <- function(tt){
#	eps <- sqrt(.Machine$double.eps)
	84*(1-tt)
}

# Optimal pricing policy assuming customers arrive singly:
X <- xsolve(S=S,lambda=lambda1,gprob=1,tmax=1,qmax=5,nstep=100,nverb=50)
lambda2 <- function(tt){
#	eps <- sqrt(.Machine$double.eps)
	36*(1-tt)
}
# Expected values if the customers actually arrive in groups, using the
# (sub-optimal) pricing policy based on the (erroneous) assumption that
# they arrive singly.  Note that the two scenarios are ``comparable'' in
# that the expected total number of customers is 42 in each case.
V <- vsolve(S=S,lambda=lambda2,gprob=(5:1)/15,x=X$x,nstep=100,
            alpha=0.5,nverb=50)
}
\keyword{ math }