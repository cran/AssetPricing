\name{vsolve}
\alias{vsolve}
\title{
  Solve for expected value of assets.
}
\description{
  Solves the system of coupled differential equations for the expected
  value of a number \code{q} of (\dQuote{perishable}) assets, with
  \code{q} running from \code{1} to \code{qmax}, \emph{given} a pricing
  policy.  Treats the system in a vectorized form and uses the method
  of Runge-Kutta.
}
\usage{
   vsolve(type, tmax, nstep, x, lambda, S, alpha, gprob,
          salval = 0, verb = FALSE)
}
\arguments{
  \item{type}{
  Scalar character string taking one of the two values \code{sip}
  (singly indexed prices) and \code{dip} (doubly indexed prices).
  In the \code{dip} case the price of the asset which is quoted
  to the arriving group is allowed to depend upon the group size
  (as well as upon [residual] time and the number of assets remaining
  for sale.  In the \code{sip} case the quoted price does not depend
  upon group size.
}
  \item{tmax}{
  The maximum residual time; think of this as being the initial time
  at which the assets go on sale (with time \emph{decreasing} to zero,
  at which point the value of each asset drops to the \dQuote{salvage
  value} (\code{salval}), usually \code{0}).  The system of differential
  equations is solved over the time interval [0,\code{tmax}].
}
  \item{nstep}{
  The number of (equal) sub-intervals into which the interval
  [0,\code{tmax}] is divided in the process of solving the system
  of differential equations numerically.
}
  \item{x}{
  An object of class \code{funlist} (see \code{\link{xsolve.cont}()})
  specifying the (given) pricing policy.  It has the form of a
  list of functions \code{x_i(t)}, with \code{i} running from 1 to
  \code{qmax} if \code{type=="sip"} or \code{x_ij(t)}, with \code{i}
  running from 1 to \code{qmax} and \code{j} running from 1 to the
  maximum group size if \code{type=="dip"}.  Note that \code{x} has
  (must have) an \emph{attribute} \code{qmax} specifying the maximum
  number of assets available for sale, i.e. the number of assets
  available at the starting (residual) time \code{tmax}.
}
  \item{lambda}{
  A function (of \code{t}) specifying the intensity of the inhomogeneous
  Poisson process of arrival times of groups of potential customers.
}
  \item{S}{
  An expression or a list of expressions specifying the price
  sensitivity functions \code{S_j(x,t)}.  If \code{S} is a single
  expression then \code{S_j(x,t)} is taken to be this expression
  raised to the power \code{j}.  If \code{S} is a list, then
  \code{S_j(x,t)} is taken to be its \code{j}-th entry.

  Note that in general the expression or expressions will depend upon
  a number of \emph{parameters} as well as upon \code{x} and \code{t}.
  The values of these parameters are specified via an attribute or
  attributes.  If \code{S} is a (single) expression it has (must have)
  an attribute called \code{parlist} which is a \emph{named} vector of
  parameter values.  If \code{S} is a list of expressions each entry
  of the list has (must have) such an attribute.
}
  \item{alpha}{
  A numeric scalar between 0 and 1 specifying the probability that
  an arriving group of size \code{j > q} (where \code{q} is the number
  of assets remaining for sale) will \emph{consider} purchasing (all of)
  these remaining assets.
}
  \item{gprob}{
  A function (to calculate probabilities) or a numeric vector of 
  of probabilities determining the distribution of the size of
  an arriving group of customers.
}
  \item{salval}{
  A (non-negative) numeric scalar specifying the \dQuote{salvage
  value} of an asset --- i.e. the quantity to which the value of
  an asset drops at residual time \code{t=0}.  Usually \code{salval}
  is equal to 0.
}
  \item{verb}{
  Logical scalar; if \code{TRUE} then a \dQuote{progress report}
  (actually consisting solely of the step number) is printed out
  every 50 steps of the Runge-Kutta procedure.
}
}
\value{
A list with components
\item{x}{The argument \code{x} which was passed to \code{vsolve};
   an object of class \code{funlist}.}
\item{v}{An object of class \code{funlist} whose entries are
   (spline) functions \code{v_q(t)} specifying the expected value of
   \code{q} assets at time \code{t} as determined by numerically
   solving the coupled system of differential equations.
}
\item{vdot}{An object of class \code{funlist} whose entries are
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
\code{\link{xsolve.cont}()}, \code{\link{plot.funlist}()}
}
\examples{
S <- expression(exp(-kappa*x/(1+gamma*exp(-beta*t))))
attr(S,"parlist") <- c(kappa=10/1.5,gamma=9,beta=1)
lambda1 <- function(tt){
	eps <- sqrt(.Machine$double.eps)
	if(tt < -eps | tt > 1+eps) stop("Time out of range.\n")
	84*(1-tt)
}

# Optimal pricing policy assuming customers arrive singly:
X <- xsolve.cont(type="sip",tmax=1,nstep=1000,qmax=5,lambda=lambda1,
                 S=S,alpha=0.5, gprob=1,verb=TRUE)
lambda2 <- function(tt){
	eps <- sqrt(.Machine$double.eps)
	if(tt < -eps | tt > 1+eps) stop("Time out of range.\n")
	36*(1-tt)
}
# Expected values if the customers actually arrive in groups, using the
# (sub-optimal) pricing policy based on the (erroneous) assumption that
# they arrive singly.  Note that the two scenarios are ``comparable'' in
# that the expected total number of customers is 42 in each case.
V <- vsolve(type="sip",tmax=1,nstep=1000,x=X$x,lambda=lambda2,
                 S=S,alpha=0.5,gprob=(5:1)/15,verb=TRUE)
}
\keyword{ math }
