\name{xsolve.cont}
\alias{xsolve.cont}
\title{
  Solve for an optimal continuous pricing policy
}
\description{
  Solves the system of coupled differential equations for the
  optimal prices of a number \code{q} of (\dQuote{perishable})
  assets, with \code{q} running from \code{1} to \code{qmax} in
  the case of \emph{continuous} prices.  Treats the system in a
  vectorized form and uses the method of Runge-Kutta.
}
\usage{
   xsolve.cont(type, tmax, nstep, qmax, lambda, S, alpha, gprob,
               salval = 0, verb = FALSE, n.out = max(500, round(nstep/200)))
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
  \item{qmax}{
  The maximum number of assets available for sale, i.e. the number
  of assets available at the starting (residual) time \code{tmax}.
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
  \item{n.out}{
  The number of points at which to evaluate a spline fitted to
  the numerical results; these points are those on which the
  spline \emph{functions} which form the entries of the output
  lists are based.
}
}
\details{
  No idea what the use of the \code{n.out} argument is, or why
  I wrote the code to depend on it.  I think I will probably
  get rid of it.
}
\value{
A list with components
\item{x}{The optimal pricing policy, chosen to maximize the expected
   value of the remaining assets at any given time; an object of
   class \code{funlist}.  It has the form of a list of (spline)
   functions \code{x_i(t)}, with \code{i} running from 1 to
   \code{qmax} if \code{type=="sip"} or \code{x_ij(t)}, with \code{i}
   running from 1 to \code{qmax} and \code{j} running from 1 to the
   maximum group size if \code{type=="dip"}.  Note that \code{x}
   has (must have) an \emph{attribute} \code{qmax} specifying the
   maximum number of assets available for sale, i.e. the number of
   assets available at the starting (residual) time \code{tmax}.
   Note that if \code{type == "dip"} then the entry \code{x[[i]]} is
   equal to the function \code{x_qj(t)} where \code{i = (j-1)*(qmax -
   j/2) + q}.
}

\item{v}{An object of class \code{funlist} whose entries are (spline)
   functions \code{v_q(t)} specifying the (optimal) expected value
   of \code{q} assets at time \code{t} corresponding to the (optimal)
   pricing policy \code{x}.
}
\item{vdot}{An object of class \code{funlist} whose entries are
the derivatives (with respect to \code{t}) of the functions \code{v_q(t)}
described above.  The values of these derivatives are determined
sequentially in the process of solving the system of differential
equations for the optimal pricing policy.
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
\code{\link{vsolve}()}, \code{\link{plot.funlist}()}
}
\examples{
S <- expression(exp(-kappa*x/(1+gamma*exp(-beta*t))))
attr(S,"parlist") <- c(kappa=10/1.5,gamma=9,beta=1)

# Optimal pricing policy assuming customers arrive singly:
lambda1 <- function(tt){
	if(tt<0 | tt> 1) stop("Time out of range.\n")
	84*(1-tt)
}
X1 <- xsolve.cont(type="sip",tmax=1,nstep=1000,qmax=5,lambda=lambda1,
                 S=S,alpha=0.5, gprob=1,verb=TRUE)
# Optimal pricing policy assuming customers arrive in groups of
# size up to 5, with group size probabilities 1/3, 4/15, 1/5, 2/15,
# and 1/15 respectively.
lambda2 <- function(tt){
	if(tt<0 | tt> 1) stop("Time out of range.\n")
	84*(1-tt)
}
X2 <- xsolve.cont(type="sip",tmax=1,nstep=1000,qmax=5,lambda=lambda2,
                 S=S,alpha=0.5, gprob=(5:1)/15,verb=TRUE)

# Note that the intensity functions lambda1() and lambda2() are
# such that the expected total number of customers is 42 in each case.
}
\keyword{ math }
