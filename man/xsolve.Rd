\name{xsolve}
\alias{xsolve}
\title{
  Optimal pricing policy
}
\description{
  Determines (by solving a coupled system of differential equations)
  the optimal prices as functions of (residual) time for a
  number perishable assets. Prices may be discrete or continuous.
  For continuous prices, the price sensitivity function may either
  be a smooth (twice differentiable) function or a function which
  is piecewise linear in price.
}
\usage{
   xsolve(S, lambda, gprob=1, tmax=NULL, qmax, prices=NULL, nstep, type="sip",
                   alpha=NULL, salval=0, maxFac=1000, nverb=0)
}
\arguments{
  \item{S}{
  An expression, or list of expressions, or a function or list
  of functions, specifying the price sensitivity functions
  \code{S_j(x,t)}.  See \bold{Details}.
}
  \item{lambda}{
  A function (of residual time \code{t} --- see \code{tmax}) or
  a positive constant, specifying the intensity of the (generally
  inhomogeneous) Poisson process of arrival times of groups of
  potential customers.
}
  \item{gprob}{
  A function (to calculate probabilities) or a numeric vector of 
  probabilities determining the distribution of the size of
  an arriving group of customers.
}
  \item{tmax}{
  The maximum residual time; think of this as being the initial time
  at which the assets go on sale (with time \emph{decreasing} to zero,
  at which point the value of each asset drops to the \dQuote{salvage
  value} (\code{salval}), usually \code{0}).  The system of differential
  equations is solved over the time interval [0,\code{tmax}].
  See \bold{Details}.
}
  \item{qmax}{
  The maximum number of assets available for sale, i.e. the number
  of assets available at the starting (residual) time \code{tmax}.
}
  \item{prices}{
  A numeric vector (with positive values) listing the possible
  prices at which items may be offered for sale in the discrete
  pricing scenario.
}
  \item{nstep}{
  The number of (equal) sub-intervals into which the interval
  [0,\code{tmax}] is divided in the process of solving the system
  of differential equations numerically.
}
  \item{type}{
  Scalar character string taking one of the two values \code{"sip"}
  (singly indexed prices) and \code{"dip"} (doubly indexed prices).
  In the \code{"dip"} case the price of the asset which is quoted
  to the arriving group is allowed to depend upon the group size
  (as well as upon (residual) time as well as on the number of
  assets remaining for sale.  In the \code{"sip"} case the quoted
  price does not depend upon group size.
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
  \item{maxFac}{
  A (positive) numeric scalar used in determining the optimal price
  in settings in which this involves maximizing over a discrete set.
  See \bold{Details}.  Ignored if the price sensitivity function is
  smooth.
}
  \item{nverb}{
  Non-negative integer; if \code{nverb > 0} then a \dQuote{progress
  report} (actually consisting solely of the step number) is
  printed out every \code{nverb} steps of the Runge-Kutta procedure.
  (If \code{nverb==0} the procedure runs silently.)
}
}
\details{
  If prices are modelled as being continuous, and if the price
  sensitivity function is differentiable, a coupled system of
  differential equations for the optimal prices is solved.  If the
  prices are modelled as being discrete or if the price sensitivity
  function is piecewise linear in price, then a coupled system of
  differential equations for the expected value of the stock is
  solved, with the optimal price being determined at each step by
  maximizing over an appropriate finite discrete set. The systems
  of differential equations are treated in vectorised form and are
  solved via the method of Runge-Kutta.

  The components of the argument \code{S} should be provided as
  expressions when the price sensitivity functions are assumed to
  be smooth, and these should be amenable to differentiation with
  respect to \code{x} and \code{t} via the function \code{deriv3()}.

  Note that in general the expression or expressions will depend upon
  a number of \emph{parameters} as well as upon \code{x} and \code{t}.
  The values of these parameters are specified via an attribute or
  attributes.  If \code{S} is a (single) expression it has (must have)
  an attribute called \code{parvec} which is a \emph{named} vector of
  parameter values.  If \code{S} is a list of expressions each entry
  of the list has (must have) such an attribute.

  In the \dQuote{piecewise linear} context \code{S} can be
  specified \emph{only} as a single function.  It is then assumed
  that the price sensitivity function for a group of size \code{j}
  is given by \code{S_j(x,t) = S(x,t)^j}.  Such piecewise linear
  price sensitivity functions should be built using the function
  \code{\link{buildS}()}.

  In the case of discrete prices the argument \code{S} must be a
  function or list of functions specifying the price sensitivity
  functions \code{S_j(x,t)}.  These functions need only be defined
  for the prices listed in the \code{prices} argument.

  If \code{S} is a single expression or function, then
  \code{S_j(x,t)} is taken to be this expression or function raised
  to the power \code{j}.  If \code{S} is a list, then \code{S_j(x,t)}
  is taken to be its \code{j}-th entry.

  For discrete prices and for piecewise linear price sensitivity
  functions, determining the optimal price involves maximizing
  expected values over finite discrete sets.  It can happen that
  there are a number of values \dQuote{close} to the maximum, and
  as a consequence the location of the maximum can make a sudden
  \dQuote{jump} from one time step to the next, causing anomalous
  looking discontinuities in the optimal price functions.  To avoid
  this, the \dQuote{optimal} price is chosen to be that one from
  among a suitable subset of the candidate prices which is closest
  to the price at the previous time step.

  The subset consists of those prices yielding expected values
  which are within \eqn{\epsilon}{epsilon} of the maximum.  Here
  \eqn{\epsilon}{epsilon} is equal to the range (i.e. \dQuote{max
  minus min}) of the possible expected values, divided by
  \code{maxFac}.  Setting \code{maxFac=Inf} has the effect of
  simply picking the maximal value (or the first maximal value if
  there are a number of such values that are numerically equal).
  Using different values of \code{maxFac} may have an impact on
  the optimal price functions but will generally have little or no
  discernible impact on the corresponding expected values of stocks.

  In the case where argument \code{S} is a piecewise linear
  price sensitivity function, the argument \code{tmax} is, if not
  specified, taken to be the value of the corresponding attribute
  of \code{S}.  In this setting, if \code{tmax} \emph{is} specified
  it must be less than or equal to the corresponding attribute
  of \code{S}.
}

\value{
A list with components:
\item{x}{The optimal pricing policy, chosen to maximize the expected
   value of the remaining assets at any given time; an object of
   class \code{flap} (\dQuote{function list for asset pricing}).  (In
   the case of discrete prices it also inherits from \code{pwc.flap}
   (\code{pwc} stands for \dQuote{piecewise constant}), and if
   \code{type=="dip"} it also inherits from \code{di.flap}.)
   It has the form of a list of functions \code{x_q(t)}, with
   \code{q} running from 1 to \code{qmax} if \code{type=="sip"}.
   If \code{type=="dip"} if has the form of a list of functions
   \code{x_qj(t)} with \code{q} running from 1 to \code{qmax} and
   \code{j} running from 1 to \code{min(q,jmax)} where \code{jmax}
   is the maximum group size.

   In the case of continuous prices these functions will be
   continuous functions created by \code{\link{splinefun}()}.  In the
   case of discrete prices these functions will be piecewise constant
   (of class \code{stepfun}) created by \code{\link{stepfun}()}.

   Note that \code{x} has an \emph{attribute} \code{qmax} specifying
   the maximum number of assets available for sale, i.e. the number
   of assets available at the starting (residual) time \code{tmax}.
   Of course if \code{type=="sip"} then this attribute is simply
   the length of \code{x}.  Note that if \code{type == "dip"} then
   the entry \code{x[[i]]} is equal to the function \code{x_qj(t)}
   where \code{i = (j-1)*(qmax - j/2) + q}.
}

\item{v}{An object of class \code{flap} whose entries are (spline)
   functions \code{v_q(t)} specifying the (optimal) expected value
   of \code{q} assets at time \code{t} corresponding to the (optimal)
   pricing policy \code{x}.
}
\item{vdot}{An object of class \code{flap} whose entries are
   the derivatives (with respect to \code{t}) of the functions
   \code{v_q(t)} described above.  The values of these derivatives
   are determined sequentially in the process of solving the system
   of differential equations for the optimal pricing policy.
}
}
\references{
Banerjee, P. K. and Turner, T. R.
A flexible model for the pricing of perishable assets.
Omega, vol. 40, number 5, 2012, pages 533 -- 540,
doi: 10.1016/j.omega.2011.10.001.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}
  \url{http://www.math.unb.ca/~rolf}
}
\seealso{
\code{\link{vsolve}()}, \code{\link{plot.flap}()}, \code{\link{buildS}()}
}
\examples{
#
# In these examples "nstep" has been set equal to 100
# to reduce the time for package checking on CRAN.
# A more "realistic" value is nstep=1000.  Likewise "qmax"
# has been set equal to 5 which is an unrealistically low
# value for the total number of assets.
#
# Smooth price sensitivity function.
S <- expression(exp(-kappa*x/(1+gamma*exp(-beta*t))))
attr(S,"parvec") <- c(kappa=10/1.5,gamma=9,beta=1)

# Optimal pricing policy assuming customers arrive singly:
lambda1 <- function(tt){
	84*(1-tt)
}
X1 <- xsolve(S=S,lambda=lambda1,gprob=1,tmax=1,qmax=5,nstep=100,
             type="sip",nverb=50)
# Optimal pricing policy assuming customers arrive in groups of
# size up to 5, with group size probabilities 1/3, 4/15, 1/5, 2/15,
# and 1/15 respectively.
lambda2 <- function(tt){
	36*(1-tt)
}
X2 <- xsolve(S=S,lambda=lambda2,gprob=(5:1)/15,tmax=1,qmax=5,nstep=100,
             type="sip", alpha=0.5,nverb=50)

# Note that the intensity functions lambda1() and lambda2() are
# such that the expected total number of customers is 42 in each case.

# Discrete prices:
lambda3 <- function(t){42}
S      <- function(x,t){
                e <- numeric(2)
                e[x==1] <- exp(-2*t)
                e[x==0.6] <- 1.0
                e
          }
X3 <- xsolve(S=S,lambda=lambda3,gprob=1,tmax=1,qmax=5,prices=c(1,0.6),
             nstep=100,type="sip",nverb=50)

# Piecewise linear price sensitivity function.
#
# Take S as in the example for buildS.
# This takes a loonnngggg time; the procedure is slow
# in the piecewise linear setting.
\dontrun{
l0 <- get("lambda",envir=environment(get("alpha",envir=environment(S))[[1]]))
lambda4 <- function(t){apply(l0(t),1,sum)}
X4 <- xsolve(S=S,lambda=lambda4,gprob=(5:1)/15,qmax=30,nstep=1000,type="sip",
                   alpha=0.5,nverb=50)
}
}
\keyword{ math }
