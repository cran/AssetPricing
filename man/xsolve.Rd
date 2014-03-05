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
   xsolve(S, lambda, gprob=1, tmax=NULL, qmax, prices=NULL, nout=300,
          type="sip", alpha=NULL, salval=0, epsilon=NULL,
          method="lsoda",verbInt=0)
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
  \item{nout}{
  The number of points at which values of the solution are to be
  provided.  These are taken to be equispaced on
  [0,\code{tmax}].
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
  \item{epsilon}{
  A numeric scalar used in determining the optimal price
  in settings in which this involves maximizing over a discrete set.
  See \bold{Details}.  It defaults to \code{.Machine$double.eps^0.25}
  in the case of discrete prices and to \code{.Machine$double.eps^0.5}
  when the price sensitivity function is piecewise linear.
  It is ignored if the price sensitivity function is smooth.
}
  \item{method}{
  Character string specified the solution method to be used
  by the differential equation solver \code{ode()}.  There is
  a fairly large number of possible methods.  See
  \code{\link[deSolve]{ode}()} for details.
}
  \item{verbInt}{
  A scalar value which controls \dQuote{verbosity}.
  If \code{verbInt > 0} then a \dQuote{progress report} is
  printed every \code{verbInt} \emph{seconds} (roughly).  That is
  if the current value of \code{Sys.time()} is greater than or
  equal to the value stored at the time of the last report, plus
  \code{verbInt} seconds, then a new \dQuote{report} is printed out.
  If \code{verbInt} is less than or equal to 0 then the solution
  process runs \dQuote{silently}.  See section \bold{Progress
  Reports} for a bit more detail.
}
}
\details{
  If prices are modelled as being continuous, and if the price
  sensitivity function is differentiable, a coupled system
  of differential equations for the optimal prices is solved.
  If the prices are modelled as being discrete or if the price
  sensitivity function is piecewise linear in price, then a
  coupled system of differential equations for the expected
  value of the stock is solved, with the optimal price being
  determined at each step by maximizing over an appropriate finite
  discrete set.  These differential equations are solved by the
  \code{\link[deSolve]{ode}()} function from the \code{deSolve}
  package.

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

  In the case where argument \code{S} is a piecewise linear
  price sensitivity function, the argument \code{tmax} is, if not
  specified, taken to be the value of the corresponding attribute
  of \code{S}.  In this setting, if \code{tmax} \emph{is} specified
  it must be less than or equal to the corresponding attribute
  of \code{S}.

  For discrete prices and for piecewise linear price sensitivity
  functions, determining the optimal price involves maximizing
  expected values over finite discrete sets.  It can happen that the
  location of the maximum can make a sudden \dQuote{jump} from one
  time step to the next, causing anomalous looking discontinuities
  in the optimal price functions.  To avoid this, we check on the
  change in the expected value at each of the possible new prices
  as compared with that at the \dQuote{previous} price.

  If the maximal \dQuote{improvement} in expected value is less
  than or equal to \code{epsilon} then the \dQuote{new} price is
  set equal to the previous value.  If the maximal improvement
  is greater than \code{epsilon} then those values of price,
  where the expected value is greater than the maximum value minus
  \code{epsilon}, are considered and the one which is closest to
  the previous price is chosen.

  If \code{epsilon} is set equal to a value less than or equal to
  \code{0} then the smoothing strategy described above is dispensed
  with. In this case the maximum is taken to be the first of the
  (possibly) multiple maxima of the expected value.
}
\section{Note}{
A substantial change was made to this package as of the change
of version number from 0.0-11 to 0.1-0.  Previously the differential
equations which arise were solved via a \dQuote{locally produced}
roll-your-own Runge-Kutta procedure.  Now they are solved (in a
more sophisticated manner) using the package \code{deSolve}.  This
increases the solution speed by a factor of about 7.  There will
be (minor, it is to be hoped) numerical differences in solutions
produced from the same input.
}

\section{Progress Reports}{
  The \dQuote{progress reports} produced when \code{verbInt > 0}
  consist of rough estimates of the percentage of \code{[0,tmax]}
  (the interval over which the differential equation is being solved)
  remaining to be covered.  A rough estimate of the total elapsed
  time since the solution process started is also printed out.

  Having \dQuote{progress reports} printed out appears to have no
  (or at worst negligible) impact on computation time.`
}

\value{
A list with components:
\item{x}{The optimal pricing policy, chosen to maximize the expected
   value of the remaining assets at any given time; an object of
   class \code{flap} (\dQuote{function list for asset pricing}).  (In
   the case of discrete prices it also inherits from \code{pwc.flap}
   (\code{pwc} stands for \dQuote{piecewise constant}), and if
   \code{type=="dip"} it also inherits from \code{di.flap}.)
   If \code{type=="sip"} it has the form of a list of functions
   \code{x_q(t)}, with \code{q} running from 1 to \code{qmax}
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
  \url{http://www.stat.auckland.ac.nz/~rolf}
}
\seealso{
\code{\link{vsolve}()}, \code{\link{plot.AssetPricing}()},
\code{\link{buildS}()}
}
\examples{
#
# In these examples "qmax" has been set equal to 5 which is
# an unrealistically low value for the total number of assets.
# This is done so as to reduce the time for package checking on CRAN.
#
# Smooth price sensitivity function.
S <- expression(exp(-kappa*x/(1+gamma*exp(-beta*t))))
attr(S,"parvec") <- c(kappa=10/1.5,gamma=9,beta=1)

# Optimal pricing policy assuming customers arrive singly:
lambda1 <- function(tt){
	84*(1-tt)
}
X1 <- xsolve(S=S,lambda=lambda1,gprob=1,tmax=1,qmax=5,
             type="sip",verbInt=5)
# Optimal pricing policy assuming customers arrive in groups of
# size up to 5, with group size probabilities 1/3, 4/15, 1/5, 2/15,
# and 1/15 respectively.
lambda2 <- function(tt){
	36*(1-tt)
}
X2 <- xsolve(S=S,lambda=lambda2,gprob=(5:1)/15,tmax=1,qmax=5,
             type="sip", alpha=0.5,verbInt=5)

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
             type="sip",verbInt=5)

# Piecewise linear price sensitivity function.
#
# Take S as in the example for buildS.
# This takes a loonnngggg time; the procedure is slow
# in the piecewise linear setting.
\dontrun{
l0 <- get("lambda",envir=environment(get("alpha",envir=environment(S))[[1]]))
lambda4 <- function(t){apply(l0(t),1,sum)}
X4 <- xsolve(S=S,lambda=lambda4,gprob=(5:1)/15,qmax=30,type="sip",
                   alpha=0.5,verbInt=20)
}
}
\keyword{ math }
