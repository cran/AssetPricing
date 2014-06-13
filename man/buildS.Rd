\name{buildS}
\alias{buildS}
\title{
   Build a piecewise linear price sensitivity function
}
\description{
   Builds a price sensitivity function which is piecewise linear
   in price, in an automated manner, with built-in checks for
   possible infelicities.
}
\usage{
   buildS(alpha, beta, kn, tmax)
}
\arguments{
  \item{alpha}{
  A list of functions of \code{t} giving the constant terms of
  the linear functions comprising the price sensitivity function.
}
  \item{beta}{
  A list of functions of \code{t} giving the slopes of the linear
  functions comprising the price sensitivity function.
}
  \item{kn}{
  The knots (with respect to price) of the piecewise linear
  price sensitivity function.  The zero knot (which is always
  the first knot) is \emph{not} included in \code{kn}.
}
  \item{tmax}{
  The maximum time value to which the price sensitivity function
  is to be applied.  Needed for internal consistency checks.
}
}
\details{
  The price sensitivity function is assumed to be of the form
  \deqn{S(x,t) = \alpha_k(t) + \beta_k(t)x}{S(x,t) = alpha_k(t) +
  beta_k(t)*x} for \eqn{x_{k-1} \leq x \leq x_k}{x_{k-1} <= x <=
  x_k} where \eqn{x_1, x_2, \ldots, x_K}{x_1, x_2, ..., x_K} are the
  (non-zero) knots of the function.  It is assumed that \eqn{x_0
  = 0}.  The variable \code{x} represents price and the variable
  \code{t} represents residual time.

  The function is defined over the rectangle
  \eqn{[0,x_K] \times [0,t_{\rm max}]}{[0,x_K] x [0,tmax]}.

  Checks are done to make sure that
  \describe{
  \item{\eqn{\bullet}{*} \eqn{S(x,t)} is continuous}{}
  \item{\eqn{\bullet}{*} \eqn{S(0,t) = 1} for all \eqn{t}}{}
  \item{\eqn{\bullet}{*} \eqn{S(x,t)} is non-increasing in \eqn{x} for all t}{}
  \item{\eqn{\bullet}{*} \eqn{S(x,t) \geq 0}{S(x,t) >= 0} for all \eqn{x} and \eqn{t}}{}
  }
}
\value{
   A function of two variables \code{x} and \code{t}, which is
   a price sensitivity function.  The argument \code{x} represents
   price and the argument \code{t} represents (residual) time.
   The value of the function is interpreted as the probability
   that a customer \dQuote{arriving} at time \code{t} will purchase
   an item offered at price \code{x}.
}
\references{
P. K. Banerjee, and T. R. Turner (2012).
A flexible model for the pricing of perishable assets.
\emph{Omega} \bold{40}:5, 533--540,
doi: 10.1016/j.omega.2011.10.001.

Rolf Turner, Pradeep Banerjee and Rayomand Shahlori (2014).  Optimal
Asset Pricing.  \emph{Journal of Statistical Software} \bold{58}:11, 1--25.
URL http://www.jstatsoft.org/v58/i11/.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}
  \url{http://www.stat.auckland.ac.nz/~rolf}
}

\seealso{
\code{\link{xsolve}()}
}
\examples{
   lambda <- function(t) {
      tn <- 1:4
      A <- matrix(c(0,12,12,12,
              0,-16,16,64,
              20,30,30,0),nrow=4)
      B <- matrix(c(12,0,0,0,
               0,16,0,-16,
               0,-10,-10,0),nrow=4)
      s <- cut(t,breaks=c(0,tn),include.lowest=TRUE,labels=tn)
      s <- as.numeric(levels(s)[s])
      M <- matrix(A[s,] + B[s,]*t,ncol=ncol(A))
      M[!is.finite(M)] <- 0
      M
   }

   alpha <- vector("list",4)
   beta  <- vector("list",4)
   alpha[[1]] <- with(list(lambda=lambda),
	function(t) {
	A <- c(1,1,1)
        lll <- lambda(t)
        dnm <- apply(lll,1,sum)
        dnm[dnm==0] <- 1
        lll\%*\%A/dnm
   })
   beta[[1]] <- with(list(lambda=lambda),
	function(t) {
	B <- c(0,0,0)
        lll <- lambda(t)
        dnm <- apply(lll,1,sum)
        dnm[dnm==0] <- 1
        lll\%*\%B/dnm
   })
   alpha[[2]] <- with(list(lambda=lambda),
	function(t) {
	A <- c(1.495,1,1)
        lll <- lambda(t)
        dnm <- apply(lll,1,sum)
        dnm[dnm==0] <- 1
        lll\%*\%A/dnm
   })
   beta[[2]] <- with(list(lambda=lambda),
	function(t) {
	B <- c(-0.2475,0,0)
        lll <- lambda(t)
        dnm <- apply(lll,1,sum)
        dnm[dnm==0] <- 1
        lll\%*\%B/dnm
   })
   alpha[[3]] <- with(list(lambda=lambda),
	function(t) {
	A <- c(0.01,2.485,1)
        lll <- lambda(t)
        dnm <- apply(lll,1,sum)
        dnm[dnm==0] <- 1
        lll\%*\%A/dnm
   })
   beta[[3]] <- with(list(lambda=lambda),
	function(t) {
	B <- c(0,-0.2475,0)
        lll <- lambda(t)
        dnm <- apply(lll,1,sum)
        dnm[dnm==0] <- 1
        lll\%*\%B/dnm
   })
   alpha[[4]] <- with(list(lambda=lambda),
	function(t) {
	A <- c(0.01,0.01,3.475)
        lll <- lambda(t)
        dnm <- apply(lll,1,sum)
        dnm[dnm==0] <- 1
        lll\%*\%A/dnm
   })
   beta[[4]] <- with(list(lambda=lambda),
	function(t) {
	B <- c(0,0,-0.2475)
        lll <- lambda(t)
        dnm <- apply(lll,1,sum)
        dnm[dnm==0] <- 1
        lll\%*\%B/dnm
   })
   kn <- c(2,6,10,14)
   S  <- buildS(alpha,beta,kn,4)
   x  <- seq(0,14,length=41)
   t  <- seq(0,4,length=41)
   z  <- S(x,t)
\dontrun{
   persp(x,t,z,theta=150,phi=40,d=4,xlab="price",ylab="time",
         zlab="probability",ticktype="detailed")
}
}
\keyword{ utilities }
