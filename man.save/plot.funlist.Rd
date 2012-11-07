\name{plot.funlist}
\alias{plot.funlist}
\title{
   Plot a list of functions.
}
\description{
   Plot a list of functions --- in particular expected value
   functions or optimal price functions (which are functions
   of residual time), the index of the list being the number
   of items available for sale.
}
\usage{
\method{plot}{funlist}(x,xlim=NULL,ylim=NULL,lty=NULL,cols=NULL,
             xlab=NULL,ylab=NULL,main=NULL,indices=1:length(x),
             add=FALSE,aao=TRUE,gloss=NULL,glind=NULL,
             col.gloss=1,mfrow=NULL,...)
}
\arguments{
  \item{x}{
  An object of class \code{funlist} such as the components of
  the objects produced by \code{\link{vsolve}()},% \code{\link{xsolve.disc}()}
  or \code{\link{xsolve.cont}()}.%, or \code{\link{xsolve.pwl}()}.
}
  \item{xlim}{
  The \code{x} limits of the plot.  Defaults to the \code{xlim}
  attribute of the object \code{x}.  If this attribute does not
  exist and \code{xlim} is not supplied then an error is given.
}
  \item{ylim}{
  The \code{y} limits of the plot.  Defaults to the \code{ylim}
  attribute of the object \code{x}.  If this attribute does not
  exist and \code{ylim} is not supplied then an error is given.
}
\item{lty}{
  Vector of line types for the plotted lines.  Replicated to yield a
  vector of the same length as \code{indices}.  Defaults to \code{1},
  i.e. solid lines.
}
  \item{cols}{
  Vector of colours for the plotted lines.  Replicated to yield a
  vector of the same length as \code{indices}.  Defaults to \code{1},
  i.e. black.
}
  \item{xlab}{
  A text string giving a label for the \code{x} axis.  Defaults to
  the null string.  Ignored if \code{add} is \code{TRUE}.
}
  \item{ylab}{
  A text string giving a label for the \code{y} axis.  Defaults to
  the null string.  Ignored if \code{add} is \code{TRUE}.
}
  \item{main}{
  A text string giving an overall title for the plot.  Defaults to
  the null string.  Ignored if \code{add} is \code{TRUE}.
}
  \item{indices}{
  A vector of indices indicating which functions in the list
  are to be plotted.  Should be a subset of \code{1:n} where
  \code{n} is the length of the list of functions.
}
  \item{add}{
  Logical scalar; should the plot be added to an existing plot?
}
  \item{aao}{
  Logical scalar; should the functions be plotted ``all at once'',
  i.e. all superimposed in the same plot.  If \code{aao} is \code{FALSE}
  each function is plotted in a separate panel.  The panels are
  arranged in arrays whose dimension is determined by \code{mfrow}.
}
  \item{gloss}{
  Either (1) A logical scalar (should a ``marginal gloss'' be
  added to the plot?) in which case a gloss (consisting of the
  index of the function) is added in the right hand margin for
  every function whose index is in \code{indices[glind]}.  Or (2)
  A vector of character strings of which the marginal gloss is to
  consist.  This vector must be of the same length as \code{glind}.
  If \code{gloss} is \code{NULL} then no marginal gloss is added.
  The \code{gloss} argument has an effect only if \code{aao} is
  \code{TRUE}, and is otherwise ignored.  Note that if \code{add}
  is \code{TRUE} then the gloss may not actually appear in the
  plot, since it is placed in the right hand margin and may
  consequently be outside of the plotting region.  Thus if you
  wish to use a gloss when adding to an existing plot you will
  probably need to take steps to ensure that there is room in
  the right hand margin for the plot to appear, or possibly
  set \code{par(xpd=NA)}.
}
  \item{glind}{
  A vector of integers which are a subset of
  \code{1:length(indices)}.  These integers indicate the indices
  for which marginal gloss is to be added.  I.e. marginal gloss
  is added for the graphs of functions whose indices are in
  \code{indices[glind]}.
}
  \item{col.gloss}{
  Scalar specifying the colour in which the marginal gloss
  is to be added, e.g. \dQuote{red} (or equivalently 2).
  Defaults to 1, i.e. \dQuote{black}.
}
  \item{mfrow}{
  The dimensions of the array(s) of panels in which the functions
  are plotted if \code{aao} is \code{FALSE}.
}
  \item{\dots}{
  Extra arguments to be passed to \code{plot} (effectively
  to \code{plot.function} or to \code{plot.stepfunction}).
}
}
\value{
None.  This function exists only for its side effect, i.e.
the production of a plot or plots.
}
\references{
Banerjee, P. K., and Turner, Rolf (2009).  A flexible
model for the pricing of perishable assets.  Submitted
for publication.
}
\author{Rolf Turner
   \email{r.turner@auckland.ac.nz}
   \url{http://www.math.unb.ca/~rolf}
}
\seealso{
%\code{\link{xsolve.disc}()},
\code{\link{xsolve.cont}()},
%\code{\link{xsolve.pwl}()}
\code{\link{vsolve}()},
}
\examples{
\dontrun{
S       <- expression(exp(-kappa*x/(1+gamma*exp(-beta*t))))
attr(S,"parlist") <- c(kappa=10/1.5,gamma=9,beta=1)
LAMBDA  <- function(tt){
               if(tt<0 | tt> 1) stop("Time out of range.\n")
               36*(1-tt)
}
OUT     <- xsolve.cont(type="sip",tmax=1,nstep=1000,qmax=30,
                       lambda=LAMBDA,S=S,alpha=0.5,
                       gprob=(5:1)/15,verb=TRUE)
plot(OUT$v,indices=5*(1:6),gloss=TRUE,
     xlab="residual time",ylab="expected revenue")
plot(OUT$v,indices=5*(1:6),gloss=c("five","ten","fifteen","twenty",
     "twenty-five","thirty"),xlab="residual time",ylab="expected revenue")
}
}
\keyword{ hplot }
