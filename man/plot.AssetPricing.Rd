\name{plot.AssetPricing}
\alias{plot.AssetPricing}
\title{
   Plot a list of asset pricing functions.
}
\description{
   Plot a list of functions --- in particular optimal price
   functions or expected value functions or derivatives of the
   expected value functions.  Such a list is assumed to occur as
   a component of an object produced by \code{\link{xsolve}()} of
   \code{\link{vsolve}()}.  The functions in the list are functions
   of residual time.  The indices of the list correspond to the
   number of items available for sale and possibly (for optimal
   price functions) the size of the arriving group of customers.
}
\usage{
\method{plot}{AssetPricing}(x,witch=c("price","expVal","vdot"),
             xlim=NULL,ylim=NULL,lty=NULL,cols=NULL,xlab=NULL,
             ylab=NULL,main=NULL,main.panel= NULL,groups=NULL,
             add=FALSE,gloss=FALSE,glind=NULL,extend=0.3,col.gloss=1,
             cex.gloss=0.8,mfrow=NULL,...)
}
\arguments{
  \item{x}{
  An object of class \code{AssetPricing}, i.e. an object
  produced by \code{\link{vsolve}()}, or \code{\link{xsolve}()}.
}
  \item{witch}{
  A text string indicating which of the three possible components
  of \code{x} should be plotted.  May be abbreviated, e.g. to
  \code{p}, \code{e} or \code{v}.
}
  \item{xlim}{
  The \code{x} limits of the plot.  Defaults to the \code{tlim}
  attribute of the object \code{x[[witch]]}.  If this attribute does not
  exist and \code{xlim} is not supplied then an error is given.
}
  \item{ylim}{
  The \code{y} limits of the plot.  Defaults to the \code{ylim}
  attribute of the object \code{x[[witch]]}.  If this attribute does not
  exist and \code{ylim} is not supplied then an error is given.
}
\item{lty}{
  A vector of line types.  It will be replicated to have a length
  equal to the number of rows of \code{groups} (see below).
  Defaults to having all entries of the vector equal to \code{1},
  i.e. solid lines.
}
  \item{cols}{
  A vector of colours for the plotted lines.  It will be replicated
  to have a length equal to the number of rows of \code{groups}
  (see below).  Defaults to having all entries of the vector equal
  to \code{1}, i.e. black.
}
  \item{xlab}{
  A text string giving a label for the \code{x} axis (or axes).
  Defaults to the null string.  Ignored if \code{add} is \code{TRUE}.
}
  \item{ylab}{
  A text string giving a label for the \code{y} axis (or axes).
  Defaults to the null string.  Ignored if \code{add} is \code{TRUE}.
}
  \item{main}{
  A text string giving an overall title for the plot or for each
  page of plots if there is more than one.  Defaults to the null
  string and is ignored if \code{add} is \code{TRUE}.
}
  \item{main.panel}{
  A text string which is replicated \dQuote{np} times (where
  \dQuote{np} is the total number of panels) or a vector of text
  strings of length equal to \dQuote{np}.  Note that \dQuote{np} will
  be equal to the number of unique entries of \code{groups$group}.
  (See below.)  The \eqn{i^{th}}{i-th} entry of the vector is
  used as the title of the \eqn{i^{th}}{i-th} panel of the plots
  that are created.  If \code{main.panel} is left \code{NULL}
  the \eqn{i^{th}}{i-th} entry of the vector is set equal to
  \code{paste("group",i)}.  This argument is ignored if there
  is only a single panel.
}
  \item{groups}{
  A data frame with one, two or three columns, named \code{group},
  \code{q} and \code{j}.  The total number of rows should be less
  than or equal to the total number of entries of the function
  list \code{x[[witch]]}.  Only those function traces corresponding to a
  row of \code{groups} are plotted.  The traces corresponding to
  an individual value in the \code{group} column are plotted in the
  same panel of a multi-panel array of plots.  See \bold{Details}.
}
  \item{add}{ Logical scalar; should the plot be added to an
  existing plot?
}
  \item{gloss}{ Either a logical scalar (should a ``marginal
  gloss'' be added to the plot? --- if \code{TRUE} then the gloss
  is constructed internally; see \bold{Details})  or a vector of
  character strings of which the marginal gloss is to consist.
}
  \item{glind}{ A logical vector indicating which entries of
  \code{gloss} should actually be used (plotted).  I.e. marginal
  gloss is added for the graphs of functions whose corresponding
  values in the entries of \code{glind} are \code{TRUE}.  Ignored if
  \code{gloss} is \code{FALSE}.  If \code{gloss} is \code{TRUE} or
  is explicitly provided, then if \code{glind} is not specified it
  defaults to a vector, of the same length as \code{gloss} all of
  whose entries are \code{TRUE}.
}
  \item{extend}{ A scalar, between 0 and 1, indicating how much
  the \code{x}-axis should be extended (to the right) in order to
  accommodate the marginal gloss.
}
  \item{col.gloss}{ Scalar specifying the colour in which the
  marginal gloss is to be added, e.g. \code{"red"} (or equivalently
  \code{2}).  The default, i.e. \code{1}, is black.
}
  \item{cex.gloss}{ Character expansion (\code{cex}) specifier for
  the marginal gloss.
}
  \item{mfrow}{ The dimensions of the array(s) of panels in which
  the functions are plotted.
}
  \item{\dots}{ Extra arguments to be passed to \code{plot}
  (effectively to \code{\link{plot.function}()} or to
  \code{\link{plot.stepfun}()}).
} } \details{
  If the argument \code{groups} is specified then:
  \itemize{
  \item it must \emph{always} have a column \code{q}.  The values
  in this column should be integers between \code{1} and \code{qmax}
  (see below).

  \item if \code{jmax} (see below) is greater than 1 it must also
  have a column \code{j}.  The entries of this column should be
  integers between \code{1} and \code{jmax}.

  \item if \code{jmax} is equal to \code{1} then column \code{j}
  need not be present.  In this case, it is internally set equal
  to a column of \code{1}-s.

  \item if the \code{group} column is present its entries
  should be (consecutive) positive integers running from \code{1}
  to the total number of groups.

  \item if the \code{group} column is not present then this column
  is internally set equal to a column of \code{1}-s i.e. there is
  a single group of traces.
  }

  The value of \code{qmax} is the maximum number of items that
  are available for sale in the time period under consideration.
  It may be obtained as \code{attr(x,"qmax")}.

  The value of \code{jmax} is, when \dQuote{double indexing}
  applies, the maximum size of an arriving group of customers,
  and is otherwise equal to \code{1}.  It may be obtained as
  \code{attr(x,"jmax")}.  Note that \dQuote{double indexing} can
  \emph{only} apply when \code{x[[witch]]} is a list of \emph{price}
  functions, i.e. when \code{witch} is equal to \code{price}.
  Hence \dQuote{double indexing} does not apply when \code{witch}
  is equal to \code{expVal} or to \code{vdot}.  In these cases
  \code{jmax} is equal to \code{1}.

  If \code{groups} is not specified then it defaults to a data
  frame with number of rows equal to the length of \code{x[[witch]]}, The
  \code{group} column has entries all equal to 1, i.e. there is
  a single group of traces.  The \code{q} and \code{j} columns
  contain all possible (valid) combinations of stock size and
  customer group size.

  If \code{gloss} is \code{FALSE} then no marginal gloss is plotted.
  If \code{gloss} is \code{TRUE} then the marginal gloss is created
  from the values of the \code{q} and \code{j} entries in the
  columns of \code{groups} using \code{paste()}.

  Note that if \code{add} is \code{TRUE} then the gloss may not
  actually appear in the plot, since it is placed at the right
  hand edge of the plot and may consequently be outside of the
  plotting region.  Thus if you wish to use a gloss when adding to
  an existing plot you will probably need to take steps to ensure
  that there is room in the right hand margin for the plot to appear,
  or possibly set \code{par(xpd=NA)}.

  If \dQuote{double indexing} applies then \code{x[[i]]} corresponds
  to a stock size of \code{q} and a customer group size of \code{j}
  where \code{i = (j-1)*(qmax - j/2) + q}.

  To get traces plotted in individual panels (one trace per panel)
  set the \code{group} column of \code{groups} to be \code{1:n}
  where \code{n} is the total number of traces being plotted.

  This function (i.e. \code{plot.AssetPricing()} calls upon an
  \dQuote{internal} function \code{plot.flap()} to do the hard
  yakka. (Note that \code{flap} stands for dQuote{function list
  for asset pricing}.)

  The function \code{plot.flap()} makes use of a modified
  version of \code{plot.stepfun()}, rather than the one which
  appears in \code{package:stats}.  The modification causes
  \code{plot.stepfun()} to treat the \code{xlim} argument in a manner
  similar to the way in which it is treated by \code{plot.function}.
  Note that \code{plot.stepfun()} is \emph{not} exported from
  this package.
}
\value{
None.  This function exists only for its side effect, i.e.
the production of a plot or plots.
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
\code{\link{xsolve}()}
\code{\link{vsolve}()},
}
\examples{
\dontrun{
S <- expression(exp(-kappa*x/(1+gamma*exp(-beta*t))))
attr(S,"parvec") <- c(kappa=10/1.5,gamma=9,beta=1)
LAMBDA <- function(tt){
    if(tt<0 | tt> 1) stop("Time out of range.\n")
    36*(1-tt)
}
OUT <- xsolve(S=S,lambda=LAMBDA,gprob=(5:1)/15,tmax=1,qmax=30,
                  nstep=1000,alpha=0.5,type="dip",nverb=20)
GLND <- rep(FALSE,30)
GLND[c(1:5,10,15,20,30)] <- TRUE
plot(OUT$v,xlab="residual time",ylab="expected revenue",
     gloss=TRUE,glind=GLND)
GRPS <- data.frame(group=rep(1:6,each=5),q=1:30)
GLND <- c(TRUE,FALSE,TRUE,FALSE,TRUE,rep(c(rep(FALSE,4),TRUE),5))
plot(OUT$v,groups=GRPS,xlab="residual time",ylab="expected revenue",
     gloss=TRUE,glind=GLND)
GRPS <- data.frame(group=rep(1:5,each=6),j=rep(1:5,each=6))
GRPS$q <- with(GRPS,pmax(j,rep(c(1,6,11,16,21,26),5)))
GLND <- rep(c(TRUE,TRUE,rep(FALSE,3),TRUE),5)
plot(OUT$x,groups=GRPS,mfrow=c(3,2),gloss=TRUE,glind=GLND)
# Pretty messy looking:
GRPS$group <- 1
GLND <- unlist(lapply(1:5,function(k){(1:6)==k}))
plot(OUT$x,groups=GRPS,gloss=TRUE,glind=GLND,cols=GRPS$j)
}
}
\keyword{ hplot }
