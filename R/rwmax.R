rwmax <- function(x,eps=sqrt(.Machine$double.eps)) {
#
# Local version of which.max() which is robust to
# numerical noise.
#
i <- which.max(x)[1]
which(x >= x[i] - eps)
}
