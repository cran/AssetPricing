turnPts <- function(a,b,v,Kpa,xlo,xhi) {
#
# Construct the polynomial.
    q <- length(v)
    jmax  <- suppressWarnings(max(which(Kpa > 0)))
    if(jmax < 1) return(as.polynomial(0))
    Kpa <- Kpa[1:jmax]
    vq <- v[q]
    ply <- 0
    for(j in 1:jmax) {
        vqmj <- if(j < q) v[q-j] else 0
        p1   <- polynomial(c(vqmj-vq,j))
        p2   <- polynomial(c(a,b))
        ply  <- ply + Kpa[j]*p1*p2^j
    }
#
# Take its derivative.
    dply <- deriv(as.polylist(ply))[[1]]
#
# Find the zeroes.
    zzz  <- polyroot(dply)
    rrr  <- sapply(zzz,function(z){isTRUE(all.equal(Im(z),0))})
    if(!any(rrr)) return(numeric(0))
    zzz  <- Re(zzz[rrr])
#
# Return those zeroes that lie in the interval [xlo,xhi]
    ok   <- xlo <= zzz & zzz <= xhi
    zzz[ok]
}
