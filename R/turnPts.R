turnPts <- function(a,b,v,Kpa,xlo,xhi) {
#
# Construct the polynomial.
    q <- length(v)
    rmax  <- suppressWarnings(max(which(Kpa > 0)))
    if(rmax < 1) return(as.polynomial(0))
    Kpa <- Kpa[1:rmax]
    vq <- v[q]
    ply <- 0
    for(r in 1:rmax) {
        vqmr <- if(r < q) v[q-r] else 0
        p1   <- polynomial(c(vqmr-vq,r))
        p2   <- try(polynomial(c(a,b)))
        ply  <- ply + Kpa[r]*p1*p2^r
    }
#
# Take its derivative.
    dply <- deriv(ply)
#
# Find the zeroes.
    zzz  <- try(polyroot(dply))
    rrr  <- sapply(zzz,function(z){isTRUE(all.equal(Im(z),0))})
    if(!any(rrr)) return(numeric(0))
    zzz  <- Re(zzz[rrr])
#
# Return those zeroes that lie in the interval [xlo,xhi]
    ok   <- xlo <= zzz & zzz <= xhi
    zzz[ok]
}
