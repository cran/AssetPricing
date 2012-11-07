getPossPrices <- function(v,t,alpha,beta,kn,Kpa) {
K <- length(kn)
kn <- c(0,kn)
rslt <- vector("list",K)
for(k in 1:K) {
	a <- alpha[[k]](t)
	b <- beta[[k]](t)
	xlo <- kn[k]
	xhi <- kn[k+1]
	rslt[[k]] <- turnPts(a,b,v,Kpa,xlo,xhi)
}
sort(unique(c(unlist(rslt),kn)))
}
