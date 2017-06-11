library(mHG)

k = apply(ATAC.clean[c(1:500),],1, function(x) mHG.test(as.numeric(x> 0))$p.value  )

l = mHG.test(c(rep(1,1),rep(0,200)))

as.numeric(ATAC.clean[3,] > 0)

plot(y = ATAC.clean[which.min(k),],x = c(1:576),col = factor(Cell.type))
