library(mHG)

system.time(( k = apply(ATAC.clean[c(1:51405),],1, function(x) mHG.test(as.numeric(x> 1))$p.value  ) ) )

#l = mHG.test(c(rep(1,1),rep(0,200)))

#as.numeric(ATAC.clean[3,] > 0)

plot(y = ATAC.clean[which.min(k),],x = c(1:576),col = factor(Cell.type))

plot(y = as.numeric(ATAC.clean[which.min(k),] > 0),x = c(1:576),col = factor(Cell.type))

plot(y = ATAC.clean[which.min(k),],x = c(1:576),col = factor(Cell.type))

dist(t(ATAC.clean))


for(i in which(k< (0.05/51405)  )){
  plot(y = ATAC.clean[i,],x = c(1:576),col = factor(Cell.type))
  
}

plot(y = colSums(ATAC.clean),x = c(1:576),col = factor(Cell.type))

