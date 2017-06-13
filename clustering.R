
ATAC.clean.filter = ATAC.clean[apply(ATAC.clean,1,function(x) sum(x>2))>=10,]


#plot(density(colSums(ATAC.clean)))

pca.res = PCA(t(   log2(ATAC.clean.norm.filter+0.5)   ),graph = F)
PCs <- data.frame(PC1=as.numeric(pca.res$ind$coord[,1]),
                  PC2=as.numeric(pca.res$ind$coord[,2]),
                  sample.names = Cell.type)

P<-ggplot(PCs, aes(PC1,PC2,colour = sample.names)) 
P +geom_point()

ATAC.clean.norm = sweep(ATAC.clean, 2, colSums(ATAC.clean), `/`)*1000000

ATAC.clean.norm.filter = ATAC.clean.norm[apply(ATAC.clean.norm,1,function(x) sum(x>20))>=10,]

library(Rtsne)



ATAC.clean.norm.filter.unique = unique(  t(ATAC.clean.norm.filter)  )

b = (  log2( ATAC.clean.norm.filter.unique +0.5))



tsne.obj = Rtsne(    b , perplexity =10)

plot.data = data.frame(Y1 = tsne.obj$Y[,1],Y2 = tsne.obj$Y[,2],nam = rownames(b))
plot(plot.data$Y1,plot.data$Y2,col = factor(plot.data$nam))

plot(plot.data$Y1,plot.data$Y2)
