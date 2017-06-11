# load 
library(data.table)
library(FactoMineR)
library(ggplot2)
library(ggrepel)
library(Matrix)
library('irlba')
library(bigpca)


ATAC = fread("GSE74310_scATACseq_All_Counts.txt",header = T)
ATAC = as.data.frame(ATAC)

# prepare cell type tag
Cell.type =   sapply(strsplit(colnames(ATAC), split='-',fixed=TRUE),function(x) paste(x[2],x[3],sep="-")  )
Cell.type = Cell.type[4:length(Cell.type)]


# prepare ATAC reads matrix
ATAC.clean = ATAC[,  c(4:  (dim(ATAC)[2])   )]
colnames(ATAC.clean) = colnames(ATAC) [c(4:  (dim(ATAC)[2])   )]
rownames(ATAC.clean) = paste(ATAC$Chr,ATAC$Start,ATAC$End)



## Sparse matrix PCA

ATAC.clean.sparse <- (as.big.matrix(ATAC.clean)) 
x = big.PCA(ATAC.clean.sparse,pcs.to.keep = 2,verbose=T)

PC.plot.data = data.frame(PC1 = x$PCs[,1], PC2 = x$PCs[,2],Type=Cell.type)

P<-ggplot(PC.plot.data, aes(PC1,PC2,colour = Type)) 
P = P +geom_point() 
ggsave(P,filename = "PCA.jpeg")
# PCA cannot separate samples at all

