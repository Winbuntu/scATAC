library(GenomicRanges)
library(data.table)
library(FactoMineR)
library(ggplot2)
####

Serum.read.count = read.table("simutaneous/GSE74534_RNA-seq_normalized_counts.txt",row.names = 1,header = T)
Serum.read.count.filtered = Serum.read.count[apply(Serum.read.count,1,function(x) sum(x>0))>=6,] # 10% cells express a gene, then keep
write.csv(Serum.read.count.filtered,file = "Filtered.normlaized.expression.csv",quote = F)


# first PCA with expressed genes
pca.res = PCA( t(Serum.read.count.filtered),graph = F)

PC.plot.data = data.frame(PC1 = as.numeric(pca.res$ind$coord[,1]), 
                          PC2 = as.numeric(pca.res$ind$coord[,2]))
P<-ggplot(PC.plot.data, aes(PC1,PC2)) 
P = P +geom_point() 
print(P)
save(P,file = "First PCA with expressed genes.jpeg")


####################################################################################################################################

# here we finished locally linear embedding analysis in python.
# compute promoter level methylation

# first extract promoter regions into a GenomicRange object


library(GenomicFeatures)

txdb <- makeTxDbFromGFF("genes.gtf", format="gtf") # mm10 genome annotation

mm10.genes <- genes(txdb)

mm10.promoters.500bp <- promoters(mm10.genes,upstream = 500,downstream =0)

save(mm10.promoters.500bp,file = "mm10.promoters.500bp.RData")
#write.table(as.data.frame(genes)[,-4], file="Just_genes.txt", colnames=F, sep="\t")

####################################################################################################################################

mESC_exp_SamFil_GenFil = read.csv("mESC_exp_SamFil_GenFil.csv",header = T,row.names = 1)



## process sample CpG methylation information # loop starts at here

sample.name = substr(colnames(mESC_exp_SamFil_GenFil),13,100)
WGBS.name = paste(sample.name,".CpG.txt.gz",sep ="")
Promoter.meth.big.table = data.frame(order = c(1:37334))



for(one.WGBS.lib in WGBS.name) {
#one.WGBS.lib = "A02.CpG.txt.gz"
system(paste("gzip -d ./simutaneous/GSE68642_RAW/" ,one.WGBS.lib,sep =""))

meth.info = fread(paste("./simutaneous/GSE68642_RAW/" ,substr(one.WGBS.lib,1,11),sep =""))

# <chromosome> <position> <strand> <count methylated> <count non-methylated> <C-context> <trinucleotide context>

meth.info.GR <- GRanges(
  
       seqnames = meth.info$V1,
       ranges = IRanges(start = meth.info$V2,end = meth.info$V2),
       strand = meth.info$V3,
       CpG.coverage = (meth.info$V4 + meth.info$V5),
       Meth = meth.info$V4/(meth.info$V4 + meth.info$V5)
)

overlapped.sites = findOverlaps(meth.info.GR,mm10.promoters.500bp) # query, subject

counts = as.data.frame(table(mcols(mm10.promoters.500bp)$gene_id[subjectHits(overlapped.sites)]))
hist(counts$Freq,main = "CpG Frequency on Promoter (500bp Upstream)") # frequency of CpG on promoters

promoter.CpGs = data.frame(gene = mcols(mm10.promoters.500bp)$gene_id[subjectHits(overlapped.sites)],
           CpG.Meth = mcols(meth.info.GR)$Meth[queryHits(overlapped.sites)])

promoter.meth.level = aggregate(x = promoter.CpGs$CpG.Meth, by = list(promoter.CpGs$gene), FUN = mean , na.rm = T)

Promoter.meth.big.table = cbind(Promoter.meth.big.table, promoter.meth.level)

gc() # collect garbale to free memory

system(paste("gzip ./simutaneous/GSE68642_RAW/" ,substr(one.WGBS.lib,1,11) ,sep =""))

}
