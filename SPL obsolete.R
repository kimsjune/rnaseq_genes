# Total splenocyte gene expression data is not shown in the paper

chooseCRANmirror()
install.packages("BiocManager")
library(DBI)
library(AnnotationDbi)
library(GenomicFeatures)
library(GO.db)
library(org.Mm.eg.db)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#install.packages("Rtools")
install.packages("edgeR")
library(edgeR)
#install.packages("extrafont")

library("edgeR")
setwd("~/Lab/Spleen_RNA_seq/star/counts/")
# read count table
SD1<-read.delim("SPL_D1.csv",header=F)
SD2<-read.delim("SPL_D2.csv",header=F)
SD3<-read.delim("SPL_D3.csv",header=F)
SG1<-read.delim("SPL_G1.csv",header=F)
SG2<-read.delim("SPL_G2.csv",header=F)
SG3<-read.delim("SPL_G3.csv",header=F)

# cbind all counts
z.counts<- data.frame(
  row.names=gsub("\\..*","",SD1[,1]),
  SD1=SD1[,2],SD2=SD2[,2],SD3=SD3[,2],SG1=SG1[,2],SG2=SG2[,2],SG3=SG3[,2]
)


# information about columns
z.meta <-data.frame(
  row.names=colnames(z.counts),
  treatment=c("DMSO","DMSO","DMSO","GSK","GSK","GSK"),
  #libsize=c(26713498,27068674,21813303,27818411,27275411,28664320),
  replicate=c("1","2","3","1","2","3")
  
)

treatment<-factor(z.meta$treatment)
replicate<-factor(z.meta$replicate)

# paired t-test model: common effect of treatment considering paired samples
design<-model.matrix(~replicate+treatment)
rownames(design)<-colnames(z.counts)
#colnames(design)<-levels(meta$group,meta$replicate)



z<-DGEList(counts=z.counts, group=treatment)
keep<- filterByExpr(z)
z <- z[keep, , keep.lib.sizes=F]
z<- calcNormFactors(z)
z$samples
# PCA plot equivalent
plotMDS(z)
# general linear model
z.glm<-estimateDisp(z,design)
# how well does the model fit?
plotBCV(z.glm)
# determine differentially expressed genes
z.glm.fit<-glmFit(z.glm,design)
# testing the last coefficient in the linear model (GSK vs DMSO)
z.glm.lrt<-glmLRT(z.glm.fit)
topTags(z.glm.lrt)
summary(decideTests(z.glm.lrt))
# subset DGEs
z.dge<- topTags(z.glm.lrt,n=405, sort.by="PValue")


tiff("SPL_MD.tiff",units='in',width=5,height=5,res=1200)
z.md<-plotMD(z.glm.lrt, xlim=c(-5,15),main="GSK343 treated WT splenocytes")
dev.off()


######################################################################################
######################################################################################
# heatmap
#install.packages("gplots")
#install.packages("RColorBrewer")
library(gplots)
library(RColorBrewer)
# coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(8))
colfunc<-colorRampPalette(c("darkblue","red"))
colfunc(10)
# sort DGEs (by P value) now by log FC
z.dge.fcsort <- z.dge[order(z.dge$table$logFC),]
# side column of repeat type
# match rownames of fcsorted dge table with that of the unsorted count table
# outputs a vector of positions in the count table where those rownames are found
# search for the names of repeat types at those positions from the repeat type column from WB1.txt
types<-SD1$V2[match(rownames(z.dge.fcsort),rownames(z.counts))]
#install.packages("stringr")
library(stringr)

types<-str_replace(types,"tRNA","pink2")
types<-str_replace(types,"LINE","green3")
types<-str_replace(types,"DNA","black")
types<-str_replace(types,"LTR","slateblue3")
types<-str_replace(types,"SINE","green3")
types<-str_replace(types,"Satellite","yellow")
types<-str_replace(types,"snRNA","white")
types<-str_replace(types,"srpRNA","white")
types<-str_replace(types,"rRNA","white")
types<-str_replace(types,"RNA","white")
types<-str_replace(types,"Other","white")










# sort DGEs (by P value) now by log FC
z.dge.fcsort <- z.dge[order(z.dge$table$logFC),]

tiff(filename="SPL_heatmap.tiff",units='in',width=7,height=6,res=600)
heatmap.2(cpm(z)[rownames(z.dge.fcsort),],
          col=colfunc(10), scale="row", Rowv=NA, Colv=NA, 
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE,
          lhei=c(1.2,5), lwid=c( 4, 10),margins=c(5,10),cexRow=1)

dev.off()

install.packages("GO.db")
install.packages("AnnotationDbi")
install.packages("org.Mm.eg.db")
install.packages("DBI")
install.packages("GenomicFeatures")
install.packages("AnnotationDbi")








z.glm.lrt$symbol<- unname(mapIds(org.Mm.eg.db,rownames(z.glm.lrt),"SYMBOL","ENSEMBL"))
z.glm.lrt$genes<- unname(mapIds(org.Mm.eg.db,rownames(z.glm.lrt),"ENTREZID","ENSEMBL"))
z.go<-goana(z.glm.lrt,species="Mm", geneid=z.glm.lrt$genes)

topGO(z.go, sort="up")
topGO(z.go, sort="down")
#####################################################################


o<-order(y.glm.lrt$table$PValue, y.glm.lrt$table$logFC)

degenes<-cpm(y)[o[1:82],]
write.table(degenes,file="WT_degenes.txt",quote=F,sep="\t")


logcpm<-cpm(y,log=TRUE,lib.size=libsize)
logcpm<-as.data.frame(logcpm)
colnames(logcpm)<-factor(meta$group)
cpm<-cpm(y,lib.size=libsize)
cpm<-as.data.frame(cpm     )
colnames(cpm)<-factor(meta$group)
write.table(cpm,file="cpm.txt",quote=F,sep="\t")


yfit<-glmQLFit(y,design)
qlf<-glmQLFTest(yfit)
topTags(qlf)
results <- matrix(nrow=dim(counts)[1],ncol=0)
logfc <- matrix(nrow=dim(counts)[1],ncol=0)

et<-exactTest(x,pair=c("DMSO","GSK"))
