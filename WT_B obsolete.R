chooseCRANmirror()
install.packages("BiocManager")



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#install.packages("Rtools")
install.packages("edgeR")
library(edgeR)
#install.packages("extrafont")


setwd("~/Lab/Spleen_RNA_seq/star/counts/")
# read count table
WD1<-read.delim("WT_B_D1.csv",header=F)
WD2<-read.delim("WT_B_D2.csv",header=F)
WD3<-read.delim("WT_B_D3.csv",header=F)
WG1<-read.delim("WT_B_G1.csv",header=F)
WG2<-read.delim("WT_B_G2.csv",header=F)
WG3<-read.delim("WT_B_G3.csv",header=F)

# cbind all counts
y.counts<- data.frame(
  row.names=gsub("\\..*","",WD1[,1]),
  WD1=WD1[,2],WD2=WD2[,2],WD3=WD3[,2],WG1=WG1[,2],WG2=WG2[,2],WG3=WG3[,2]
)


# information about columns
y.meta <-data.frame(
  row.names=colnames(y.counts),
  treatment=c("DMSO","DMSO","DMSO","GSK","GSK","GSK"),
  #libsize=c(26713498,27068674,21813303,27818411,27275411,28664320),
  replicate=c("1","2","3","1","2","3")
  
)
#libsize=z.meta$libsize
treatment<-factor(y.meta$treatment)
replicate<-factor(y.meta$replicate)

# paired t-test model: common effect of treatment considering paired samples
design<-model.matrix(~replicate+treatment)
rownames(design)<-colnames(y.counts)
#colnames(design)<-levels(meta$group,meta$replicate)



y<-DGEList(counts=y.counts, group=treatment)
keep<- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=F]
y<- calcNormFactors(y)
y$samples
# PCA plot equivalent
plotMDS(y)
# general linear model
y.glm<-estimateDisp(y,design)
# how well does the model fit?
plotBCV(y.glm)
# determine differentially expressed genes
y.glm.fit<-glmFit(y.glm,design)
# testing the last coefficient in the linear model (GSK vs DMSO)
y.glm.lrt<-glmLRT(y.glm.fit)
topTags(y.glm.lrt)
summary(decideTests(y.glm.lrt))
# subset DGEs
y.dge<- topTags(y.glm.lrt,n=2681, sort.by="PValue")


tiff("WT_B_MD.tiff",units='in',width=5,height=5,res=1200)
y.md<-plotMD(y.glm.lrt, xlim=c(-5,15),main="GSK343 treated WT B cells")
dev.off()

################################################################################
################################################################################
# https://github.com/nerettilab/RepEnrich2/tree/d76151a4c9d42332f38e67c086f8608a668906fb
# 
# Initialize result matrices to contain the results of the GLM
results <- matrix(nrow=dim(counts)[1],ncol=0)
logfc <- matrix(nrow=dim(counts)[1],ncol=0)





res <- topTags(lrt,n=dim(c)[1],sort.by="none")$table
results <- cbind(results,res[,c(1,5)])
logfc <- cbind(logfc,res[c(1)])


# Add the repeat types back into the results.
# We should still have the same order as the input data
results$class <- WD1[,2]
results$type <- WD1[,3]



# Sort the results table by the logFC
results <- results[with(results, order(-abs(logFC))), ]
logFC <- results[, paste0("logFC")]

classes <- with(results, reorder(class, -logFC, median))

par(mar=c(6,10,4,1))
boxplot(logFC ~ classes, data=results, outline=FALSE, horizontal=TRUE,
        las=2, xlab="log(Fold Change)", main="GSK-DMSO")
abline(v=0)
# Plot the repeat types
types <- with(results, reorder(type, -logFC, median))
boxplot(logFC ~ types, data=results, outline=FALSE, horizontal=TRUE,
        las=2, xlab="log(Fold Change)", main="GSK-DMSO")
abline(v=0)
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
y.dge.fcsort <- y.dge[order(y.dge$table$logFC),]
# side column of repeat type
# match rownames of fcsorted dge table with that of the unsorted count table
# outputs a vector of positions in the count table where those rownames are found
# search for the names of repeat types at those positions from the repeat type column from WB1.txt
types<-WD1$V2[match(rownames(y.dge.fcsort),rownames(y.counts))]
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
y.dge.fcsort <- y.dge[order(y.dge$table$logFC),]

tiff(filename="WT_B_heatmap.tiff",units='in',width=7,height=6,res=600)
heatmap.2(cpm(y)[rownames(y.dge.fcsort),],
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

library(DBI)
library(AnnotationDbi)
library(GenomicFeatures)
library(GO.db)
library(org.Mm.eg.db)






#select(org.Mm.eg.db,row.names(y.counts),"ENSEMBL")
y.glm.lrt$symbol<- unname(mapIds(org.Mm.eg.db,rownames(y.glm.lrt),"SYMBOL","ENSEMBL"))
y.glm.lrt$genes<- unname(mapIds(org.Mm.eg.db,rownames(y.glm.lrt),"ENTREZID","ENSEMBL"))
y.go<-goana(y.glm.lrt,species="Mm",geneid=y.glm.lrt$genes)

topGO(y.go, sort="up", ontology = "BP", n =60)
topGO(y.go, sort="down")
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
