# Differential gene expression analysis in wildtype B or ric mutant B cells with edgeR
# Author: Seung June Kim
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

# list of all packages needed
packages <- c("edgeR","gplots","ggplot2","ggrepel","RColorBrewer","org.Mm.eg.db","EnhancedVolcano","AnnotationDbi","VennDiagram", "DBI")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


set.seed(0)
setwd("C:/Users/jk/Documents/Lab/Spleen_RNA_seq/star/counts")

# Make a dataframe of all read count tables
# gsub removes everything after "." in ENSEMBL gene names
b.counts<- data.frame(
  row.names=gsub("\\..*","",RD1[,1]),
  WD1=WD1[,2],WD2=WD2[,2],WD3=WD3[,2],WG1=WG1[,2],WG2=WG2[,2],WG3=WG3[,2],
  RD1=RD1[,2],RD2=RD2[,2],RD3=RD3[,2],RG1=RG1[,2],RG2=RG2[,2],RG3=RG3[,2]
)

# Remove the three genes that are KO in the mutant
row.names.remove <- c("ENSMUSG00000026896","ENSMUSG00000040296","ENSMUSG00000032344")
b.counts.remove.KO.genes <- b.counts[!(row.names(b.counts)%in% row.names.remove),]

# Remove last five rows which are mapping stats
row.names.remove.mapstats <- rownames(tail(b.counts.remove.KO.genes, n=5))
b.counts.remove.KO.genes.mapstats <- b.counts.remove.KO.genes[!(row.names(b.counts.remove.KO.genes)%in% row.names.remove.mapstats),]


# sample info
b.meta <-data.frame(
  row.names=colnames(b.counts.remove.KO.genes),
  treatment=c("DMSO","DMSO","DMSO","GSK","GSK","GSK","DMSO","DMSO","DMSO","GSK","GSK","GSK"),
  replicate=c("1","2","3","1","2","3","4","5","6","4","5","6"),
  genotype=c(rep("wt",6),rep("ric",6))

)
group<- factor(paste(b.meta$treatment,b.meta$genotype,sep="."))
cbind(b.meta,group=group)
treatment<-factor(b.meta$treatment)
replicate<-factor(b.meta$replicate)
genotype<-factor(b.meta$genotype)

# paired t-test model: common effect of treatment considering paired samples
design<-model.matrix(~0+group)
colnames(design)<-levels(group)


# First step of edgeR
b<-DGEList(counts=b.counts.remove.KO.genes.mapstats, group=group)
keep<- filterByExpr(b)
b <- b[keep, , keep.lib.sizes=F]
b<- calcNormFactors(b)


# Create a PCA plot equivalent
tiff("PCA.tiff")
plotMDS(b)
dev.off()

# changing ENSEMBL gene names to to SYMBOLS
ens <- rownames(b$counts)
symbols <- mapIds(org.Mm.eg.db, keys=ens, column='SYMBOL', keytype='ENSEMBL', multiVals="filter")
symbols <- symbols[match(rownames(b$counts),names(symbols))]
rownames(b$counts)<-symbols
notNA <- !is.na(rownames(b$counts))
b <- b[notNA,]

# general linear model
b.glm<-estimateDisp(b,design)

# how well does the model fit?
plotBCV(b.glm)


# Determine DGEs using Likelihood Ratio test
b.glm.fit<-glmFit(b.glm,design)

my.contrasts<- makeContrasts(
  wt.GvD= GSK.wt-DMSO.wt,
  ric.GvD= GSK.ric-DMSO.ric,
  GSK.wtvric=GSK.wt-GSK.ric,
  DMSO.wtvric=DMSO.wt-DMSO.ric,
  levels=design
)

# DGE from a given comparison
b.gsk.wt.ric.lrt <- glmLRT(b.glm.fit, contrast=my.contrasts[,"GSK.wtvric"])
b.dmso.wt.ric.lrt <- glmLRT(b.glm.fit, contrast=my.contrasts[,"DMSO.wtvric"])
b.wt.gskvdmso.lrt <- glmLRT(b.glm.fit,contrast=my.contrasts[,"wt.GvD"])
b.ric.gskvdmso.lrt <- glmLRT(b.glm.fit,contrast=my.contrasts[,"ric.GvD"])




# plot MD
tiff("combined_wt_DMSOvGSK.tiff",units='in',width=5,height=5,res=1200)
b.wt.DvG.md<-plotMD(b.wt.gskvdmso.lrt, xlim=c(-5,15),cex=0.5,main=NULL)
dev.off()
tiff("combined_ric_DMSOvGSK.tiff",units='in',width=5,height=5,res=1200)
b.mut.DvG.md<-plotMD(b.ric.gskvdmso.lrt, xlim=c(-5,15),cex=0.5,main=NULL)
dev.off()




# find out the exact number of DGEs and write to an object
summary(decideTests(b.wt.gskvdmso.lrt))
b.wt.gskvdmso.dge <- topTags(b.wt.gskvdmso.lrt, n=1871, sort.by="PValue")
b.wt.gskvdmso.dge.fcsort <- b.wt.gskvdmso.dge[order(b.wt.gskvdmso.dge$table$logFC,decreasing=TRUE),]
summary(decideTests(b.ric.gskvdmso.lrt))
b.ric.gskvdmso.dge <- topTags(b.ric.gskvdmso.lrt, n=595,sort.by="PValue" )
b.ric.gskvdmso.dge.fcsort <- b.ric.gskvdmso.dge[order(b.ric.gskvdmso.dge$table$logFC,decreasing=T),]

# subset positive FC in wt
b.wt.gskvdmso.dge.up<-b.wt.gskvdmso.dge.fcsort[b.wt.gskvdmso.dge.fcsort$table$logFC>0,]
# subset positive FC in mutant
b.ric.gskvdmso.dge.up <- b.ric.gskvdmso.dge.fcsort[b.ric.gskvdmso.dge.fcsort$table$logFC>0,]
# write rownames of each slice of the VennDiagram -- unique to wt, unique to Mut, intersect
b.up.in.wt.not.ric <- setdiff(rownames(b.wt.gskvdmso.dge.up),rownames(b.ric.gskvdmso.dge.up))
b.up.in.ric.not.wt <- setdiff(rownames(b.ric.gskvdmso.dge.up),rownames(b.wt.gskvdmso.dge.fcsort))
b.up.in.both <- intersect(rownames(b.wt.gskvdmso.dge.up),rownames(b.ric.gskvdmso.dge.up))
# negative FC in wt
b.wt.gskvdmso.dge.dn<-b.wt.gskvdmso.dge.fcsort[b.wt.gskvdmso.dge.fcsort$table$logFC<0,]
# negative FC in mutant
b.ric.gskvdmso.dge.dn <- b.ric.gskvdmso.dge.fcsort[b.ric.gskvdmso.dge.fcsort$table$logFC<0,]
# write rownames of each slice of the VennDiagram -- unique to wt, unique to Mut, intersect
b.dn.in.wt.not.ric <- setdiff(rownames(b.wt.gskvdmso.dge.dn),rownames(b.ric.gskvdmso.dge.dn))
b.dn.in.ric.not.wt <- setdiff(rownames(b.ric.gskvdmso.dge.dn),rownames(b.wt.gskvdmso.dge.fcsort))
b.dn.in.both <- intersect(rownames(b.wt.gskvdmso.dge.dn),rownames(b.ric.gskvdmso.dge.dn))



# up or down regulated genes together and add rowsep ## 3 colors
b.up.combined <- c(b.up.in.wt.not.ric, b.up.in.both, b.up.in.ric.not.wt)
b.dn.combined <- c(b.dn.in.wt.not.ric, b.dn.in.both, b.dn.in.ric.not.wt)


#### Change ENSEMBL to ENTREZID for CAMERA/GSEA
rownames(b.glm.fit) <- unname(mapIds(org.Mm.eg.db, keys=rownames(b.glm.fit$counts), keytype = "SYMBOL", column="ENTREZID"))
# some ENSEMBL gene IDs are not mapped to ENTREZ ID, so they are NA. 
# These have to be removed.
omitNArows <- na.omit(rownames(b.glm.fit))
b.glm.fit.noNA <- b.glm.fit[omitNArows,]

#rownames(b.wt.gskvdmso.lrt) <- unname(mapIds(org.Mm.eg.db, keys=rownames(b.wt.gskvdmso.lrt), keytype = "ENSEMBL", column="ENTREZID"))


#####

# For GSEA with CAMERA, go to c*.R 
# https://www.bioconductor.org/packages/devel/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html#camera-gene-set-enrichment-analysis

# All of the comparisons
GSK.wtvric<- makeContrasts(GSK.wt-GSK.ric, levels=design)
GSK.ricvwt<- makeContrasts(GSK.ric-GSK.wt, levels=design)
DMSO.wtvric<- makeContrasts(DMSO.wt-DMSO.ric, levels=design)
wt.GSKvDMSO <- makeContrasts(GSK.wt-DMSO.wt, levels=design)
ric.GSKvDMSO <- makeContrasts(GSK.ric-DMSO.ric, levels=design)


### VennDiagram of common vs. uniquely up/downregulated genes
# commonly up and downregulated genes upon GSK343 treatment
#wt.gvd.lrt<-glmLRT(b.glm.fit.noNA, contrast=my.contrasts[,"wt.GvD"])
dif.wt.gvd<-topTags(b.wt.gskvdmso.lrt,n=Inf,p=0.05)$table
up.dif.wt.gvd <- row.names(dif.wt.gvd[dif.wt.gvd$logFC>0,])
down.dif.wt.gvd <- row.names(dif.wt.gvd[dif.wt.gvd$logFC<0,])

#ric.gvd.lrt<-glmLRT(b.glm.fit.noNA, contrast=my.contrasts[,"ric.GvD"])
dif.ric.gvd<-topTags(b.ric.gskvdmso.lrt,n=Inf,p=0.05)$table
up.dif.ric.gvd <- row.names(dif.ric.gvd[dif.ric.gvd$logFC>0,])
down.dif.ric.gvd <- row.names(dif.ric.gvd[dif.ric.gvd$logFC<0,])

venn.up<- venn.diagram(
  
  x=list(up.dif.wt.gvd,up.dif.ric.gvd),
  category.names=c("",""),
  #resolution=600,
  filename= NULL,
  euler.d=T,
  output=T,
  scaled=T,
  cex=2.5,
  lwd=4,
  col="transparent",
 
  fill=c("blue","red"),
  ext.dist=0.05,
  main.fontfamily = 'sans'
)
ggsave(venn.up, file="venn_up.svg",device="svg",height=4,width=4) 
dev.off()

venn.down <- venn.diagram(
  x=list(down.dif.wt.gvd,down.dif.ric.gvd),
  category.names=c("",""),
  #resolution=600,
  filename= NULL,
  euler.d=T,
  output=T,
  scaled=T,
  cex=2.5,
  lwd=4,
  col="transparent",
  
  fill=c("blue","red"),
  ext.dist=0.05,
  main.fontfamily='sans'
)

ggsave(venn.down, file="venn_down.svg",device="svg",height=4,width=4)
dev.off()






# Creating files for GSEA (Broad) 
#####
# gene enrichment analysis uses the whole normalized count table, not just the DGE list

cpm_gsk_wtvric<-cpm(b[,c(4,5,6,10,11,12)])
#rownames(table)<- unname(mapIds(org.Mm.eg.db,rownames(table),"SYMBOL","ENSEMBL"))
write.table(cpm_gsk_wtvric,file="cpm_gsk_wtvric.txt",quote=F,sep="\t")
cpm_wt_gskvdmso<-cpm(b[,c(1,2,3,4,5,6)])
write.table(cpm_wt_gskvdmso,file="cpm_wt_gskvdmso.txt",quote=F,sep="\t")
cpm_ric_gskvdmso<-cpm(b[,c(7,8,9,10,11,12)])
write.table(cpm_ric_gskvdmso,file="cpm_ric_gskvdmso.txt",quote=F,sep="\t")

cpm_dmso_wtvric<-cpm(b[,c(1,2,3,7,8,9)])
write.table(cpm_dmso_wtvric,file="cpm_dmso_wtvric.txt",quote=F,sep="\t")
#####


## Expression of Cgas, Ifih1 and Ddx58 
b.target.not.removed<-DGEList(counts=b.counts, group=group)
b.target.not.removed<- calcNormFactors(b.target.not.removed)
cpm(b.target.not.removed)[c("ENSMUSG00000026896","ENSMUSG00000040296","ENSMUSG00000032344"),]

# fold-change of monocyte chemotaxis genes

monocyte <- c("Ccl22","Ccl3","Ccl5","Pdgfb","Ccl1","Cx3cl1",
              "Slamf8","Lgals3","Tnfrsf11a","Ccl6","Lgmn","Ccr1")

b.wt.gskvdmso.lrt$table[monocyte,]


# Volcano plots

res <-(topTags(b.wt.gskvdmso.lrt,n="inf", sort.by="PValue"))
#labs.wt <- topTags(labs.wt, sort.by="PValue",n=10)
source("C:/users/jk/Documents/Lab/EnhancedVolcano.R")
library(ggplot2)
library(ggrepel)
svg("WT_volcano.svg", family='sans', width=7, height=5)

EnhancedVolcano(data.frame(res)
                ,
                lab=rownames(res),
                selectLab= c("Ahr","Tagln2","Ahrr",
                             "Rgs10","Flrt3","Scin","Prkar1b",
                             "Cebpa","Nme4"),
                x='logFC',
                y='FDR',
                # xlim=c(-2.5,2.5),
                # ylim=c(0,5),
                
                pCutoff = 0.05,
                FCcutoff=0,
                gridlines.major=F,
                gridlines.minor=F,
                
                legendPosition='none',
                title=NULL,
                subtitle=NULL,
                labSize=9,
                axisLabSize=24,
                labFace='bold',
                drawConnectors=T,
                widthConnectors=0.5,
                col=c('grey','grey','grey','red'),
                colAlpha=1,
                shadeAlpha=1,
                xlab =" ",
                ylab= " ",
                caption=NULL,
                boxedLabels=F)
dev.off()                





r.res <-(topTags(b.rlr.gskvdmso.lrt,n="inf", sort.by="PValue"))
#labs.wt <- topTags(labs.wt, sort.by="PValue",n=10)
source("C:/users/jk/Documents/Lab/EnhancedVolcano.R")
library(ggplot2)
svg("RLR_volcano.svg", family='sans', width=7, height=5)

EnhancedVolcano(data.frame(r.res)
                ,
                lab=rownames(r.res),
                selectLab= c("Ahr","Tagln2","Ahrr",
                             "Rgs10","Flrt3","Scin","Prkar1b",
                             "Cebpa","Nme4"),
                x='logFC',
                y='FDR',
                # xlim=c(-2.5,2.5),
                ylim=c(0,15),
                
                pCutoff = 0.05,
                FCcutoff=0,
                gridlines.major=F,
                gridlines.minor=F,
                
                legendPosition='none',
                title=NULL,
                subtitle=NULL,
                labSize=9,
                axisLabSize=24,
                labFace='bold',
                drawConnectors=T,
                widthConnectors=0.5,
                col=c('grey','grey','grey','red'),
                colAlpha=1,
                shadeAlpha=1,
                xlab =" ",
                ylab= " ",
                caption=NULL,
                boxedLabels=F)
dev.off()               

# heatmaps
# 
heatmap.2(cpm(b)[(b.up.in.wt.not.ric),],
          col=colfunc(10), scale="row", Rowv=NA, Colv=NA, 
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE,
          lhei=c(1.2,5), lwid=c( 4, 10),margins=c(5,10),cexRow=1)
