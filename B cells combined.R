# edgeR: analyzing WT and RLR (triple mutant) B cells together
# Author: Seung June Kim
library(edgeR)
library(gplots)
library(RColorBrewer)
library(org.Mm.eg.db)
library(GO.db)
library(AnnotationDbi)
library(DBI)
library(VennDiagram)
library(EnhancedVolcano)
source("./replotGSEA.R")

setwd("C:/Users/jk/Documents/Lab/Spleen_RNA_seq/star/counts")

# Make a dataframe of all read count tables
b.counts<- data.frame(
  row.names=gsub("\\..*","",RD1[,1]),
  WD1=WD1[,2],WD2=WD2[,2],WD3=WD3[,2],WG1=WG1[,2],WG2=WG2[,2],WG3=WG3[,2],
  RD1=RD1[,2],RD2=RD2[,2],RD3=RD3[,2],RG1=RG1[,2],RG2=RG2[,2],RG3=RG3[,2]
)
# Remove the three genes that are KO in the mutant
row.names.remove <- c("ENSMUSG00000026896","ENSMUSG00000040296","ENSMUSG00000032344")
b.counts.remove.KO.genes <- b.counts[!(row.names(b.counts)%in% row.names.remove),]
b.counts.KO.genes <-  b.counts[(row.names(b.counts)%in% row.names.remove),]



# ENSEMBL is kind of useless with downstream funcitons. Here they are converted
# into SYMBOL for Volcano plots. Convert into ENTREZID for gene set enrichment.


# sample info
b.meta <-data.frame(
  row.names=colnames(b.counts.remove.KO.genes),
  treatment=c("DMSO","DMSO","DMSO","GSK","GSK","GSK","DMSO","DMSO","DMSO","GSK","GSK","GSK"),
  replicate=c("1","2","3","1","2","3","4","5","6","4","5","6"),
  genotype=c(rep("WT",6),rep("RLR",6))

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
b<-DGEList(counts=b.counts.remove.KO.genes, group=group)
keep<- filterByExpr(b)
b <- b[keep, , keep.lib.sizes=F]
b<- calcNormFactors(b)


# Create a PCA plot equivalent
tiff("PCA.tiff")
plotMDS(b)
dev.off()

# changing ENSEMBL to SYMBOLS
# makes you wonder why I bothered to use ENSEMBL in the first place
ens <- rownames(b$counts)
symbols <- mapIds(org.Mm.eg.db, keys=ens, column='SYMBOL', keytype='ENSEMBL',multiVals="filter")
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
  WT.GvD= GSK.WT-DMSO.WT,
  RLR.GvD= GSK.RLR-DMSO.RLR,
  GSK.WTvRLR=GSK.WT-GSK.RLR,
  DMSO.WTvRLR=DMSO.WT-DMSO.RLR,
  levels=design
)

# DGE from a given comparison
b.gsk.wt.rlr.lrt <- glmLRT(b.glm.fit, contrast=my.contrasts[,"GSK.WTvRLR"])
b.dmso.wt.rlr.lrt <- glmLRT(b.glm.fit, contrast=my.contrasts[,"DMSO.WTvRLR"])
b.wt.gskvdmso.lrt <- glmLRT(b.glm.fit,contrast=my.contrasts[,"WT.GvD"])
b.rlr.gskvdmso.lrt <- glmLRT(b.glm.fit,contrast=my.contrasts[,"RLR.GvD"])




# plot MD
tiff("combined_WT_DMSOvGSK.tiff",units='in',width=5,height=5,res=1200)
b.wt.DvG.md<-plotMD(b.wt.gskvdmso.lrt, xlim=c(-5,15),cex=0.5,main=NULL)
dev.off()
tiff("combined_RLR_DMSOvGSK.tiff",units='in',width=5,height=5,res=1200)
b.mut.DvG.md<-plotMD(b.rlr.gskvdmso.lrt, xlim=c(-5,15),cex=0.5,main=NULL)
dev.off()




# find out the exact number of DGEs and write to an object
summary(decideTests(b.wt.gskvdmso.lrt))
b.wt.gskvdmso.dge <- topTags(b.wt.gskvdmso.lrt, n=1785, sort.by="PValue")
b.wt.gskvdmso.dge.fcsort <- b.wt.gskvdmso.dge[order(b.wt.gskvdmso.dge$table$logFC,decreasing=TRUE),]
summary(decideTests(b.rlr.gskvdmso.lrt))
b.rlr.gskvdmso.dge <- topTags(b.rlr.gskvdmso.lrt, n=581,sort.by="PValue" )
b.rlr.gskvdmso.dge.fcsort <- b.rlr.gskvdmso.dge[order(b.rlr.gskvdmso.dge$table$logFC,decreasing=T),]

# subset positive FC in WT
b.wt.gskvdmso.dge.up<-b.wt.gskvdmso.dge.fcsort[b.wt.gskvdmso.dge.fcsort$table$logFC>0,]
# subset positive FC in mutant
b.rlr.gskvdmso.dge.up <- b.rlr.gskvdmso.dge.fcsort[b.rlr.gskvdmso.dge.fcsort$table$logFC>0,]
# write rownames of each slice of the VennDiagram -- unique to WT, unique to Mut, intersect
b.up.in.wt.not.rlr <- setdiff(rownames(b.wt.gskvdmso.dge.up),rownames(b.rlr.gskvdmso.dge.up))
b.up.in.rlr.not.wt <- setdiff(rownames(b.rlr.gskvdmso.dge.up),rownames(b.wt.gskvdmso.dge.fcsort))
b.up.in.both <- intersect(rownames(b.wt.gskvdmso.dge.up),rownames(b.rlr.gskvdmso.dge.up))
# negative FC in WT
b.wt.gskvdmso.dge.dn<-b.wt.gskvdmso.dge.fcsort[b.wt.gskvdmso.dge.fcsort$table$logFC<0,]
# negative FC in mutant
b.rlr.gskvdmso.dge.dn <- b.rlr.gskvdmso.dge.fcsort[b.rlr.gskvdmso.dge.fcsort$table$logFC<0,]
# write rownames of each slice of the VennDiagram -- unique to WT, unique to Mut, intersect
b.dn.in.wt.not.rlr <- setdiff(rownames(b.wt.gskvdmso.dge.dn),rownames(b.rlr.gskvdmso.dge.dn))
b.dn.in.rlr.not.wt <- setdiff(rownames(b.rlr.gskvdmso.dge.dn),rownames(b.wt.gskvdmso.dge.fcsort))
b.dn.in.both <- intersect(rownames(b.wt.gskvdmso.dge.dn),rownames(b.rlr.gskvdmso.dge.dn))



# up or down regulated genes together and add rowsep ## 3 colors
b.up.combined <- c(b.up.in.wt.not.rlr, b.up.in.both, b.up.in.rlr.not.wt)
b.dn.combined <- c(b.dn.in.wt.not.rlr, b.dn.in.both, b.dn.in.rlr.not.wt)

# set colour
colfunc<-colorRampPalette(c("blue","black","red"))
blank <- rep("white",1925)
combined.symbols <- rownames(cpm(b[,c(1:12)])[c(b.up.combined,b.dn.combined),])
match <-  match(g.IL6.symbol ,combined.symbols)
sidecol <- replace(blank,match,"blue")
# plot both up and down together!
svglite::svglite(file="combined_heatmap-finalwsidecol.svg",width=6.5,height=9,
                 system_fonts='Arial',bg='white')
heatmap.2(cpm(b[,c(1:12)])[c(b.up.combined,b.dn.combined),],
          col=colfunc(8), scale="row", Rowv=NA, Colv=NA, 
          cexCol=1,srtCol=0, adjCol=c(0.5,0), labRow  = NA, labCol=NA,
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, 
          RowSideColors=sidecol, lmat=rbind(c(5,0,4),c(3,1,2)),
          lhei=c(1,7), lwid=c( 5.5, 0.5,10),margins=c(1,1),cexRow=1,
          rowsep=c(883,1145,1229,1692,1870))
dev.off()

# heatmap of cytokines/chemokines that were not filtered out... 
IL6 <- c("Cxcl13","Csf1","Csf3r","Hmox1","Cd9","Socs1","Cd44","Fas","Ltbr",
         "Csf2ra","Ebi3","Bak1","Tnfrsf1b","Il2ra","Tnfrsf1a",
        "Il6st","Csf2rb","Il13ra1","Tnf","Il3ra","Ccr1")

ribosome <- c("Rpl23a","Rpl28","Rps12","Rplp2","Rps23",
              "Rpl36","Rps27a","Rpl22l1","Rpl14","Rpl35a",
              "Rps16","Rpl19",
              "Rpsa","Rps17","Rps11","Rpl35","Rps18","Rplp1",
              "Rpl10","Rpl10a","Ybx1","Pabpc1","Eif3d")
cytochemo <- c("Cx3cl1","Ccl22","Ccl3","Ccl5","Ccl1","Pla2g7","Pdgfb",
              "Slamf8","Lgals3","Tnfrsf11a","Dusp1","Rps19",
               "Lgmn","Ccr1")
t1inf <- c("Adar","Cactin","Cdc37","Cnot7","Dcst1","Fadd",
           "Hsp90ab1","Ifi27","Ifitm1","Ifitm2","Ifnar1",
           "Ifnar2","Ikbke","Irak1","Irf3","Irf7","Isg15",
           "Jak1","Lsm14a","Mavs","Mettl3","Mmp12","Mul1","Mx1",
           "Myd88","Nlrc5","Oas2","Oas3","Ptpn1","Ptpn11",
           "Ptpn2","Ptpn6","Rnf185","Samhd1","Shfl","Shmt2","Smpd1",
           "Sp100","Stat1","Stat2","Tbk1","Trim41",
           "Trim56","Ttll12","Tyk2","Ube2k","Usp18","Usp27x",
           "Ythdf2","Ythdf3","Zbp1")
t1infshort <- c(
           "Hsp90ab1","Ifi27","Ifitm2","Ifnar1",
           "Ifnar2","Ikbke","Irak1","Irf3","Irf7","Isg15",
           "Mx1",
           "Myd88","Nlrc5","Oas2","Oas3","Ptpn1","Ptpn11",
           "Ptpn2","Ptpn6","Rnf185","Samhd1","Shfl","Shmt2","Smpd1",
           "Sp100","Stat1","Stat2","Tbk1","Trim41",
           "Trim56","Ttll12","Tyk2","Ube2k","Usp18","Usp27x",
           "Ythdf2","Ythdf3","Zbp1")
infla <- c("blue","blue",rep("white",16),"blue","white","white", rep("blue",5),"white","white","white",rep("white",10)
           )
#svglite::svglite(file="gsea_heatmap.svg",width=8,height=7,
           #      system_fonts='Arial',bg='white')
pdf(file="gsea_heatmap.pdf")
par(font=2, font.axis=2, font.lab=2,cex.lab=1.5, cex.axis=1.5,
    lty=1)
heatmap.2(cpm(b[,c(1:12)])[as.character(c(ribosome,cytochemo)),],
          col=colfunc(100), scale="row", Rowv=FALSE, Colv=FALSE, 
          cexCol=1,srtCol=0, adjCol=c(0.5,0), labCol=NA,
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, 
          RowSideColors=c(rep("white",23),c(rep("blue",5),rep("white",9))), 
          lmat=rbind(c(5,0,4),c(3,1,2)),
          lhei=c(1,4), lwid=c( 3, 0.5,8),margins=c(5,7),
          cexRow=1.2,
          colsep=c(6),
          rowsep=c(23))
dev.off()

# what about another cytokine?
heatmap.2(cpm(b[,c(1:12)])[as.character(c(t1infshort)),],
          col=colfunc(8), scale="row", Rowv=TRUE
          , Colv=FALSE, 
          cexCol=1,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors=infla, 
          lmat=rbind(c(5,0,4),c(3,1,2)),
          lhei=c(2,4), lwid=c( 5, 0.5,10),margins=c(5,7),cexRow=0.6,
          colsep=c(6))




# Volcano plots
# show IL6-STAT gene set genes?
IL6 <- read.delim("./../gsea_hallmark_svg/WT_gvd_IL6_geneset.txt",sep='\t')
E2F <- read.delim("./../gsea_hallmark_svg//WT_gvd_E2F_geneset.txt",sep='\t')
library(stringr)
IL6.symbol <- str_to_title(IL6$SYMBOL)
E2F.symbol <- str_to_title(E2F$SYMBOL)
w.res <- topTags(b.wt.gskvdmso.lrt, sort.by="logFC", n="Inf")
w.lab <- rownames(topTags(b.wt.gskvdmso.lrt,sort.by="PVal",n=50))
svg("WT_volcano.svg", family='sans', width=11, height=7)
EnhancedVolcano(data.frame(w.res),
                lab=rownames(w.res),
                x='logFC',
                y='FDR',
                selectLab=c("Ccl22"),
                pCutoff=0.05,
                FCcutoff=F,
                gridlines.minor=F,
                gridlines.major=F,
                legendPosition='none',
                ylab="-logFDR",
                title=NULL,
                subtitle=NULL,
                labSize=4,
                xlim=c(-2.25,2.25),
                ylim=c(0,8),
                axisLabSize=18,
                drawConnectors=T,
                col=c('gray','gray','gray','red'),
                colAlpha=0.3
)
dev.off()

r.res <- topTags(b.rlr.gskvdmso.lrt, sort.by="logFC", n="Inf")
#r.lab <- rownames(topTags(b.rlr.gskvdmso.lrt,sort.by="PVal",n=50))
svg("Mut_volcano.svg", family='sans', width=11, height=7)
EnhancedVolcano(data.frame(r.res),
                lab=rownames(r.res),
                x='logFC',
                y='FDR',
                selectLab=c("Ccl22"),
                pCutoff=0.05,
                FCcutoff=F,
                gridlines.minor=F,
                gridlines.major=F,
                legendPosition='none',
                ylab="-logFDR",
                title=NULL,
                subtitle=NULL,
                labSize=4,
                xlim=c(-2.25,2.25),
                ylim=c(0,8),
                axisLabSize=18,
                drawConnectors=F,
                col=c('gray','gray','gray','red'),
                colAlpha=0.3
)
dev.off()
g.IL6 <- read.delim("./../gsea_hallmark_svg/gsk_wvr_Il6.txt",sep='\t')
g.tnf <- read.delim("./../gsea_hallmark_svg/gsk_wvr_TNF.txt",sep='\t')

g.IL6.symbol <- str_to_title(g.IL6$SYMBOL)
g.tnf.symbol <- str_to_title(g.tnf$SYMBOL)
g.res <- topTags(b.gsk.wt.rlr.lrt, sort.by="logFC", n="Inf")
#r.lab <- rownames(topTags(b.rlr.gskvdmso.lrt,sort.by="PVal",n=50))
svg("GSK_volcano.svg", family='sans', width=11, height=7)
EnhancedVolcano(data.frame(g.res),
                lab=rownames(g.res),
                x='logFC',
                y='FDR',
                selectLab=g.tnf.symbol,
                pCutoff=0.05,
                FCcutoff=F,
                gridlines.minor=F,
                gridlines.major=F,
                legendPosition='none',
                ylab="-logFDR",
                title=NULL,
                subtitle=NULL,
                labSize=3,
                xlim=c(-1,1),
               ylim=c(0,15),
                axisLabSize=18,
                drawConnectors=F,
                col=c('gray','gray','gray','red'),
                colAlpha=0.3
)
dev.off()

## deprecated - down regulated genes separately
#####
# dn regulated genes together and add rowsep 3 colors
b.dn.combined <- c(b.dn.in.wt.not.rlr, b.dn.in.both, b.dn.in.rlr.not.wt)
svglite::svglite(file="dn.combined_heatmap-final.svg",width=7,height=6,
                 system_fonts='Arial',bg='white')
heatmap.2(cpm(b[,c(1:12)])[b.dn.combined,],
          col=colfunc2(8), scale="row", Rowv=NA, Colv=NA, 
          cexCol=1,srtCol=0, adjCol=c(0.5,0), labRow  = NA,
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE,
          lmat=rbind(c(5,0,4),c(3,1,2)),RowSideColors = sidecol,
          lhei=c(1.2,5), lwid=c( 2,0.2, 5),margins=c(1,1),cexRow=1,
          rowsep=c(464,642))
dev.off()
#####
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
GSK.WTvRLR<- makeContrasts(GSK.WT-GSK.RLR, levels=design)
GSK.RLRvWT<- makeContrasts(GSK.RLR-GSK.WT, levels=design)
DMSO.WTvRLR<- makeContrasts(DMSO.WT-DMSO.RLR, levels=design)
WT.GSKvDMSO <- makeContrasts(GSK.WT-DMSO.WT, levels=design)
RLR.GSKvDMSO <- makeContrasts(GSK.RLR-DMSO.RLR, levels=design)


### VennDiagram of common vs. uniquely up/downregulated genes
# commonly up and downregulated genes upon GSK343 treatment
#wt.gvd.lrt<-glmLRT(b.glm.fit.noNA, contrast=my.contrasts[,"WT.GvD"])
dif.wt.gvd<-topTags(b.wt.gskvdmso.lrt,n=Inf,p=0.05)$table
up.dif.wt.gvd <- row.names(dif.wt.gvd[dif.wt.gvd$logFC>0,])
down.dif.wt.gvd <- row.names(dif.wt.gvd[dif.wt.gvd$logFC<0,])

#rlr.gvd.lrt<-glmLRT(b.glm.fit.noNA, contrast=my.contrasts[,"RLR.GvD"])
dif.rlr.gvd<-topTags(b.rlr.gskvdmso.lrt,n=Inf,p=0.05)$table
up.dif.rlr.gvd <- row.names(dif.rlr.gvd[dif.rlr.gvd$logFC>0,])
down.dif.rlr.gvd <- row.names(dif.rlr.gvd[dif.rlr.gvd$logFC<0,])

venn.up<- venn.diagram(
  
  x=list(up.dif.wt.gvd,up.dif.rlr.gvd),
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
  x=list(down.dif.wt.gvd,down.dif.rlr.gvd),
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
# a revelation: gene enrichment analysis uses the whole normalized count table, not just the DGE list

cpm_gsk_wtvrlr<-cpm(b[,c(4,5,6,10,11,12)])
#rownames(table)<- unname(mapIds(org.Mm.eg.db,rownames(table),"SYMBOL","ENSEMBL"))
write.table(cpm_gsk_wtvrlr,file="cpm_gsk_wtvrlr.txt",quote=F,sep="\t")
cpm_wt_gskvdmso<-cpm(b[,c(1,2,3,4,5,6)])
write.table(cpm_wt_gskvdmso,file="cpm_wt_gskvdmso.txt",quote=F,sep="\t")
cpm_rlr_gskvdmso<-cpm(b[,c(7,8,9,10,11,12)])
write.table(cpm_rlr_gskvdmso,file="cpm_rlr_gskvdmso.txt",quote=F,sep="\t")

cpm_dmso_wtvrlr<-cpm(b[,c(1,2,3,7,8,9)])
write.table(cpm_dmso_wtvrlr,file="cpm_dmso_wtvrlr.txt",quote=F,sep="\t")
#####


## Expression of Cgas, Ifih1 and Ddx58 
b.target.not.removed<-DGEList(counts=b.counts, group=group)
b.target.not.removed<- calcNormFactors(b.target.not.removed)
cpm(b.target.not.removed)[c("ENSMUSG00000026896","ENSMUSG00000040296","ENSMUSG00000032344"),]

# fold-change of monocyte chemotaxis genes

monocyte <- c("Ccl22","Ccl3","Ccl5","Pdgfb","Ccl1","Cx3cl1",
              "Slamf8","Lgals3","Tnfrsf11a","Ccl6","Lgmn","Ccr1")

b.wt.gskvdmso.lrt$table[monocyte,]
