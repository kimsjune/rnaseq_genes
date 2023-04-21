
#########Option for adding more than one RowSideColors
#####
#####Kind of obsolete?

library(devtools)
#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#### PART A: all the upregulated genes
# change rownames of upregulated genes into ENTREZID from ENSEMBL
up.gene.entrezid<-unname(mapIds(org.Mm.eg.db,rownames(cpm(b[,c(1:12)])[b.up.combined,]),"ENTREZID","ENSEMBL"))

# blank slate that matches total number of rows (number of upregulated genes)
up.white <- rep("white",1141)

up.match.dnameth <- match(Mm.c2$REACTOME_DNA_METHYLATION ,up.gene.entrezid)
dnameth <- replace(up.white,up.match.dnameth,"darkblue")

up.match.prc2 <- match(Mm.c2$REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA,up.gene.entrezid)
prc2<- replace(up.white,up.match.prc2,"grey50")

up.match.nonsense <- match(Mm.c2$REACTOME_NONSENSE_MEDIATED_DECAY_NMD,up.gene.entrezid)
nonsense<- replace(up.white, up.match.nonsense,"red")

up.match.flu <- match(Mm.c2$REACTOME_INFLUENZA_INFECTION,up.gene.entrezid)
flu<-replace(up.white,up.match.flu,"black")

up.match.pkn1 <- match(Mm.c2$REACTOME_ACTIVATED_PKN1_STIMULATES_TRANSCRIPTION_OF_AR_ANDROGEN_RECEPTOR_REGULATED_GENES_KLK2_AND_KLK3, up.gene.entrezid)
pkn1 <- replace(up.white, up.match.pkn1, "pink")

up.match.condense <- match(Mm.c2$REACTOME_CONDENSATION_OF_PROPHASE_CHROMOSOMES, up.gene.entrezid)
condense <- replace(up.white, up.match.condense, "black")
# start with colors matching VennDiagram
sideColVennD <- c(rep("#8080FF",884),rep("#BF4080",203),rep("#FF8080",54))
# cbind all together
sideColFinal <- cbind(sideColVennD,  dnameth, nonsense)


tiff(filename="up.combined_heatmap3_c2.tiff",units='in',width=7,height=6,res=600)
heatmap.3(cpm(b[,c(1:12)])[b.up.combined,],
          col=colfunc2(12), scale="row", Rowv=NA, Colv=NA, 
          cexCol=0.7, labRow  = NA,
          density.info = "none",trace="none", dendrogram = "none",
          symkey=TRUE,symbreaks=TRUE,revC = FALSE,
          lhei=c(1.2,5), lwid=c(4, 7),margins=c(5,10),cexRow=0.5,
          rowsep=c(884,1087), RowSideColors=t(sideColVennD),
          RowSideColorsSize = 2.5)
dev.off()


#### PART B: all the downregulated genes
# change rownames of downregulated genes into ENTREZID from ENSEMBL
dn.gene.entrezid<-unname(mapIds(org.Mm.eg.db,rownames(cpm(b[,c(1:12)])[b.dn.combined,]),"ENTREZID","ENSEMBL"))

# blank slate that matches total number of rows (number of upregulated genes)
dn.white <- rep("white",1145)

dn.match.dnameth <- match(Mm.c2$REACTOME_DNA_METHYLATION ,dn.gene.entrezid)
dn.dnameth <- replace(dn.white,dn.match.dnameth,"darkblue")

dn.match.prc2 <- match(Mm.c2$REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA,dn.gene.entrezid)
dn.prc2<- replace(dn.white,dn.match.prc2,"grey50")

dn.match.nonsense <- match(Mm.c2$REACTOME_NONSENSE_MEDIATED_DECAY_NMD,dn.gene.entrezid)
dn.nonsense<- replace(dn.white, dn.match.nonsense,"red")

dn.match.flu <- match(Mm.c2$REACTOME_INFLUENZA_INFECTION,dn.gene.entrezid)
dn.flu<-replace(dn.white,dn.match.flu,"black")



# start with colors matching VennDiagram
dn.sideColVennD <- c(rep("#8080FF",922),rep("#BF4080",166),rep("#FF8080",57))
# cbind all together
dn.sideColFinal <- cbind(dn.sideColVennD,  dn.prc2,  dn.dnameth, dn.nonsense)



tiff(filename="dn.combined_heatmap3_c2.tiff",units='in',width=7,height=6,res=600)
heatmap.3(cpm(b[,c(1:12)])[b.dn.combined,],
          col=colfunc2(12), scale="row", Rowv=NA, Colv=NA, 
          cexCol=0.7, labRow  = NA,
          density.info = "none",trace="none", dendrogram = "none",
          symkey=TRUE,symbreaks=TRUE,revC = FALSE,
          lhei=c(1.2,5), lwid=c(4, 7),margins=c(5,10),cexRow=0.5,
          rowsep=c(922,1088), RowSideColors=t(dn.sideColVennD),
          RowSideColorsSize = 2.5)
dev.off()





#######








dn.symbol<-tolower(unname(mapIds(org.Mm.eg.db,rownames(cpm(b[,c(1:12)])[b.dn.combined,]),"SYMBOL","ENSEMBL")))

white <- rep("white",1147)

match.ifna <- match(gene.set.ifna,dn.symbol)
ifna <- replace(white,match.ifna,"black")

match.allograft <- match(gene.set.allograft,dn.symbol)
allograft<- replace(white,match.allograft,"black")

match.ifnresp <- match(gene.set.ifnresp,dn.symbol)
ifnresp<- replace(white,match.ifnresp,"black")

match.ifng <- match(gene.set.ifng,dn.symbol)
ifng<-replace(white,match.ifng,"black")

match.myc <- match(gene.set.myc,dn.symbol)
myc <- replace(white,match.myc,"black")

match.e2f <- match(gene.set.e2f, dn.symbol)
e2f<- replace(white,match.e2f,"black")

match.mtor<-match(gene.set.mtor,dn.symbol)
mtor<-replace(white,match.mtor,"black")
sideColOne <- c(rep("#8080FF",924),rep("#BF4080",166),rep("#FF8080",57))

sideColOne <- cbind(sideColOne,allograft,ifng,ifna)


tiff(filename="dn.combined_heatmap-final-heatmap3.tiff",units='in',width=7,height=6,res=600)
heatmap.3(cpm(b[,c(1:12)])[b.dn.combined,],
          col=colfunc2(12), scale="row", Rowv=NA, Colv=NA, 
          cexCol=0.7, labRow  = NA,
          density.info = "none",trace="none", dendrogram = "none",
          symkey=TRUE,symbreaks=TRUE,revC = FALSE,
          lhei=c(1.2,5), lwid=c(4, 7),margins=c(5,10),cexRow=0.5,
          rowsep=c(924,1090), RowSideColors=t(sideColOne),
          RowSideColorsSize = 2)
dev.off()