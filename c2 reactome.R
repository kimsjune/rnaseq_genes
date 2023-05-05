# requires variables from B cells combined.R

Mm.c2 <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c2.cp.reactome.v7.1.entrez.rds"))
# from B cells combined.R, b.glm.fit.rownames already has ENTREZID swapped for ENSEMBL
idx.c2 <- ids2indices(Mm.c2,id=rownames(b.glm.fit))

cam.c2.GSK.WTvRLR <- camera(b.glm, idx.c2, design, contrast=GSK.WTvRLR)
cam.c2.GSK.RLRvWT<- camera(b.glm,idx.c2,design,contrast=GSK.RLRvWT)
cam.c2.DMSO.WTvRLR<- camera(b.glm,idx.c2,design, contrast=DMSO.WTvRLR)
cam.c2.WT.GSKvDMSO <- camera(b.glm,idx.c2,design, contrast=WT.GSKvDMSO)
cam.c2.RLR.GSKvDMSO <- camera(b.glm,idx.c2,design, contrast=RLR.GSKvDMSO)
cam.c2.RLR.GSKvDMSO2 <- camera(b.glm, idx.c2, design, contrast=RLR.GSKvDMSO, inter.gene.cor = 0.01)
cam.c2.WT.GSKvDMSO2 <- camera(b.glm, idx.c2, design, contrast=WT.GSKvDMSO, inter.gene.cor = 0.01)

options(digits=2)
head(cam.c2.GSK.WTvRLR,10)
head(cam.c2.GSK.RLRvWT,10)
head(cam.c2.DMSO.WTvRLR,15)
head(cam.c2.WT.GSKvDMSO,10)
head(cam.c2.RLR.GSKvDMSO,10)
head(cam.c2.RLR.GSKvDMSO2,10)
head(cam.c2.WT.GSKvDMSO2,10)
write.table(head(cam.c2.GSK.WTvRLR,30),file="GSK.WTvRLR_gsea_c2.txt",quote=F,sep = "\t")
write.table(head(cam.c2.DMSO.WTvRLR,30),file="DMSO.WTvRLR_gsea_c2.txt",quote=F,sep = "\t")
write.table(head(cam.c2.WT.GSKvDMSO,30),file="WT.GSKvDMSO_gsea_c2.txt",quote=F,sep = "\t")
write.table(head(cam.c2.RLR.GSKvDMSO,30),file="RLR.GSKvDMSO_gsea_c2.txt",quote=F,sep = "\t")


## barcodeplot
tiff(filename="WT.GSKvDMSO_ReactomeInfluenza.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx.c2[["REACTOME_INFLUENZA_INFECTION"]],
            
            
            main="REACTOME_INFLUENZA_INFECTION",
            alpha=1,
            
)
dev.off()

## barcodeplot
tiff(filename="RLR.GSKvDMSO_ReactomeInfluenza.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx.c2[["REACTOME_INFLUENZA_INFECTION"]],
            
            
            main="REACTOME_INFLUENZA_INFECTION",
            alpha=1,
            
)
dev.off()

## barcodeplot
tiff(filename="WT.GSKvDMSO_ReactomeNonSense.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx.c2[["REACTOME_NONSENSE_MEDIATED_DECAY_NMD"]],
            
            
            main="REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
            alpha=1,
            
)
dev.off()

## barcodeplot
tiff(filename="RLR.GSKvDMSO_ReactomeNonSense.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx.c2[["REACTOME_NONSENSE_MEDIATED_DECAY_NMD"]],
            
            
            main="REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
            alpha=1,
            
)
dev.off()

tiff(filename="WT.GSKvDMSO_ReactomeDNAmeth.tiff", units="in",width=5,height=4,res=600)
barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx.c2[["REACTOME_DNA_METHYLATION"]],
            main="REACTOME_DNA_METHYLATION",
            alpha=1)
dev.off()
tiff(filename="RLR.GSKvDMSO_ReactomeDNAmeth.tiff", units="in",width=5,height=4,res=600)
barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx.c2[["REACTOME_DNA_METHYLATION"]],
            main="REACTOME_DNA_METHYLATION",
            alpha=1)
dev.off()
tiff(filename="WT.GSKvDMSO_ReactomePRC2.tiff", units="in",width=5,height=4,res=600)
barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx.c2[["REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA"]],
            main="REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA",
            alpha=1)
dev.off()
tiff(filename="RLR.GSKvDMSO_ReactomePRC2.tiff", units="in",width=5,height=4,res=600)
barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx.c2[["REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA"]],
            main="REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA",
            alpha=1)
dev.off()













## barcodeplot
tiff(filename="WT.GSKvDMSO_ReactomePKN1.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx.c2[["REACTOME_ACTIVATED_PKN1_STIMULATES_TRANSCRIPTION_OF_AR_ANDROGEN_RECEPTOR_REGULATED_GENES_KLK2_AND_KLK3"]],
            
            
            main="REACTOME_PKN1",
            alpha=1,
            
)
dev.off()

## barcodeplot
tiff(filename="RLR.GSKvDMSO_ReactomePKN1.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx.c2[["REACTOME_ACTIVATED_PKN1_STIMULATES_TRANSCRIPTION_OF_AR_ANDROGEN_RECEPTOR_REGULATED_GENES_KLK2_AND_KLK3"]],
            
            
            main="REACTOME_PKN1",
            alpha=1,
            
)
dev.off()
