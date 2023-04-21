# requires variables from B cells combined.R

Mm.c2total <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c2.cp.v7.1.entrez.rds"))
# from B cells combined.R, b.glm.fit.rownames already has ENTREZID swapped for ENSEMBL
idx.c2total <- ids2indices(Mm.c2total,id=rownames(b.glm.fit))

cam.c2total.GSK.WTvRLR <- camera(b.glm, idx.c2total, design, contrast=GSK.WTvRLR)
cam.c2total.GSK.RLRvWT<- camera(b.glm,idx.c2total,design,contrast=GSK.RLRvWT)
cam.c2total.DMSO.WTvRLR<- camera(b.glm,idx.c2total,design, contrast=DMSO.WTvRLR)
cam.c2total.WT.GSKvDMSO <- camera(b.glm,idx.c2total,design, contrast=WT.GSKvDMSO)
cam.c2total.RLR.GSKvDMSO <- camera(b.glm,idx.c2total,design, contrast=RLR.GSKvDMSO)
cam.c2total.RLR.GSKvDMSO2 <- camera(b.glm, idx.c2total, design, contrast=RLR.GSKvDMSO, inter.gene.cor = 0.01)
cam.c2total.WT.GSKvDMSO2 <- camera(b.glm, idx.c2total, design, contrast=WT.GSKvDMSO, inter.gene.cor = 0.01)

options(digits=2)
head(cam.c2total.GSK.WTvRLR,10)
head(cam.c2total.GSK.RLRvWT,10)
head(cam.c2total.DMSO.WTvRLR,15)
head(cam.c2total.WT.GSKvDMSO,20)
head(cam.c2total.RLR.GSKvDMSO,20)
head(cam.c2total.RLR.GSKvDMSO2,10)
head(cam.c2total.WT.GSKvDMSO2,10)
write.table(head(cam.c2total.GSK.WTvRLR,30),file="GSK.WTvRLR_gsea_c2total.txt",quote=F,sep = "\t")
write.table(head(cam.c2total.DMSO.WTvRLR,30),file="DMSO.WTvRLR_gsea_c2total.txt",quote=F,sep = "\t")
write.table(head(cam.c2total.WT.GSKvDMSO,30),file="WT.GSKvDMSO_gsea_c2total.txt",quote=F,sep = "\t")
write.table(head(cam.c2total.RLR.GSKvDMSO,30),file="RLR.GSKvDMSO_gsea_c2total.txt",quote=F,sep = "\t")


## barcodeplot
tiff(filename="WT.GSKvDMSO_ReactomeTotalInfluenza.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx.c2total[["REACTOME_INFLUENZA_INFECTION"]],
            
            
            main="REACTOME_INFLUENZA_INFECTION",
            alpha=1,
            
)
dev.off()

## barcodeplot
tiff(filename="RLR.GSKvDMSO_ReactomeTotalInfluenza.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx.c2total[["REACTOME_INFLUENZA_INFECTION"]],
            
            
            main="REACTOME_INFLUENZA_INFECTION",
            alpha=1,
            
)
dev.off()

## barcodeplot
tiff(filename="WT.GSKvDMSO_ReactomeTotalNonSense.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx.c2total[["REACTOME_NONSENSE_MEDIATED_DECAY_NMD"]],
            
            
            main="REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
            alpha=1,
            
)
dev.off()

## barcodeplot
tiff(filename="RLR.GSKvDMSO_ReactomeTotalNonSense.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx.c2total[["REACTOME_NONSENSE_MEDIATED_DECAY_NMD"]],
            
            
            main="REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
            alpha=1,
            
)
dev.off()

tiff(filename="WT.GSKvDMSO_ReactomeTotalDNAmeth.tiff", units="in",width=5,height=4,res=600)
barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx.c2total[["REACTOME_DNA_METHYLATION"]],
            main="REACTOME_DNA_METHYLATION",
            alpha=1)
dev.off()
tiff(filename="RLR.GSKvDMSO_ReactomeTotalDNAmeth.tiff", units="in",width=5,height=4,res=600)
barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx.c2total[["REACTOME_DNA_METHYLATION"]],
            main="REACTOME_DNA_METHYLATION",
            alpha=1)
dev.off()
tiff(filename="WT.GSKvDMSO_ReactomeTotalPRc2total.tiff", units="in",width=5,height=4,res=600)
barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx.c2total[["REACTOME_PRc2total_METHYLATES_HISTONES_AND_DNA"]],
            main="REACTOME_PRc2total_METHYLATES_HISTONES_AND_DNA",
            alpha=1)
dev.off()
tiff(filename="RLR.GSKvDMSO_ReactomeTotalPRc2total.tiff", units="in",width=5,height=4,res=600)
barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx.c2total[["REACTOME_PRc2total_METHYLATES_HISTONES_AND_DNA"]],
            main="REACTOME_PRc2total_METHYLATES_HISTONES_AND_DNA",
            alpha=1)
dev.off()













## barcodeplot
tiff(filename="WT.GSKvDMSO_ReactomeTotalPKN1.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx.c2total[["REACTOME_ACTIVATED_PKN1_STIMULATES_TRANSCRIPTION_OF_AR_ANDROGEN_RECEPTOR_REGULATED_GENES_KLK2_AND_KLK3"]],
            
            
            main="REACTOME_PKN1",
            alpha=1,
            
)
dev.off()

## barcodeplot
tiff(filename="RLR.GSKvDMSO_ReactomeTotalPKN1.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx.c2total[["REACTOME_ACTIVATED_PKN1_STIMULATES_TRANSCRIPTION_OF_AR_ANDROGEN_RECEPTOR_REGULATED_GENES_KLK2_AND_KLK3"]],
            
            
            main="REACTOME_PKN1",
            alpha=1,
            
)
dev.off()
