# follows from B cells combined.R
Mm.c5 <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c5.bp.v7.1.entrez.rds"))

idx.c5 <- ids2indices(Mm.c5, id=row.names(b.glm.fit.rownames))

# comparing GSK treated cells to each other
GSK.WTvRLR<- makeContrasts(GSK.WT-GSK.RLR, levels=design)
GSK.RLRvWT<- makeContrasts(GSK.RLR-GSK.WT, levels=design)
DMSO.WTvRLR<- makeContrasts(DMSO.WT-DMSO.RLR, levels=design)
WT.GSKvDMSO <- makeContrasts(GSK.WT-DMSO.WT, levels=design)
RLR.GSKvDMSO <- makeContrasts(GSK.RLR-DMSO.RLR, levels=design)
cam.c5.GSK.WTvRLR <- camera(b.glm, idx.c5, design, contrast=GSK.WTvRLR)
cam.c5.GSK.RLRvWT<- camera(b.glm,idx.c5,design,contrast=GSK.RLRvWT)
cam.c5.DMSO.WTvRLR<- camera(b.glm,idx.c5,design, contrast=DMSO.WTvRLR)
cam.c5.WT.GSKvDMSO <- camera(b.glm,idx.c5,design, contrast=WT.GSKvDMSO)
cam.c5.RLR.GSKvDMSO <- camera(b.glm,idx.c5,design, contrast=RLR.GSKvDMSO)
options(digits=2)
head(cam.c5.GSK.WTvRLR,10)
head(cam.c5.GSK.RLRvWT,10)
head(cam.c5.DMSO.WTvRLR,15)
head(cam.c5.WT.GSKvDMSO,55)
head(cam.c5.RLR.GSKvDMSO,15)
write.table(head(cam.c5.GSK.WTvRLR,20),file="GSK.WTvRLR_gsea_c5.txt",quote=F,sep = "\t")
write.table(head(cam.c5.DMSO.WTvRLR,20),file="DMSO.WTvRLR_gsea_c5.txt",quote=F,sep = "\t")
write.table(head(cam.c5.WT.GSKvDMSO,20),file="WT.GSKvDMSO_gsea_c5.txt",quote=F,sep = "\t")
write.table(head(cam.c5.RLR.GSKvDMSO,20),file="RLR.GSKvDMSO_gsea_c5.txt",quote=F,sep = "\t")

