# follows from B cells combined.R
Mm.c7 <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c7.all.v7.1.entrez.rds"))

idx.c7 <- ids2indices(Mm.c7, id=row.names(b.glm.fit.rownames))

# comparing GSK treated cells to each other

cam.c7.GSK.WTvRLR <- camera(b.glm, idx.c7, design, contrast=GSK.WTvRLR)
cam.c7.GSK.RLRvWT<- camera(b.glm,idx.c7,design,contrast=GSK.RLRvWT)
cam.c7.DMSO.WTvRLR<- camera(b.glm,idx.c7,design, contrast=DMSO.WTvRLR)
cam.c7.WT.GSKvDMSO <- camera(b.glm,idx.c7,design, contrast=WT.GSKvDMSO)
cam.c7.RLR.GSKvDMSO <- camera(b.glm,idx.c7,design, contrast=RLR.GSKvDMSO)
options(digits=2)
head(cam.c7.GSK.WTvRLR,10)
head(cam.c7.GSK.RLRvWT,10)
head(cam.c7.DMSO.WTvRLR,15)
head(cam.c7.WT.GSKvDMSO,15)
head(cam.c7.RLR.GSKvDMSO,15)
write.table(head(cam.c7.GSK.WTvRLR,10),file="GSK.WTvRLR_gsea_c7.txt",quote=F,sep = "\t")
write.table(head(cam.c7.DMSO.WTvRLR,10),file="DMSO.WTvRLR_gsea_c7.txt",quote=F,sep = "\t")
write.table(head(cam.c7.WT.GSKvDMSO,10),file="WT.GSKvDMSO_gsea_c7.txt",quote=F,sep = "\t")
write.table(head(cam.c7.RLR.GSKvDMSO,10),file="RLR.GSKvDMSO_gsea_c7.txt",quote=F,sep = "\t")

