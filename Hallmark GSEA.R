# GSEA Hallmark gene set analysis
# a list of gene sets, each a vector of NCBI GeneIDs
Mm.hallmark <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.h.all.v7.1.entrez.rds"))

# create a table of vectors that contains row numbers that matches a gene in the gene set
idx <- ids2indices(Mm.hallmark,id=row.names(b.glm.fit.rownames))
# actual enrichment analysis
cam.GSK.WTvRLR <- camera(b.glm, idx, design, contrast=GSK.WTvRLR)
cam.GSK.RLRvWT<- camera(b.glm,idx,design,contrast=GSK.RLRvWT)
cam.DMSO.WTvRLR<- camera(b.glm,idx,design, contrast=DMSO.WTvRLR)
cam.WT.GSKvDMSO <- camera(b.glm,idx,design, contrast=WT.GSKvDMSO)
cam.RLR.GSKvDMSO <- camera(b.glm,idx,design, contrast=RLR.GSKvDMSO)
write.table(head(cam.GSK.WTvRLR,10),file="GSK.WTvRLR_gsea_hallmark_top10.txt",quote=F,sep = "\t")
write.table(head(cam.DMSO.WTvRLR,10),file="DMSO.WTvRLR_gsea_hallmark_top10.txt",quote=F,sep = "\t")
write.table(head(cam.WT.GSKvDMSO,10),file="WT.GSKvDMSO_gsea_hallmark_top10.txt",quote=F,sep = "\t")
write.table(head(cam.RLR.GSKvDMSO,10),file="RLR.GSKvDMSO_gsea_hallmark_top10.txt",quote=F,sep = "\t")

# common effect of GSK treatment in each genotype
cam.WT.GSKvDMSO <-camera(b.glm, idx, design, contrast=WT.GSKvDMSO)
head(cam.WT.GSKvDMSO,35)
write.table(head(cam.WT.GSKvDMSO,50), file="cam.WT.GSKvDMSO.hallmark.txt", quote=F,sep="\t")
cam.RLR.GSKvDMSO<- camera(b.glm,idx,design, contrast=RLR.GSKvDMSO, inter.gene.cor=0.01)
head(cam.RLR.GSKvDMSO, 35)
write.table(head(cam.RLR.GSKvDMSO,50),file="cam.RLR.GSKvDMSO.hallmark.txt",quote=F,sep="\t")

## barcodeplot
tiff(filename="GSK_WTvRLR_ALLOGRAFT.tiff",units="in",width=5,height=4,res=600)
barcodeplot(b.gsk.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_ALLOGRAFT_REJECTION"]],
            
            
            main="HALLMARK_ALLOGRAFT_REJECTION",
            alpha=1,
            
)
dev.off()
tiff(filename="GSK_WTvRLR_IFNG2.tiff",units='in',width=5,height=4,res=600)
barcodeplot(b.gsk.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
            
            
            main="HALLMARK_INTERFERON_GAMMA_RESPONSE",
            alpha=1
)
dev.off()
barcodeplot(b.gsk.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_SPERMATOGENESIS"]],
            
            
            main="HALLMARK_SPERMATOGENESIS",
            alpha=1
)
tiff(filename="GSK_WTvRLR_IFNa2.tiff",units='in',width=5,height=4,res=600)

barcodeplot(b.gsk.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
            
            
            main="HALLMARK_INTERFERON_ALPHA_RESPONSE",
            alpha=1
)
dev.off()
tiff(filename="DMSO_WTvRLR_ALLOGRAFT.tiff",units='in',width=5,height=4,res=600)
barcodeplot(b.dmso.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_ALLOGRAFT_REJECTION"]],
            
            
            main="HALLMARK_ALLOGRAFT_REJECTION",
            alpha=1
)
dev.off()
tiff(filename="DMSO_WTvRLR_IFNG.tiff",units='in',width=5,height=4,res=600)
barcodeplot(b.dmso.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
            
            
            main="HALLMARK_INTERFERON_GAMMA",
            alpha=1
)
dev.off()
tiff(filename="DMSO_WTvRLR_IFNa2.tiff",units='in',width=5,height=4,res=600)
barcodeplot(b.dmso.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
            
            
            main="HALLMARK_INTERFERON_ALPHA_RESPONSE",
            alpha=1
)
dev.off()



barcodeplot(b.dmso.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_E2F_TARGETS"]],
            
            
            main="HALLMARK_E2F_TARGETS",
            alpha=1
)
tiff(filename="DMSO_WTvRLR_inflammatory.tiff",units='in',width=5,height=4,res=600)

barcodeplot(b.dmso.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_INFLAMMATORY_RESPONSE"]],
            
            
            main="HALLMARK_INFLAMMATORY_RESPONSE",
            alpha=1
)
dev.off()
tiff(filename="GSK_WTvRLR_inflammatory.tiff",units='in',width=5,height=4,res=600)

barcodeplot(b.gsk.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_INFLAMMATORY_RESPONSE"]],
            
            
            main="HALLMARK_INFLAMMATORY_RESPONSE",
            alpha=1
)
dev.off()
tiff(filename="DMSO_WTvRLR_ifna.tiff",units='in',width=5,height=4,res=600)

barcodeplot(b.dmso.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
            
            
            main="HALLMARK_INTERFERON_ALPHA_RESPONSE",
            alpha=1
)
dev.off()
tiff(filename="GSK_WTvRLR_ifna.tiff",units='in',width=5,height=4,res=600)

barcodeplot(b.gsk.wt.rlr.lrt$table$logFC,
            index=idx[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
            
            
            main="HALLMARK_INTERFERON_ALPHA_RESPONSE",
            alpha=1
)
dev.off()

tiff(filename="WT_GSKvDMSO_E2F.tiff",units='in',width=5, height=4,res=600)

barcodeplot(b.wt.gskvdmso.lrt$table$logFC,
            index=idx[["HALLMARK_E2F_TARGETS"]],
            main="HALLMARK_E2F_TARGETS",
            alpha=1)
dev.off()
tiff(filename="RLR_GSKvDMSO_E2F.tiff",units='in',width=5, height=4,res=600)

barcodeplot(b.rlr.gskvdmso.lrt$table$logFC,
            index=idx[["HALLMARK_E2F_TARGETS"]],
            main="HALLMARK_E2F_TARGETS",
            alpha=1)
dev.off()
