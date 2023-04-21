# Needs B cells combined
# ENSEMBL rownames converted into ENTREZ ID
# Then NA rows (no ENTREZ ID) have been removed

c.counts.remove.KO.genes <- b.counts.remove.KO.genes
c<-DGEList(counts=c.counts.remove.KO.genes, group=group)
keep.c<- filterByExpr(c)
c <- c[keep.c, , keep.lib.sizes=F]
c<- calcNormFactors(c)


ens.c <- rownames(c$counts)
symbols.c <- mapIds(org.Mm.eg.db, keys=ens.c, column='ENTREZID', keytype='ENSEMBL',multiVals="filter")
symbols.c <- symbols.c[match(rownames(c$counts),names(symbols.c))]
rownames(c$counts)<-symbols.c
notNA.c <- !is.na(rownames(c$counts))
c <- c[notNA.c,]

# general linear model
c.glm<-estimateDisp(c,design)
c.glm.fit<-glmFit(c.glm,design)
c.gsk.wt.rlr.lrt <- glmLRT(c.glm.fit, contrast=my.contrasts[,"GSK.WTvRLR"])
c.dmso.wt.rlr.lrt <- glmLRT(c.glm.fit, contrast=my.contrasts[,"DMSO.WTvRLR"])
c.wt.gskvdmso.lrt <- glmLRT(c.glm.fit,contrast=my.contrasts[,"WT.GvD"])
c.rlr.gskvdmso.lrt <- glmLRT(c.glm.fit,contrast=my.contrasts[,"RLR.GvD"])


go.wt.gskvdmso <- goana(c.wt.gskvdmso.lrt, species="Mm", FDR=0.05)
go.rlr.gskvdmso <- goana(c.rlr.gskvdmso.lrt, species="Mm", FDR=0.05)
go.dmso.wtvrlr <- goana(c.dmso.wt.rlr.lrt, species="Mm", FDR=0.05)
go.gsk.wtvrlr <- goana(c.gsk.wt.rlr.lrt, species="Mm", FDR=0.05)



topGO(go.wt.gskvdmso,n=10,truncate=50,sort="up")
topGO(go.rlr.gskvdmso,n=10, truncate=50,sort="up")

gsk.sum<-topGO(go.gsk.wtvrlr,n=10)
dmso.sum <-topGO(go.dmso.wtvrlr,n=10)

write.table(wt.sum, file="topGO.gsk.wtvrlr.txt",sep='\t',quote=F)
write.table(rlr.sum, file="topGO.dmso.wtvrlr.txt",sep='\t',quote=F)
dmso.topGO<-topGO(go.dmso.wtvrlr, ontology="BP",n=10, truncate=50)
gsk.topGO<-topGO(go.gsk.wtvrlr, ontology="BP",n=10, truncate=50)


dotplot(data.frame(dmso.topGO))






topKEGG(kegg.wt.gvd, n=50, truncate=50)
go.rlr.gskvdmso <- goana(rlr.gvd.lrt, species="Mm", FDR=0.05)
kegg.rlr.gvd <- kegga(rlr.gvd.lrt, species="Mm")
topGO(go.rlr.gskvdmso, n=50, truncate=50)
topKEGG(kegg.rlr.gvd, n=50, truncate=50)

b.dmso.wt.rlr.lrt.go <- b.dmso.wt..rlr.lrt

rownames(b.dmso.wt.rlr.lrt.go) <- unname(mapIds(org.Mm.eg.db, keys=rownames(b.dmso.wt.rlr.lrt.go), keytype = "ENSEMBL", column="ENTREZID"))
