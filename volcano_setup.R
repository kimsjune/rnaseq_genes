# some cleanup needed
# Remove NAs 
b.glm.fit.symbol <- b.glm.fit
rownames(b.glm.fit.symbol) <- unname(mapIds(org.Mm.eg.db, keys=rownames(b.glm.fit.symbol$counts), keytype = "ENSEMBL", column="SYMBOL"))
omitNArows <- na.omit(rownames(b.glm.fit.symbol))
# there are some genes that once converted into 
b.glm.fit.symbol.omitNArows <- b.glm.fit.symbol[omitNArows,]
dups <- c("Txnl4a", "Taf9", "Pms2", "Bfar", "Bcl2l2", "Erdr1", "4933427D14Rik", "Hspa14", "Dohh")
dup <-duplicated(b.glm.fit.symbol.omitNArows)
b.glm.fit.symbol.omitNArows.rmDup <- b.glm.fit.symbol.omitNArows[
  !rownames(b.glm.fit.symbol.omitNArows) %in% dups,]



tr.wt.gskvdmso <- glmTreat(b.glm.fit.symbol.omitNArows.rmDup, contrast=my.contrasts[,"WT.GvD"],
                           lfc=log2(1.1))
summary(decideTests(tr.wt.gskvdmso))
wt.volcano.table<-topTags(tr.wt.gskvdmso, n="Inf")$table
rownames(wt.volcano.table)<-rownames(topTags(tr.wt.gskvdmso, n="Inf")$table)



tr.rlr.gskvdmso <- glmTreat(b.glm.fit.symbol.omitNArows.rmDup, contrast=my.contrasts[,"RLR.GvD"],
                            lfc=log2(1.1))
summary(decideTests(tr.rlr.gskvdmso))
rlr.volcano.table<-topTags(tr.rlr.gskvdmso, n="Inf")$table
rownames(rlr.volcano.table)<-rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table)

# maybe I want to generate volcano plots comparing GSK343 WT / GSK343 RLR, instead of
# GSK343 / DMSO for each genotype...

b.glm.fit.symbol <- b.glm.fit

rownames(b.glm.fit.symbol) <- unname(mapIds(org.Mm.eg.db, keys=rownames(b.glm.fit.symbol$counts), keytype = "ENSEMBL", column="SYMBOL"))

tr.gsk.wtvrlr <- glmLRT(b.glm.fit.symbol, contrast=my.contrasts[,"GSK.WTvRLR"])
summary(decideTests(tr.gsk.wtvrlr))
gsk.volcano.table<-topTags(tr.gsk.wtvrlr, n="Inf")$table
rownames(gsk.volcano.table)<-rownames(topTags(tr.gsk.wtvrlr, n="Inf")$table)

