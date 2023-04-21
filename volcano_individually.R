BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Mm.eg.db)
wt.names <- as.character(rownames(wt.volcano.table))
wt.fc <- as.numeric(as.character(wt.volcano.table$logFC))
wt.pval<- as.numeric(as.character(wt.volcano.table$PValue))
b.wt.gskvdmso.volcano <- data.frame(
  wt.names,
  as.numeric(wt.fc),
  as.numeric(wt.pval))

                                          
colnames(b.wt.gskvdmso.volcano) <- c("ID","Log2FC","pvalue")

selectlabs_IL6 <- unname(mapIds(org.Mm.eg.db, keys=Mm.hallmark$HALLMARK_IL6_JAK_STAT3_SIGNALING, keytype = "ENTREZID", column="SYMBOL"))
selectlabs_PI3<-unname(mapIds(org.Mm.eg.db, keys=Mm.hallmark$HALLMARK_PI3K_AKT_MTOR_SIGNALING, 
                              keytype = "ENTREZID", column="SYMBOL"))
selectlabs_E2F <- unname(mapIds(org.Mm.eg.db, keys=Mm.hallmark$HALLMARK_E2F_TARGETS, 
                              keytype = "ENTREZID", column="SYMBOL"))
               # unname(mapIds(org.Mm.eg.db, keys=Mm.c5$GO_B_cell_activation, 
                #              keytype = "ENTREZID", column="SYMBOL")))
               #innateDB.viral.gene.names)
key.vals.wt.IL6 <- ifelse(rownames(topTags(tr.wt.gskvdmso, n="Inf")$table) 
                      %in% unname(mapIds(org.Mm.eg.db,
                      keys=Mm.hallmark$HALLMARK_IL6_JAK_STAT3_SIGNALING, 
                      keytype = "ENTREZID", column="SYMBOL")),'blue','snow4')
key.vals.wt.PI3<-ifelse(rownames(topTags(tr.wt.gskvdmso, n="Inf")$table) 
                       %in% unname(mapIds(org.Mm.eg.db,
                        keys=Mm.hallmark$HALLMARK_PI3K_AKT_MTOR_SIGNALING, 
                         keytype = "ENTREZID", column="SYMBOL")),'blue','snow4')
key.vals.wt.E2F<-ifelse(rownames(topTags(tr.wt.gskvdmso, n="Inf")$table)
                       %in% unname(mapIds(org.Mm.eg.db,
                            keys=Mm.hallmark$HALLMARK_E2F_TARGETS, 
                            keytype = "ENTREZID", column="SYMBOL")),'blue','snow4')
                #ifelse(rownames(topTags(tr.wt.gskvdmso, n="Inf")$table) 
                #%in% unname(mapIds(org.Mm.eg.db,
                #keys=Mm.c5$GO_B_cell_activation, 
                #keytype = "ENTREZID", column="SYMBOL")),'darkorchid','snow4')))
                       
#key.vals.wt[is.na(key.vals.wt)] <- 'white'
names(key.vals.wt.IL6)[key.vals.wt=='snow4'] <- "NotAnnotated"
names(key.vals.wt.IL6)[key.vals.wt=='blue'] <-"IL6"
#names(key.vals.wt)[key.vals.wt=='purple'] <- "Myc"
#names(key.vals.wt)[key.vals.wt=='red'] <- "PI3K"
#names(key.vals.wt)[key.vals.wt=='darkorchid'] <- "E2F"

tiff("wt_gskvdmso_volcano_HALLMARK.tiff",units="in",width=3,height=3,res=600)
EnhancedVolcano(b.wt.gskvdmso.volcano,
                lab=rownames(topTags(tr.wt.gskvdmso, n="Inf")$table),
                x="Log2FC",
                y="pvalue",
                selectLab =  selectlabs_IL6,
                pCutoff = 0.05,
                FCcutoff= 0,
                ylim=c(1,5),
                xlim=c(-2,2),
                axisLabSize = 12,
                title=NULL,
                subtitle=NULL,
                caption=NULL,
                #transcriptPointSize = 0.5,
                legendPosition = 0,
                pointSize = 
                  ifelse(rownames(topTags(tr.wt.gskvdmso, n="Inf")$table) %in% 
                  unname(mapIds(org.Mm.eg.db,keys=Mm.hallmark$HALLMARK_IL6_JAK_STAT3_SIGNALING , 
                  keytype = "ENTREZID", column="SYMBOL")),2,0.001),
                
                      # ifelse(rownames(topTags(tr.wt.gskvdmso, n="Inf")$table)%in% 
                      # unname(mapIds(org.Mm.eg.db,
                      # keys=Mm.hallmark$HALLMARK_PI3K_AKT_MTOR_SIGNALING,
                      # keytype = "ENTREZID", column="SYMBOL")),2,
                      #   ifelse(rownames(topTags(tr.wt.gskvdmso, n="Inf")$table)%in% 
                      #   unname(mapIds(org.Mm.eg.db,keys=Mm.hallmark$HALLMARK_E2F_TARGETS,
                      #   keytype = "ENTREZID", column="SYMBOL")),2,0.2))),
                          #ifelse(rownames(topTags(tr.wt.gskvdmso, n="Inf")$table)%in% 
                          #innateDB.viral.gene.names,2.5,0.4))),
                labSize = 0,
                colAlpha = 0.6,
                gridlines.major = F,
                gridlines.minor = F,
                cutoffLineWidth = 0.5,
                borderWidth = 0.5,
                colCustom = key.vals.wt.IL6)
dev.off()
#####INSET idea
tiff("wt_gskvdmso_volcano_inset.tiff",units="in",width=7,height=7,res=600)
EnhancedVolcano(b.wt.gskvdmso.volcano,
                lab=rownames(topTags(tr.wt.gskvdmso, n="Inf")$table),
                x="Log2FC",
                y="pvalue",
                selectLab =  selectlabs,
                pCutoff = 0.05,
                FCcutoff= 0,
                ylim=c(0,5),
                #xlim=c(0.2,0.5),
                axisLabSize = 10,
                title=NULL,
                subtitle=NULL,
                #transcriptPointSize = 0.5,
                legendPosition = 0,
                pointSize = 2.5,
                labSize = 5,
                colAlpha = 1,
                gridlines.major = F,
                gridlines.minor = F,
                cutoffLineWidth = 0.5,
                borderWidth = 0.5,
                colCustom = key.vals.wt)
dev.off()
#####

























#### RLRs
rlr.names <- as.character(rownames(rlr.volcano.table))
rlr.fc <- as.numeric(as.character(rlr.volcano.table$logFC))
rlr.pval<- as.numeric(as.character(rlr.volcano.table$PValue))
b.rlr.gskvdmso.volcano <- data.frame(
  rlr.names,
  as.numeric(rlr.fc),
  as.numeric(rlr.pval))


colnames(b.rlr.gskvdmso.volcano) <- c("ID","Log2FC",
                                     "pvalue")


key.vals.rlr <- ifelse(rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table) 
  %in% unname(mapIds(org.Mm.eg.db, keys=Mm.hallmark$HALLMARK_IL6_JAK_STAT3_SIGNALING, 
  keytype = "ENTREZID", column="SYMBOL")),'blue',
                ifelse(rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table) 
                  %in% unname(mapIds(org.Mm.eg.db,
                      keys=Mm.hallmark$HALLMARK_PI3K_AKT_MTOR_SIGNALING, 
                      keytype = "ENTREZID", column="SYMBOL")),'red',
                             ifelse(rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table) %in% 
                              unname(mapIds(org.Mm.eg.db, 
                              keys=Mm.hallmark$HALLMARK_E2F_TARGETS, 
                              keytype = "ENTREZID", column="SYMBOL")),'darkorchid','snow4')))
#ifelse(rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table) 
#       %in% unname(mapIds(org.Mm.eg.db,
#      keys=Mm.c5$GO_B_cell_activation, 
#     keytype = "ENTREZID", column="SYMBOL")),'darkorchid','snow4')))

#key.vals.rlr[is.na(key.vals.rlr)] <- 'white'
names(key.vals.rlr)[key.vals.rlr=='snow4'] <- "Not Annotated"
names(key.vals.rlr)[key.vals.rlr=='blue'] <-"IL6"
#names(key.vals.rlr)[key.vals.rlr=='purple'] <- "Myc"
names(key.vals.rlr)[key.vals.rlr=='red'] <- "INFLA"
names(key.vals.rlr)[key.vals.rlr=='darkorchid'] <- "MYC"

tiff("rlr_gskvdmso_volcano_HALLMARK.tiff",units="in",width=6,height=6,res=600)
EnhancedVolcano(b.rlr.gskvdmso.volcano,
                lab=rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table),
                x="Log2FC",
                y="pvalue",
                selectLab =  selectlabs,
                pCutoff = 0.05,
                FCcutoff= 0,
                ylim=c(1,5),
                xlim=c(-2,2),
                axisLabSize = 12,
                title=NULL,
                subtitle=NULL,
                caption=NULL,
                #transcriptPointSize = 0.5,
                legendPosition = 0,
                pointSize = 
                  ifelse(rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table) %in% 
                           unname(mapIds(org.Mm.eg.db,
                                         keys=Mm.hallmark$HALLMARK_IL6_JAK_STAT3_SIGNALING , 
                                         keytype = "ENTREZID", column="SYMBOL")),2,
                         ifelse(rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table)%in% 
                                  unname(mapIds(org.Mm.eg.db,
                                        keys=Mm.hallmark$HALLMARK_PI3K_AKT_MTOR_SIGNALING,
                                                keytype = "ENTREZID", column="SYMBOL")),2,
                                ifelse(rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table)%in% 
                                         unname(mapIds(org.Mm.eg.db,
                                              keys=Mm.hallmark$HALLMARK_E2F_TARGETS,
                                              keytype = "ENTREZID", column="SYMBOL")),2,0.2))),
                #ifelse(rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table)%in% 
                #innateDB.viral.gene.names,2.5,0.4))),
                labSize = 0,
                
                colAlpha = 0.6,
                gridlines.major = F,
                gridlines.minor = F,
                cutoffLineWidth = 0.5,
                borderWidth = 0.5,
                colCustom = key.vals.rlr)
dev.off()
