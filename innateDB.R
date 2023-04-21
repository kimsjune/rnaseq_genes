BiocManager::install("readxl")
library(readxl)
innateDB <-read_excel("innatedb_curated_genes.xls")
innateDBviral <-dplyr::filter(as.data.frame(innateDB), grepl('virus|viral|inflamm',innateDB$Annotation))
innateDB.gene.names <- unique(innateDB$`Gene Symbol`)
innateDB.viral.gene.names <- unique(innateDBviral$`Gene Symbol`)
rownames(topTags(tr.rlr.gskvdmso, n="Inf")$table) %in% innateDB.gene.names
