library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

# test dataset
set.seed(123)
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))

Heatmap(mat)

# column split is manual
column_split = c(rep("group1",6), rep("group2",6))
Heatmap(t(scale(t(cpm(b))))[gene.set,], 
        # scale is applied by column by default, so it has to be transposed and untransposed
        col=colfunc(10),
        cluster_rows=F,
        cluster_columns=F,
        border_gp=gpar(col="black"),
        column_split=column_split,
        column_gap = unit(5, "mm"))

library(circlize)
col_fun <- colorRamp2(c(-2,0,3), c("blue","black","red"))
lgd = Legend(col_fun=col_fun, title="", border="white",
             at=c(-2,-1,0,1,2,3),
             direction="horizontal",
             grid_height = unit(1,"cm"),
             grid_width = unit(1,"csm"),
             labels_gp=gpar(col="black",fontsize=14))
draw(lgd)
