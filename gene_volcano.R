res <-(topTags(b.wt.gskvdmso.lrt,n="inf", sort.by="PValue"))
#labs.wt <- topTags(labs.wt, sort.by="PValue",n=10)
source("C:/users/jk/Documents/Lab/EnhancedVolcano.R")
library(ggplot2)
svg("WT_volcano.svg", family='sans', width=7, height=5)

EnhancedVolcano(data.frame(res)
                ,
                lab=rownames(res),
             selectLab= c("Ahr","Tagln2","Ahrr",
                           "Rgs10","Flrt3","Scin","Prkar1b",
                           "Cebpa","Nme4"),
                x='logFC',
                y='FDR',
               # xlim=c(-2.5,2.5),
              # ylim=c(0,5),
             
                pCutoff = 0.05,
                FCcutoff=0,
                gridlines.major=F,
                gridlines.minor=F,

                legendPosition='none',
                title=NULL,
                subtitle=NULL,
                labSize=9,
axisLabSize=24,
labFace='bold',
drawConnectors=T,
widthConnectors=0.5,
col=c('grey','grey','grey','red'),
colAlpha=1,
shadeAlpha=1,
xlab =" ",
ylab= " ",
caption=NULL,
boxedLabels=F)
dev.off()                
  




r.res <-(topTags(b.rlr.gskvdmso.lrt,n="inf", sort.by="PValue"))
#labs.wt <- topTags(labs.wt, sort.by="PValue",n=10)
source("C:/users/jk/Documents/Lab/EnhancedVolcano.R")
library(ggplot2)
svg("RLR_volcano.svg", family='sans', width=7, height=5)

EnhancedVolcano(data.frame(r.res)
                ,
                lab=rownames(r.res),
                selectLab= c("Ahr","Tagln2","Ahrr",
                             "Rgs10","Flrt3","Scin","Prkar1b",
                             "Cebpa","Nme4"),
                x='logFC',
                y='FDR',
                # xlim=c(-2.5,2.5),
                 ylim=c(0,15),
                
                pCutoff = 0.05,
                FCcutoff=0,
                gridlines.major=F,
                gridlines.minor=F,
                
                legendPosition='none',
                title=NULL,
                subtitle=NULL,
                labSize=9,
                axisLabSize=24,
                labFace='bold',
                drawConnectors=T,
                widthConnectors=0.5,
                col=c('grey','grey','grey','red'),
                colAlpha=1,
                shadeAlpha=1,
                xlab =" ",
                ylab= " ",
                caption=NULL,
                boxedLabels=F)
dev.off()               
