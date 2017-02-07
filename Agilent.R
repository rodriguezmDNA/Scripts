setwd("")

library(limma)

meta <- read.csv("metadata.csv",as.is = T)
meta$File <- list.files(path = "ArrayFolder",full.names = T) #Set the folder where the raw files are.
## Read data
x <- read.maimages(meta[,"File"],source="agilent",green.only=TRUE)

boxplot(x$E)
meltESET <- melt(exprs(eset))

## Background correction and normalization
y <- backgroundCorrect(x,method="normexp")
y <- normalizeBetweenArrays(y,method="quantile")
boxplot(y$E)
## -- Filter out control probes
## To get an idea of how bright expression probes should be, we compute the 95% percentile of the negative control probes on each array. We keep probes that are at least 10% brighter than the negative controls on at least four arrays (because there are four replicates):

# --
names(y) #E for expression
head(y$genes) #
table(y$genes$ControlType)
y$genes[y$genes$ControlType==-1,] #NegativeControl
y$genes[y$genes$ControlType==0,] #Actual Genes
y$genes[y$genes$ControlType==1,] #Spike in controls and Dark Corners
# More on: http://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/cyto-chip-quality-measures-tech-note-1570-2014-038.pdf
## "Oligo slides have QC measures incorporated during microarray manufacture and are supposed to have dark corners with a few bright spots in the extreme corners."
# --
neg95 <- apply(y$E[y$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95))
cutoff <- matrix(1.1*neg95,nrow(y),ncol(y),byrow=TRUE)
isexpr <- rowSums(y$E > cutoff) >= 4
table(isexpr)

# Keep only regular probes 
y0 <- y[y$genes$ControlType==0 & isexpr,]

## Filter
y0 <- y0[grep("^AT.*",y0$genes$SystematicName),]
y0$genes$SystematicName <- gsub("\\.[0-9]*$","",y0$genes$SystematicName)
yave <- avereps(y0,ID=y0$genes[,"SystematicName"])




# Now we can find genes differentially expressed for the corn oil treatments compared to the saline control
Treatment <- meta[,"Sample"]
#Treatment <- gsub("_","",Treatment)
Treatment <- factor(Treatment,levels=unique(Treatment))
design <- model.matrix(~0+Treatment) #design <- model.matrix(~Treatment)

design
colnames(design) <- levels(Treatment)
contMatrix = makeContrasts( 
  "mutVSwt"=groupB-groupCtrl, #Modify accordingly
  levels=design
)
contMatrix # Compare mutants vs wt and mutants amongst them.


fit  <- lmFit(yave, design=design)
fit = contrasts.fit(fit,contMatrix)
fit = eBayes(fit)

results = decideTests(fit)
vennDiagram(results,cex=0.7) 
#ANOVA
DE <- topTable(fit, number=nrow(fit),adjust="BH")
head(DE)

table(DE$adj.P.Val < 0.05)

DE <- DE[,c(c("GeneName","SystematicName","AveExpr","F","P.Value","adj.P.Val"),colnames(contMatrix))]

colnames(DE)[which(colnames(DE) %in% colnames(contMatrix))] <-paste("logFC",intersect(colnames(DE) , colnames(contMatrix)),sep = "")
head(DE)
load("GeneNamesSimplified.RData")
GeneNamesSimplified[intersect(GeneNamesSimplified,gsub("\\.[0-9]*$","",DE$SystematicName)),,drop=F]


rownames(DE) <- gsub("\\.[0-9]*$","",rownames(DE))

idx <- intersect(rownames(DE), rownames(GeneNamesSimplified))
DE[idx,"Symbol"] <- GeneNamesSimplified[idx,]
head(DE)
## save

normExprs <- yave$E
colnames(normExprs) <- meta$Name
colnames(normExprs) <- paste("NormExprs_",colnames(normExprs),sep="")

DE.Sign <- DE[DE$adj.P.Val < 0.05,]
dim(DE)
dim(DE.Sign)


write.csv(DE,"DifferentialExpression.csv",row.names = T)
write.csv(DE.Sign,"Significant_DifferentialExpression.csv",row.names = T)




