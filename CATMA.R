setwd()

library("limma")

meta2 <- read.csv("metadata.csv",as.is = T,stringsAsFactors = F)
head(meta2)


Groups <- meta2$Comparison
Groups = as.factor(Groups)
Groups
#
## Read data
files <- paste("RawArrays/",meta2$File,".gpr.gz",sep="") #Where the raw array files are
files
length(files)
x <- read.maimages(files,source="genepix",green.only=F) #For two channel arrays

## For annotation. Depending on the array version it could be different files (ie, CATMA v2.3 or CATMA v5)
catma2id <- read.table("../CATMA_2.2_07122011.txt",header = T,sep="\t",as.is = T,row.names = 1)
Genes <- catma2id[catma2id$GENE_TYPE != "",]
x$genes$ID <- toupper(x$genes$ID)
Genes$CATID <- rownames(Genes)
head(Genes)

shared <- x$genes$ID %in% rownames(Genes)
table(shared)

x$genes[,c("PROBE_TYPE","GENE_ID","GENE_TYPE","DESCRIPTION","CATID")] <- NA
x$genes[shared,] <- Genes[x$genes[shared,"ID"],]
head(x$genes)
head(x$genes[shared,])
#
plotDensities(x,main='Raw')

boxplot(data.frame(log2(x$Gb)),main="Green background")
boxplot(data.frame(log2(x$Rb)),main="Red background")
#
plotMD(x)

## Background correction and normalization
y <- backgroundCorrect(x,method="normexp")
plotMD(y)
plotDensities(y)
y <- normalizeBetweenArrays(y,method="quantile")
plotDensities(y)

## -- Filter out control probes
## To get an idea of how bright expression probes should be, we compute the 95% percentile of the negative control probes on each array. We keep probes that are at least 10% brighter than the negative controls on at least four arrays (because there are four replicates):

# --
names(y) #E for expression
head(y$genes) #
table(y$genes$GENE_TYPE)
head(y$genes[is.na(y$genes$GENE_TYPE),]) #NegativeControl
head(y$genes[!is.na(y$genes$GENE_TYPE),]) #Actual Genes

# More on: http://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/cyto-chip-quality-measures-tech-note-1570-2014-038.pdf
## "Oligo slides have QC measures incorporated during microarray manufacture and are supposed to have dark corners with a few bright spots in the extreme corners."
# --
neg95 <- apply(y$A[is.na(y$genes$GENE_TYPE),],2,function(x) quantile(x,p=0.95))
cutoff <- matrix(1.1*neg95,nrow(y),ncol(y),byrow=TRUE)
isexpr <- rowSums(y$A > cutoff) >= 2
table(isexpr)

# Keep only regular probes 
y0 <- y[!is.na(y$genes$GENE_TYPE) & isexpr,]
dim(y0)
show(y0)

## Average of repeated probes
yave <- avereps(y0,ID=y0$genes[,"GENE_ID"])


#
design = modelMatrix(meta2[,c("Cy3","Cy5")],ref="Col8R0.2mM")
design <- design[,-c(1,2),drop=F]
meta2
design
#colnames(design) <- paste(colnames(design),"_logFC",sep="")
#design <- design[,grep("nlp71",colnames(design))]
head(design)

#design <- design[,-3]
#
fit <- lmFit(yave, design)
fit2 <- eBayes(fit)


results = decideTests(fit2)
vennDiagram(results,cex=0.7) 
#
DE <- topTable(fit2, adjust="BH", number=Inf)

#
colnames(DE)[grep("logFC",colnames(DE))] <- paste("logFC",colnames(design),sep="_")
#
## Get Significant DE genes
DE.Sign <- DE[DE$adj.P.Val < 0.05,]
c(nrow(DE.Sign),nrow(DE))
#
head(DE.Sign)
dim(DE.Sign)

#View(DE.Sign)
write.csv(x = DE, "DiffExpressionAll.csv",row.names = T)
write.csv(x = DE.Sign, "Significant_DE.csv",row.names = T)

