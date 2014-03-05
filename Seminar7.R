source("http://bioconductor.org/biocLite.R")
biocLite("ShortRead")
biocLite("Rsamtools")
biocLite("easyRNASeq")
biocLite("BSgenome.Dmelanogaster.UCSC.dm3")
biocLite("biomaRt")
biocLite("edgeR")

dat <- read.table("~/R/STAT540/Seminar7/bottomly_count_table.tsv", header = TRUE, row.names = 1)
des <- read.table("~/R/STAT540/Seminar7/bottomly_phenodata.tsv", header = TRUE, row.names = 1)
str(dat)
class(dat)
all(rownames(des) == colnames(dat))


##Remove genes where all samples = 0
allzeroes <- dat[rowSums(dat)=="0", ]
str(allzeroes)
rmDat <- dat[rowSums(dat)>=1, ]
str(rmDat)
nrow(dat)
(nrow(allzeroes) + nrow(rmDat))

##Remove samples genes where all of one genotype are equal to zero
rmDat2  <- apply(rmDat, 1, function(row) all(row[1:nrow(subset(des, strain == "C57BL/6J"))] != 0) | all(row[nrow(subset(des, strain == "C57BL/6J")) + 1:nrow(subset(des, strain == "DBA/2J"))]) != 0)
summary(rmDat2)
rmDat <- rmDat[rmDat2, ]

with(des, table(strain))
group <- factor(c(rep("1", 10), rep("2", 11)))
group

dge.glm <- DGEList(counts = rmDat, group = group)
str(dge.glm)
names(dge.glm)
dge.glm[["samples"]]
nrow(dge.glm[[1]])

design <- model.matrix(~group)
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp)
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)
plotBCV(dge.glm.tag.disp)

fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))

lrt <- glmLRT(fit, coef=2)
topTags(lrt)

tt.glm <- topTags(lrt, n=Inf)
class(tt.glm)

nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])

interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]

sigEdgeR <- tt.glm$table[tt.glm$table$PValue<1e-05, ]
str(sigEdgeR)

summary(de.glm <- decideTestsDGE(lrt, p=0.05, adjust="BH"))

tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags=tags.glm)
abline(h=c(-2,2), col="blue")



DESeq

source("http://www.bioconductor.org/biocLite.R")
biocLite("DESeq")

deSeqDat <- newCountDataSet(rmDat, group)
head(counts(deSeqDat))

deSeqDat <- estimateSizeFactors(deSeqDat)
sizeFactors <- deSeqDat

deSeqDat <- estimateDispersions(deSeqDat)
plotDispEsts(deSeqDat)

results <- nbinomTest(deSeqDat, levels(group)[1], levels(group)[2])
str(results)

plotMA(results)

sigDESeq <- results[results$pval<1e-05, ]
rownames(sigDESeq) <- (sigDESeq[, 1])
sigDESeq <- sigDESeq[, 2:(ncol(sigDESeq))]


##Voom & Limma
library(limma)
norm.factor <- calcNormFactors(rmDat)
dat.voomed <- voom(rmDat, design, plot=TRUE, lib.size=colSums(rmDat) * norm.factor)
dat.voomed

fit <- lmFit(dat.voomed, design)
fit <- eBayes(fit)
topTable(fit)

summary(fit)
class(fit)
vfit <- data.frame(fit)

sigVoom <- vfit[vfit$p.value.group2<1e-05, ]
str(sigVoom)


sigFig <- list(rownames(sigEdgeR), rownames(sigDESeq), rownames(sigVoom))
str(sigFig)

sigFig2 <- model.matrix(~colnames, sigFig)

a <- vennCounts(sigFig)

vennDiagram(sigFig)
