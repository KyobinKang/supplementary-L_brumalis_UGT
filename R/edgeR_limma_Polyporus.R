#EdgeR_limma_combo DE analysis

# Set environment
library(limma)
library(edgeR)
setwd("/home/wonyong/Desktop/RNA_seq_Polyporus/")

# DGEList-object construction
## Load count data
files <- c("htseq_count_C1R1re.tab", "htseq_count_C1R2re.tab", "htseq_count_C1R3re.tab", "htseq_count_T1R1re.tab", "htseq_count_T1R2re.tab", "htseq_count_T1R3re.tab")

Fg <- readDGE(files, columns=c(1,3))
class(Fg)
dim(Fg)
## Load sample information (e.g. phenotype, genotype, batch information)
samplenames <- substring(colnames(Fg), 1, nchar(colnames(Fg)))
samplenames
colnames(Fg) <- samplenames

group <- as.factor(rep(c("control", "treatment"), c(3,3)))
Fg$samples$group <- group
lane <- as.factor(rep(c("rep1","rep2","rep3"), times = 2))

Fg$samples$lane <- lane
Fg$samples

## Load annotation
genes <- read.delim("gene_length_polyporus.txt", header = TRUE, sep = "\t")
head(genes)
Fg$genes <- genes

# Data pre-processing
cpm <- cpm(Fg)
lcpm <- cpm(Fg, log=TRUE)
table(rowSums(Fg$counts==0)==6)
keep.exprs <- rowSums(cpm>1)>=3
Fg <- Fg[keep.exprs,, keep.lib.sizes=FALSE]
dim(Fg)
Fg <- calcNormFactors(Fg, method = "TMM")
Fg$samples$norm.factors


# Differential expression analysis
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design
# Contrasting using a limma function
contr.matrix <- makeContrasts(controlvstreatment = treatment-control, levels = colnames(design))
contr.matrix
FgV <- voom(Fg, design, plot=TRUE)
FgV
Fg_vfit <- lmFit(FgV, design)
Fg_vfit <- contrasts.fit(Fg_vfit, contrasts=contr.matrix)
Fg_efit <- eBayes(Fg_vfit)
plotSA(Fg_efit)
summary(decideTests(Fg_efit))
# Examining the number of DE genes with an adjusted p-value < 0.05 showing greater than 2-fold difference (= log-FC 1)
Fg_tfit <- treat(Fg_vfit, lfc=1)
Fg_dt <- decideTests(Fg_tfit)
summary(Fg_dt)

control.vs.treatment <- topTreat(Fg_tfit, coef=1, n=Inf)
head(control.vs.treatment)
write.table(control.vs.treatment, file="DE_results_Polyporus.txt", sep = "\t")
write.table(cpm, file="cpm_Polyporus.txt", sep = "\t")

# RPKM calculation
rpkm <- rpkm(Fg, Fg$genes$gene_length)
write.table(rpkm, file="rpkm_polyporus.tab", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE) 

# Useful graphical representations of differential expression results
## plotMD
col4plotMD <- c('blue3','red4')
palette(col4plotMD)
plotMD(Fg_tfit, column=1, status=Fg_dt[,1], main="control.vs.treatment", xlim=c(-3,13), xlab=expression('log'[2]*'(CPM+0.25)'), ylab=expression('log'[2]*'(fold-change)'), ylim=c(-12,12), legend=FALSE, las=1, cex=0.7, bg.col="darkgrey")

## expression distribution, boxplot option 'range =2' to extend whiskers
palette1 <- c('blue3','red4')
palette(palette1)
boxplot(lcpm, las=2, col=group, main="", ylim=c(-8,16), outline=FALSE)
title(main="control.vs.treatment", ylab=expression('log'[2]*'(CPM+0.25)'))

## Unsupervised clustering of samples
group <- as.factor(rep(c('blue3','red4'), c(3,3)))
col.group <- group
col.group <- as.character(col.group)
plotMDS(lcpm, col=col.group, pch=1, xlab="PC 1", ylab="PC 2", cex=1.5)
title(main="control.vs.treatment")
￣

