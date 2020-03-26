
# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# To run this script, use
rmarkdown::render(dirname(rstudioapi::getActiveDocumentContext()$path), output_file="Report.html", knit_root_dir=wd)  

# Reading data
countsanouk <- read.table("20171025_AZAAL_countTable.txt",row.names=1,header=TRUE, sep="\t", na.strings = c("NA","#N/A")) 
dim(countsanouk)

# Show first rows; last two columns [30,31] are not part of the data (Excel leftover..) and are hence removed
head(countsanouk)

countsanouk <- countsanouk[,-(30:31)]

# Retrieve annotation columns
annot <- countsanouk[,1:5]

# Remove annotation columns and round data to obtain counts
countsanouk <- round(countsanouk[,-(1:5)])

# Read design and customize
designanouk0 <- read.table("20171025_AZAAL_design file.txt",header=TRUE, sep="\t", na.strings = c("NA","#N/A")) 
dim(designanouk0)
designanouk1 <- designanouk0[-c(1:5,30:36),]
dim(designanouk1)
rownames(designanouk1) <- designanouk1[,1]
designanouk <- designanouk1[,-1]
designanouk
colnames(designanouk) <- c("Type","Repeat")
dim(designanouk)

### Analysis

##Normalization on the ENTIRE data set (all 24 data columns and rows)
# Activate edgeR
if (!("edgeR" %in% installed.packages())) {install.packages("edgeR")}
library(edgeR)
d <- DGEList(counts = countsanouk, group = designanouk$Type)
d <- calcNormFactors(d)

# Computed normalization factors and library sizes
nfac <- d$samples$norm.factors
ls <- d$samples$lib.size
nfac
ls

# Obtaining normalized data; relative library size computed wrt geometric mean
rellibsize <- ls/exp(mean(log(ls)))

# Normalization factors (counts per sample are divided by these numbers)
nf = nfac*rellibsize
nf
normcount = round(sweep(d$counts, 2, nf, "/"))
#save(normcount,file="normcount.Rdata")

#write to .txt file; first add ID and annotation again
normcount <- data.frame(ID= rownames(normcount),annot,normcount)
#write.table(normcount,file="normdata.txt",row.names = F,sep="\t")


# ##edgeR Analysis 7: MOCK2 (control) vs FUT9

# ###Select data for the comparison (| means OR)
whdata <- which(designanouk$Type == "MOCK2" | designanouk$Type == "FUT9")
whdata

# Selecting the rows and columns
dcomp <- d[,whdata]
dcomp

# Discard tags with more than 4 zeros
nzero <- apply(dcomp$counts,1,function(rij) length(rij[rij==0]))
wh <- which(nzero > 4)
dcompf <- dcomp[-wh, ]
dim(dcompf)
ntest <- nrow(dcompf)

# ###Multi-dimensional scaling (MDS) plot
# MDS is like principle Comp Anal (PCA) but more suitable for count data
# Conclusion: shows clear separation between two groups
plotMDS(dcomp)

# ##edgeR differential expression analysis

# Estimate common dispersion (dispersion models the variability for the counts)
dcompf <- estimateCommonDisp(dcompf)

# Estimates tagwise dispersions from common and tag-specific estimates (shrinkage)
dcompf <- estimateTagwiseDisp(dcompf)

# Compute p-value for two group setting
TestMock2FUT9 <- exactTest(dcompf,pair=c("MOCK2", "FUT9"))

top<-topTags(TestMock2FUT9,n=ntest)

# ###Write results to .txt file; 
# We filtered genes so we have to do that for the annotation as well
annotfilt <- annot[-wh,]

# topTags sorts the genes w.r.t. significance. Hence, we need to match the geneIDs 
# of the annotation with thos of top (=topTags result)
ID=rownames(top)
IDannot <- rownames(annotfilt)
matchnames <- match(ID,IDannot)
annotfiltmatch <- annotfilt[matchnames,]
# check
rownames(annotfiltmatch)[1:10] == rownames(top)[1:10]

# Now write data to tab-delimited file, which can be viewed in Excel
res <- data.frame(ID,annotfiltmatch,top)
#write.table(res,file="results.txt",row.names = F,sep="\t")

## merge tables normalized counts and DEGs
normres <- merge(normcount,res,by=1,all=TRUE)
#write.table(normres,file="07_merge.txt",row.names = F,sep="\t")

# ###Extract significant genes and visualize with <0.1

# Show the genes that are given an FDR multiple testing corrected p-value less than 0.1.
# FDR = False Dicovery Rate; computed with popular Benjamini-Hochberg method
topsig <- top[top$table$FDR < 0.1,]

# Many significant effects
nrow(topsig)

# IDs sign effects
topnames <- rownames(topsig$table)

# Show top 10
topsign <- topTags(TestMock2FUT9, n=10)
topsign

# Show the data for the most significant features
topnamesn <- rownames(topsign$table)
dcompf$counts[topnamesn, ]

# ###Plotting fold changes; red dots are differentially signif

#pdf(file="myplot.pdf"); saves plot as pdf. Then also run dev.off() after plotting
plotSmear(dcompf, de.tags = topnames, main = "Fold-change MOCK2 vs FUT9",cex=0.5)
abline(h = c(-2, 2), col = "dodgerblue")
#dev.off()

# ###Extract significant genes and visualize with <0.05

# Show the genes that are given an FDR multiple testing corrected p-value less than 0.1.
# FDR = False Dicovery Rate; computed with popular Benjamini-Hochberg method
topsig2 <- top[top$table$FDR < 0.05,]

# Many significant effects
nrow(topsig2)

# IDs sign effects
topnames2 <- rownames(topsig2$table)

# Show top 10
topsign2 <- topTags(TestMock2FUT9, n=10)
topsign2

# Show the data for the most significant features
topnamesn2 <- rownames(topsign2$table)
dcompf$counts[topnamesn2, ]

# ###Plotting fold changes; red dots are differentially signif

#pdf(file="myplot2.pdf"); saves plot as pdf. Then also run dev.off() after plotting
plotSmear(dcompf, de.tags = topnames2, main = "Fold-change MOCK2 vs FUT9",cex=0.5)
abline(h = c(-2, 2), col = "dodgerblue")
#dev.off()

# ###Volcano Plot
# extract coloms and transform FDR to -log10
head(res)
dim(res)
tab = data.frame(logFC = res[,7], FDR = res[,10], negLogPval = -log(res[,10]))
head(tab)

# remove FDR, keep logFC and -log10 FDR
tab2 = data.frame(logFC = res[,7], negLogPval = -log(res[,10]))
head(tab2)

# #Generate Volcano plot:
# dimension:default is 5, 4, 4, 4
par(mar = c(5, 5, 4, 4))

# define cut-offs
lfc = 2
pval = 0.01
signGenes = (abs(tab2$logFC) > lfc & tab2$negLogPval > -log10(pval))
dim(signGenes)
head(signGenes)

# volcanoplot:
plot(tab2, pch = 16, cex = 0.6, xlab = expression (log[2]~fold~ Change), ylab = expression(-log[10]~pvalue), main = "FUT9 vs MOCK2", abline(h = -log10(pval), v = c(-lfc, lfc), col= c("green3", "blue", "blue"), lty = 2))
points(tab2[signGenes , ], pch = 16, cex = 0.8, col = "darkgreen")
legend("topright", c("P value>0.01", "logFC>2"), lty=2, col=c("green3", "blue"))
