# Version info: R 4.3.3, Biobase 2.62.0, GEOquery 2.70.0, limma 3.58.1
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# Set working directory to Output directory created
setwd("path/to/Output")

# load series and platform data from GEO
gset <- getGEO("GSE64217", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "012012"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||(qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) 
}

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("N","CIN","CC"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[3],"-",groups[1],sep=""), paste(groups[2],"-",groups[1],sep=""), paste(groups[3],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

# "SEQUENCE" is removed, #GB_ACC renamed GenBank_Accession
# tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","GI","GenBank.Accession","Gene.symbol","Gene.title"))
# write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE64217", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE64217", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# ** Skip UMAP plot and Mean-Variance Trend as requested **

################################################################
# Create an output directory
dir.create("DEG/GSE64217", recursive=TRUE)

# Volcano plots for all three contrasts in a single plot panel
# Saving the volcano plots as an image (PNG file)
png("DEG/GSE64217/GSE64217.png", width = 1200, height = 400, res = 120) 

# Defining logFC and adjusted p-value thresholds
logFC_cutoff <- 2
pval_cutoff <- 0.05

# Setting up 3 plots in one row
par(mfrow = c(1, 3), mar = c(5, 4, 4, 2))

# Looping through each contrast (assumes 3 contrasts)
for (i in 1:3) {
  # Extracting logFC and p-values
  logFC <- fit2$coefficients[, i]
  pval <- fit2$p.value[, i]
  
  # Creating logical vectors for significant up/down genes
  sig_up <- logFC > logFC_cutoff & pval < pval_cutoff
  sig_down <- logFC < -logFC_cutoff & pval < pval_cutoff
  nonsig <- !(sig_up | sig_down)
  
  # Using volcano plot from limma
  volcanoplot(fit2, coef = i, 
              col = ifelse(sig_up, "red", ifelse(sig_down, "blue", "grey")),
              main = paste("GSE64217 -", colnames(fit2)[i]), 
              xlim = c(-10, 10), ylim = c(0, max(-log10(pval), 10)),
              xlab = "log2 Fold Change", ylab = "-log10(P-value)")
  
  # Adding horizontal line for p-value threshold
  abline(h = -log10(pval_cutoff), col = "black", lty = 2)
  
  # Adding vertical lines for logFC threshold
  abline(v = c(-logFC_cutoff, logFC_cutoff), col = "black", lty = 2)
  
  # Adding legend
  legend("topright", legend=c("UP", "DOWN", "NS"), col=c("red", "blue", "grey"), pch=20, bty="n",
         x.intersp = 0.7, y.intersp = 0.7, cex = 1, pt.cex = 2)
}

# Resetting layout
par(mfrow = c(1,1))

# Closing the device
dev.off()


# Writing all significant genes of each contrast into a .csv file
for (i in 1:3) {
  # Getting topTable for each contrast
  tT2 <- topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=Inf)
  
  # Filtering significant genes
  sig_genes <- subset(tT2, abs(logFC) >= 2 & adj.P.Val < 0.05)

  # Getting upregulated and downregulated gene symbols
  upregulated <- subset(sig_genes, logFC >= 2)
  downregulated <- subset(sig_genes, logFC <= -2)
  
  # Extracting unique gene symbols for upregulated and downregulated genes
  unique_up_genes <- unique(upregulated$Gene.symbol)
  unique_down_genes <- unique(downregulated$Gene.symbol)
  
  # Writing unique gene symbols to .txt files
  writeLines(unique_up_genes, paste0("DEG/GSE64217/GSE64217_", cts[i], "_upregulated_gene_symbols.txt"))
  writeLines(unique_down_genes, paste0("DEG/GSE64217/GSE64217_", cts[i], "_downregulated_gene_symbols.txt"))
}
