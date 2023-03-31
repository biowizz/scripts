#Export counts from dd.RData
write.table(counts(dds),file="counts.txt", sep="\t")



## load edgeR
library('edgeR')

### Input raw count data into edgeR ###
### DUplicate entries exists, hence merging duplicate to add the sum of all totals
count_data <- read.table("counts.txt", header=FALSE, skip=1)
colnames(count_data) <- c("genes", "NC", "NNC", "NNT", "NT")
rawdata <- aggregate(count_data[, 2:5], by=list(names=count_data$genes), FUN=sum)


dim(rawdata)

### Converting it to a mtrix
gene_names <- rawdata[,1] # Assumes gene names are in first column
rawdata_matrix <- as.matrix(rawdata[,2:ncol(rawdata)])
rownames(rawdata_matrix) <- gene_names

## Read a design file with sample information ###
design<- read.table("design.txt", sep ="\t", row.names=1, header =TRUE)

### Make class label ###
group <- factor(design$Group)

# Make DGEList object ###
y <- DGEList(counts=rawdata_matrix, group=group)


keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]


### Normalization for composition bias
y <- calcNormFactors(y)
y$samples

# Heatmap 
# Log2-transformed counts
library(pheatmap)
logcpm <- cpm(y, log=TRUE)

# Create heatmap using pheatmap package
pheatmap(logcpm, cluster_rows=TRUE, cluster_cols=TRUE, scale="row", main="Heatmap of log2-transformed counts")


############### Estimate dispersion ###
# As no replicate we take a acceptable dispersion value
# Typical values for the common BCV (square-rootdispersion) for datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates
bcv <- 0.4
### Result for NC - NNC ####

NC_NNC <- exactTest(y, pair=c("NC","NNC"),dispersion=bcv^2)
res.NC_NNC <- topTags(NC_NNC)

con.test <- topTags(NC_NNC, n=nrow(NC_NNC$table))
write.table(con.test, file="NC_NNC_result.txt", sep = "\t")


# Subset the normalized counts for the top 50 differentially expressed genes
top50 <- rownames(con.test$table)[1:50]
nc_nnc_counts <- cpm(y, log=TRUE)[top50, c("NC", "NNC")]

# Generate the heatmap
library(gplots)
heatmap.2(nc_nnc_counts, trace="none", margins=c(10, 5), dendrogram="column", col=colorRampPalette(c("green", "red"))(50))


### Result for NC - NT ####
NC_NT <- exactTest(y, pair=c("NC","NT"),dispersion=bcv^2)
res.NC_NT <- topTags(NC_NT)

con.test <- topTags(NC_NT, n=nrow(NC_NT$table))
write.table(con.test, file="NC_NT_result.txt", sep = "\t")