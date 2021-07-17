dat1 <- read.delim("exp.txt")
ID_TYPE <- "SYMBOL"
#LOGFC_CUTOFF <- 1
#P_CUTOFF <- 0.05

library("DOSE")
library("GO.db")
library("GSEABase")
library("org.Hs.eg.db")
library("clusterProfiler")
library("dplyr")
library("tidyr")
library("ggplot2")
library("stringr")
library("RColorBrewer")
library("rWikiPathways")
library("enrichplot")
library("Rgraphviz")
library("topGO")


#up.genes <- dat1[dat1$logFC > LOGFC_CUTOFF & dat1$PValue < P_CUTOFF, 1]
#dn.genes <- dat1[dat1$logFC > -LOGFC_CUTOFF & dat1$PValue < P_CUTOFF, 1]
bkgd.genes <- dat1[,1]


#### ENtrez gene list for background, up and down
bkgd.genes.entrez <- bitr(
  bkgd.genes,
  fromType = ID_TYPE,
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db)
  

#### GO Biological process (BP) term enrichment for up and down separate
#### GO ALL term (Biological process BP, Molecular function MF and cellular Component CC) enrichment for up and down separate
##############################################################################################################################
goALL <- clusterProfiler::enrichGO(
  gene     = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "ALL",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = TRUE)
 
write.table(goALL, file="goALL.txt", sep="\t")

##### Plot generation
# Plots for up regulated
pdf("GO_ALL_BAR.pdf")
barplot(goALL, showCategory = 20)
dev.off()
pdf("GO_ALL_DOT.pdf")
dotplot(goALL, showCategory = 20)
dev.off()
pdf("GO_ALL_EMAP.pdf")
emapplot(goALL, showCategory = 20)
dev.off()
##############################################################################################################################
#### GO BP term (Biological process BP, Molecular function MF and cellular Component CC) enrichment for up and down separate

goBP <- clusterProfiler::enrichGO(
  gene     = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = TRUE)

write.table(goBP, file="goBP.txt", sep="\t")
##### Plot generation
# Plots for up regulated
pdf("GO_BP_BAR.pdf")
barplot(goBP, showCategory = 20)
dev.off()
pdf("GO_BP_DOT.pdf")
dotplot(goBP, showCategory = 20)
dev.off()
pdf("GO_BP_EMAP.pdf")
emapplot(goBP, showCategory = 20)
dev.off()

##### Plot Go Grpah using TopGO
library(Rgraphviz)

pdf("BP_GO_graph.pdf")
plotGOgraph(goBP, firstSigNodes=10, useInfo="all", sigForAll=TRUE, useFullNames=TRUE)
dev.off()

#### GO MF term (Biological process BP, Molecular function MF and cellular Component CC) enrichment for up and down separate

goMF <- clusterProfiler::enrichGO(
  gene     = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "MF",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = TRUE)

write.table(goMF, file="goMF.txt", sep="\t")

##### Plot generation
# Plots for up regulated
pdf("GO_MF_BAR.pdf")
barplot(goMF, showCategory = 20)
dev.off()
pdf("GO_MF_DOT.pdf")
dotplot(goMF, showCategory = 20)
dev.off()
pdf("GO_MF_EMAP.pdf")
emapplot(goMF, showCategory = 20)
dev.off()

pdf("MF_GO_graph.pdf")
plotGOgraph(goMF, firstSigNodes=10, useInfo="all", sigForAll=TRUE, useFullNames=TRUE)
dev.off()

#### GO CC term (Biological process BP, Molecular function MF and cellular Component CC) enrichment for up and down separate

goCC <- clusterProfiler::enrichGO(
  gene     = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = TRUE)

write.table(goCC, file="goCC.txt", sep="\t")  

##### Plot generation
# Plots for up regulated
pdf("GO_CC_BAR.pdf")
barplot(goCC, showCategory = 20)
dev.off()
pdf("GO_CC_DOT.pdf")
dotplot(goCC, showCategory = 20)
dev.off()
pdf("GO_CC_EMAP.pdf")
emapplot(goCC, showCategory = 20)
dev.off()


pdf("CC_GO_graph.pdf")
plotGOgraph(goCC, firstSigNodes=10, useInfo="all", sigForAll=TRUE, useFullNames=TRUE)
dev.off()


########### Pathway KEGG EnrichKEGG
kegg_list <- enrichKEGG(gene         = bkgd.genes.entrez[[2]],
                 organism     = 'hsa',
				 pAdjustMethod = "none",
                 pvalueCutoff = 0.05)

# Plots for up D14 enrichment KEGG
pdf("KEGG_barplot.pdf")
barplot(kegg_list, showCategory = 20)
dev.off()
pdf("KEGG_dotplot.pdf")
dotplot(kegg_list, showCategory = 20)
dev.off()
pdf("KEGG_emapplot.pdf")
emapplot(kegg_list, showCategory = 20)
dev.off()

	
##########################################################################################
# Extra plots
pdf("GO_BP_heatplot.pdf")
heatplot(goBP, foldChange = geneList)
dev.off()
pdf("GO_MF_heatplot.pdf")
heatplot(goMF, foldChange = geneList)
dev.off()
pdf("GO_CC_heatplot.pdf")
heatplot(goCC, foldChange = geneList)
dev.off()

pdf("GO_BP_upsetplot.pdf")
upsetplot(goBP)
dev.off()
pdf("GO_MF_upsetplot.pdf")
upsetplot(goMF)
dev.off()
pdf("GO_CC_upsetplot.pdf")
upsetplot(goCC)
dev.off()

################### For disese categoty #######################
x <- enrichDO(gene          = bkgd.genes.entrez[[2]],
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "none",
              minGSSize     = 5,
              maxGSSize     = 500,
              readable      = FALSE)

pdf("DiseaseOntology_BAR.pdf")
barplot(x, showCategory = 20)
dev.off()
pdf("DiseaseOntology_DOT.pdf")
dotplot(x, showCategory = 20)
dev.off()
pdf("DiseaseOntology_EMAP.pdf")
emapplot(x, showCategory = 20)
dev.off()
	
###### Disease gene Network	
dgn <- enrichDGN(bkgd.genes.entrez[[2]])

pdf("DiseaseGeneNetwork_BAR.pdf")
barplot(dgn, showCategory = 20)
dev.off()
pdf("DiseaseGeneNetwork_DOT.pdf")
dotplot(dgn, showCategory = 20)
dev.off()
pdf("DiseaseGeneNetwork_EMAP.pdf")
emapplot(dgn, showCategory = 20)
dev.off()
			  