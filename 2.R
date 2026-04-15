library(WGCNA)
library(pheatmap)
library(ggplot2)
library(factoextra)
library(limma)
        
datExpr0 = read.csv("exp.csv",row.names = 1)
sample_cor <- cor(datExpr0, method = "pearson") 

pdf("sample_correlation_heatmap.pdf", width = 10, height = 8)
pheatmap(
  sample_cor,
  annotation_legend = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 10,
  treeheight_col = 10,
  main = "Sample Correlation Heatmap (Before QC)"
)
dev.off()


sample_mean_cor <- rowMeans(sample_cor)
outlier_samples <- names(sample_mean_cor[sample_mean_cor < 0.6])
expr_t <- t(datExpr0)
pca_result <- prcomp(expr_t, scale. = FALSE)  

pca_df <- as.data.frame(pca_result$x[, 1:2])
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Sample <- rownames(pca_df)

eig_val <- pca_result$sdev^2
explained_var <- eig_val / sum(eig_val) * 100  

pdf("sample_PCA_before_QC.pdf", width = 8, height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 3) +
  geom_text(size = 2, hjust = 0.5, vjust = -1) +
  labs(
    x = paste0("PC1 (", round(explained_var[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "%)"),
    title = "Sample PCA (Before QC)"
  ) +
  theme_bw()
dev.off()
min_samples <- 0.5 * ncol(datExpr0) 
gene_expr_pass <- rowSums(datExpr0 > 1)  
datExpr1 <- datExpr0[gene_expr_pass >= min_samples, , drop = FALSE]

gene_mad <- apply(datExpr1, 1, mad, na.rm = TRUE)
keep_mad <- gene_mad >= quantile(gene_mad, 0.3, na.rm = TRUE)
datExpr2 <- datExpr1[keep_mad, , drop = FALSE]
datExpr2 = as.data.frame(t(datExpr2))
datExpr2[1:4,1:4]
dim(datExpr2)

gsg = goodSamplesGenes(datExpr2,verbose = 3)
print(gsg$allOK) 
sampleTree = hclust(dist(datExpr2),method = "average")


pdf(file = "1.pdf",width = 10,height = 8)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree,
     main = "Sample clustering to detect outliers",
     sub="",xlab="",cex.lab=1.5,
     cex.axis=1.5,cex.main=2)
dev.off()


datExpr=datExpr2
datTraits = read.csv("group.csv",row.names = 1)


sampleTree2 = hclust(dist(datExpr),method = "average")
expr_matrix <- as.matrix(datExpr)
mode(expr_matrix) <- "numeric"
traitColors = numbers2colors(expr_matrix,signed = FALSE)
pdf(file = "2.pdf",width = 10,height = 8)
plotDendroAndColors(sampleTree2,traitColors,
                    groupLabels = names(datTraits),
                    main =  "Sample dendrogram and trait heatmap")
dev.off()



powers = c(1:10,seq(from = 12,to=30,by=2))
sft = pickSoftThreshold(datExpr,powerVector = powers,verbose = 5)
sft$powerEstimat
cex1 = 0.9 
pdf(file = "3.pdf",width = 10,height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", 
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h=cex1,col="red")
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab = "Soft Threshold(power)",
     ylab = "Mean Connectivity",type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,
     cex=cex1,col="red")
dev.off()


power = sft$powerEstimate
net = blockwiseModules(datExpr,power=power,
                       TOMType = "unsigned",
                       minModuleSize = 50, 
                       reassignThreshold = 0,
                       mergeCutHeight = 0.6,
                       deepSplit = 1,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "testTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)
pdf(file = "4.pdf",width = 10,height = 5)
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,hang = 0.03,
                    addGuide = TRUE,guideHang = 0.05)
dev.off()

modulelabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
gm =data.frame(net$colors)
gm$color = moduleColors
head(gm)

genes = split(rownames(gm),gm$color)
save(genes,file = "genes.Rdata")

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs,datTraits,use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor,nSamples)

pdf(file = "5.pdf",width = 8,height = 8)

textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                   signif(moduleTraitPvalue,1),")",sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6,8.5,3,3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

text <-unique(moduleColors)
for (i in 1:length(text)){
  y=t(assign(paste(text[i],"expr",seq = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",seq = ".",quote = F))
}


modNames = substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(datExpr,MEs,use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
names(geneModuleMembership) = paste("MM",modNames,seq = "")
names(MMPvalue) = paste("p.MM",modNames,seq = "")

traitData = datTraits
i=2 
module = "turquoise"
assign(colnames(traitData)[i],traitData[i])
instrait = eval(parse(text = colnames(traitData)[i]))
geneTraitSignificance = as.data.frame(cor(datExpr,instrait,use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) = paste("p.GS",names(instrait),seq = "")
pdf(file = paste0("6.1.pdf"),width = 5,height =5)
column = match(module,modNames)
moduleGenes = moduleColors==module
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Moldule Membership in",module,"module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col = module)
dev.off()

module = "brown"
assign(colnames(traitData)[i],traitData[i])
instrait = eval(parse(text = colnames(traitData)[i]))
geneTraitSignificance = as.data.frame(cor(datExpr,instrait,use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) = paste("p.GS",names(instrait),seq = "")
pdf(file = paste0("6.2.pdf"),width = 5,height =5)
column = match(module,modNames)
moduleGenes = moduleColors==module
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Moldule Membership in",module,"module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col = module)
dev.off()


nSelect = 400
set.seed(10)
dissTOM = 1-TOMsimilarityFromExpr(datExpr,power = 6)

select = sample(nGenes,size = nSelect)
selectTOM = dissTOM[select,select]

selectTree = hclust(as.dist(selectTOM),method = "average")
selectColors = moduleColors[select]
install.packages("gplots")
library(gplots)
myheatcol = colorpanel(250,"red","orange","lemonchiffon")
pdf(file = "7.pdf",width = 15,height = 8)
plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss,selectTree,selectColors,col=myheatcol,main = "Network heatmap plot,selected genes")
dev.off()


MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
MET = orderMEs(cbind(MEs,instrait))

pdf(file = "8.pdf",width = 12,height = 8)
par(cex = 0.9)
plotEigengeneNetworks(MET,"",marDendro = c(0,4,1,2),marHeatmap = c(4,4,1,2),cex.lab = 0.8,xlabelsAngle 
                      = 90)
dev.off()

pdf(file = "9.pdf",width = 10,height = 8)
par(cex=1.0)
plotEigengeneNetworks(MET,"Eigengene dendrogram",marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

pdf(file = "10.pdf",width = 10,height = 8)
par(cex=1.0)
plotEigengeneNetworks(MET,"Eigengene adjacency heatmap",marDendro = c(4,5,2,2),
                      plotDendrograms = FALSE,xLabelsAngle = 90)
dev.off()

library(readxl)
WGCNA <- read_excel("WGCNA.xlsx")
mart_data <- read.table(
  "mart_export.txt", 
  sep = "\t", 
  header = TRUE, 
  stringsAsFactors = FALSE
)
str(mart_data)

filtered_data <- data.frame(
  Gene.name = mart_data[["Gene.name"]],
  Human.gene.name = mart_data[["Human.gene.name"]],
  Human.homology.type = mart_data[["Human.homology.type"]],
  Human.orthology.confidence = mart_data[["Human.orthology.confidence..0.low..1.high."]]
)
filtered_data <- unique(filtered_data)
WGCNA <- merge(
  WGCNA, 
  filtered_data, 
  by.x = "symbol", 
  by.y = "Gene.name", 
  all.x = TRUE
)

WGCNA$Human.gene.name[is.na(WGCNA$Human.gene.name)] <- "未找到同源基因"
gene_counts <- table(WGCNA$symbol)
multi_mapped <- names(gene_counts[gene_counts > 1])

if ("Human.orthology.confidence" %in% colnames(WGCNA)) {
  WGCNA <- WGCNA[order(WGCNA$symbol, -WGCNA$Human.orthology.confidence), ]  
  WGCNA <- WGCNA[!duplicated(WGCNA$symbol), ]  
}

non_standard <- grep("[^A-Za-z0-9_-]", WGCNA$Human.gene.name, value = TRUE)
WGCNA$Human.gene.name <- gsub("[^A-Za-z0-9_-]", "_",WGCNA$Human.gene.name)
write.table(
  WGCNA,
  file = "WGCNA.txt",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)
