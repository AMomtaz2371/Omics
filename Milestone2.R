#### used libraries ######
# Install BiocManager
# =========================
install.packages("BiocManager")

# =========================
# CRAN Packages Installation
# =========================
install.packages("readr")
install.packages("ggfortify")
install.packages("rgl")
install.packages("plot3D")
install.packages("plotly")
install.packages("scatterplot3d")
install.packages("genefilter")
install.packages("matrixStats")
install.packages("readxl")
BiocManager::install("genefilter")
# =========================
# Load Libraries
# =========================
library(readr)
library(DESeq2)
library(ggfortify)
library(rgl)
library(plot3D)
library(plotly)
library(stats)
library(scatterplot3d)
library(genefilter)
library(matrixStats)
library(ComplexHeatmap)
library(readxl)
library(EnhancedVolcano)
##### load the mRNA-Seq data #####
setwd("C:/Users/ZeroOne/Downloads/openxlsx")
exp <- read_tsv("C:/Users/ZeroOne/Downloads/GSE268366_expressed_gene_reads.txt.gz")
gene_annotation <- exp[, c(1, 12, 13, 14, 15)]
exp <- exp[, -c(1, 12, 13, 14, 15)]
exp <- as.matrix(exp)

noversion_gene_ids <- sub("\\..*", "", gene_annotation[[1]])
rownames(exp) <- noversion_gene_ids

sum(duplicated(rownames(exp)))

exp_df <- as.data.frame(exp)
exp_df$id <- rownames(exp_df)
exp.data.agg <- aggregate(. ~ id, data = exp_df, FUN = mean)
rownames(exp.data.agg) <- exp.data.agg$id
exp.data.agg <- exp.data.agg[, -which(colnames(exp.data.agg) == "id")]

pheno <- read_excel("Metadata2.xlsx")


exp.data.agg <- exp.data.agg[rowMeans(exp.data.agg) > 1, ]


boxplot(log2(exp.data.agg + 1),
        main = "RNA-Seq Box Plot",
        col = seq(1:ncol(exp.data.agg)))
plot(density(apply(exp.data.agg, 2, mean, na.rm = TRUE)),main='density plot',cex.axis=0.5)
# transpose data to prepare it to PCA
exp_t = t(exp.data.agg)
dim(exp_t)

# exp_t <- exp_t[, apply(exp_t, 2, var, na.rm = TRUE) != 0] #remove constant genes (variance = 0)
# dim(exp_t)
# exp_t.clean <- na.omit(exp_t) #remove rows which contain all zeros
# dim(exp_t.clean)

# 2. 2D PCA
exp.pca = prcomp(exp_t , center = TRUE , scale. = TRUE)
summary(exp.pca)
autoplot(exp.pca, data = pheno, colour = 'condition')

exp = t(exp_t)


# 4. 3D PCA
View(exp.pca$x)
pca_scores = as.data.frame(exp.pca$x)
mycolors = c(rep('blue' , dim(pheno)[1]))
mycolors[which(pheno$type=='WT mice control diet')]='red'
mycolors[which(pheno$type=='WT mice Etoh diet plus binge')]='black'
mycolors[which(pheno$type=='Fkbp5 KO mice control diet')]='green'
mycolors[which(pheno$type=='Fkbp5 KO mice Etoh diet plus binge')]='brown'


plot3d(pca_scores[,1:3], pch=20, col = mycolors ,radius = '2')

# ############### impute the missing values mean ##############
# # Calculate the proportion of zeros in each row
# prop_zeros <- rowSums(exp == 0) / ncol(exp)
# 
# # Identify rows with less than 40% zeros
# rows_to_fill <- which(prop_zeros < 0.4)
# imputed_genes = rows_to_fill
# 
# # Calculate the row means for these rows
# row_means <- rowMeans(exp[rows_to_fill, ], na.rm = TRUE)
# 
# # Replace the zeros in these rows with the row means and remove others rows
# exp[rows_to_fill, ][exp[rows_to_fill, ] == 0] <- row_means
# exp = exp[imputed_genes,]
# 
exp.count = round(exp)

############### Defferential expression analysis using DESeq2 #########
table(pheno$type)

dds = DESeqDataSetFromMatrix(countData = exp.count , colData = pheno , design = ~type)
dds.run = DESeq(dds)

res=results(dds.run)
res=res[complete.cases(res),]

########### make contrasts (comparisons between all conditins) ########
contrast1 = results(dds.run, contrast = c("type", "WT mice control diet",
                                          "Fkbp5 KO mice Etoh diet plus binge"))
contrast1 = contrast1[complete.cases(contrast1),]
contrast1 = as.data.frame(contrast1)

############# get the DEGs based on adj pval , LFC ############
res.df = as.data.frame(res)
res.degs = res.df[res.df$padj<0.05 & abs(res.df$log2FoldChange)>log2(1.5),]
res.degs=res.degs[order(res.degs[,6]) ,]
degs.genes= rownames(res.degs) 
exp.degs=exp[degs.genes,]
write.table(degs.genes,file = "DEGs.txt",row.names = F,col.names = F,quote = F)

############# do normalization for all exp data to further analysis ########
ntd=normTransform(dds)
exp.norm= assay(ntd)

############### creating a heatmap for the top 100 DEG genes #####
exp.degs=exp.norm[degs.genes,]
top100_DEGS = row.names(exp.degs)[1:100]
exp100_DEGS = exp.degs[top100_DEGS,]

column_ha = HeatmapAnnotation(sample.type = pheno$type)
Heatmap(exp100_DEGS,name = 'Exp', row_names_gp= gpar(fontsize=3) , column_names_gp = gpar(fontsize=10)
        , top_annotation = column_ha)


############### 2D PCA ############3
expression_t= t(exp100_DEGS)
expression.pca = prcomp(expression_t , center = TRUE , scale. = TRUE)
summary(expression.pca)
autoplot(expression.pca, data = pheno, colour = 'type')

############### 3D PCA ##########
pca_scores = as.data.frame(expression.pca$x)
mycolors = c(rep('blue' , dim(pheno)[1]))
mycolors[which(pheno$type =='WT mice control diet')]='red'
mycolors[which(pheno$type =='WT mice Etoh diet plus binge')]='green'
plot3d(pca_scores[,1:3], pch=20, col = mycolors , type = 's',radius = '0.5')


EnhancedVolcano(res.df,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                xlim =  c(-2, 2),
                ylim = c(0, 1.5),
                pCutoff = 0.05,
                pointSize = 2,
                FCcutoff = 1,
                pCutoffCol='pvalue',
                title = "(fold change cutoff = 1, adj.pvalue cutoff = 0.05)",
)

# 1. boxplot
boxplot(exp100_DEGS,main = "processed genes Box Plot",col=seq(1:15))

# 2. histogram
cols=colnames(exp100_DEGS)
hist(exp100_DEGS,main = "processed genes Histogram")
getwd()

