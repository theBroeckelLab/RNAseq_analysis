##**## BASIC EXPRESSION ANALYSIS IN R ##**##



#######################
## READ IN OUR DATA ###
#######################
##**## for this example, we're using a subset of whole transcriptome RNAseq heart tissue data from GTEx (https://gtexportal.org/home/). you can run each line by putting your curser at the beginning of the line and pressing  Ctrl+Enter##**##

##read in raw counts table (you will have to change the folder from C:/Users/jblamer to wherever your data is saved. for example, if you saved these materials to your downloads folder, you'll have to use something like "C:/Users/Downloads/2023-04-27 bioinformatics materials")
counts=read.table("Z:/Projects/Project Management/Analysis/12.16.2019/Data Analysis/2023-04-27 bioinformatics materials/gtex_rawCounts.txt", header=T, row.names=1)

##read in raw meta data (you will have to change the folder from C:/Users/jblamer to wherever your data is saved. for example, if you saved these materials to your downloads folder, you'll have to use something like "C:/Users/Downloads/2023-04-27 bioinformatics materials")
metadata=read.table("Z:/Projects/Project Management/Analysis/12.16.2019/Data Analysis/2023-04-27 bioinformatics materials/gtex_metadata.txt", header=T, row.names=NULL, stringsAsFactors=F, sep="\t")

## take a look at the counts and metadata. here's the information/metadata for the first three individuals followed by their gene expression counts 
View(metadata[1:3,])
View(counts[,1:3])

## before running any analysis, make sure that the IDs between the counts and the metadata match
colnames(counts)
metadata$Sample.ID
identical(colnames(counts), metadata$Sample.ID)





#########################
## FILTER OUR DATA ######
#########################
##**## before processing our data we sometimes want to remove certain samples or subsets of samples.  if we can get the index of these samples using 'which()', then we can easily filter them from the counts and metadata ##**##

## get indexes of samples
which(metadata$Sex=="male") ##return the indexes of all males
which(metadata$Sex=="male"&metadata$Age_Bracket=="50-59") ##return the indexes of all males aged 50-59
which(metadata$Sample.ID=="GTEX.13NYB.0226") ##return the index of sample GTEX.13NYB.0226

## remove sample GTEX.13NYB.0026
index.filter=which(metadata$Sample.ID=="GTEX.13NYB.0226") ##return the index of sample GTEX.13NYB.0226 and save to a variable called 'index.filter'
metadata<- metadata[-index.filter,]
counts=counts[,-index.filter]

## now we have 99 samples, not 100
nrow(metadata)




########################################
## REMOVE LOW ABUNDANCE READS ##########
#########################################
##**## usually, a majority of genes are not expressed in our dataset.  it's best to remove low- or zero-count genes before running analysis. there are a few ways to do this but here, we're going calculate the median  count for each gene across our 99 samples and remove any genes with a median of 10 or less ##**##

## CDC42-IT1 is an example of a low-abundence gene
which(rownames(counts)=="CDC42-IT1")
counts[651,]

## calculate gene medians using the rowMedians() function from the matrixStats package
install.packages("matrixStats")
library(matrixStats)
gene.medians=rowMedians(as.matrix(counts))

## view median counts for each gene
View(data.frame(gene=rownames(counts), median.count=gene.medians))

## get the indexes of the genes with less than 10 median read count
index.lowCountGenes=which(gene.medians<=10)

## remove these low abundence genes, now 18,712 genes are remaining
counts=counts[-index.lowCountGenes,]
dim(counts)



########################################
## NORMALIZE RAW READ COUNTS ###########
#########################################
##**## when we cluster expression data (PCAs, dendrograms, heatmaps) we usually want to first normalize our raw counts.  Normalizing the data beforehand helps reduce the impact of technical factors, outlier samples, and  very highly expressed genes. there are a few different ways to normalize raw read counts, here we'll use a  variance stabilizing transformation from the DESeq package. we'll use DESeq again to run differential expression analysis later in the script ##**##

## install DESeq2
install.packages("DESeq2")
library(DESeq2)

## normalize counts with vst (similar to a log transformation of the data)
counts.normalized=varianceStabilizingTransformation(as.matrix(counts))

## take a look at the first 5 genes for the first 3 samples. counts are now normalized, not raw
counts[1:5,1:3]              ##raw counts
counts.normalized[1:5,1:3]   ##normalized counts 



###########################################
## FILTER OUT THE LEAST VARIABLE GENES ####
###########################################
##**## even after applying our median read count of 10, we have a large number of genes (18,712) remaining. often  we can reduce this number while still capturing the characteristics of the full dataset by taking only the most  variable genes (the genes that change the most across the cell lines). this is an optional step, but can reduce computational burden of processing the full dataset ##**##

## calculate gene variances
gene.vars=rowVars(as.matrix(counts.normalized))

## view expression variability for each gene
expression.variabality_dataframe=data.frame(gene=rownames(counts), expression.variabality=gene.vars)
View(expression.variabality_dataframe)

## BMP10 is the most variable gene, meaning it's expression changes dramatically sample-to-sample
index.most_var_gene=which.max(expression.variabality_dataframe$expression.variabality)
expression.variabality_dataframe[index.most_var_gene,]
plot(counts.normalized[which(rownames(counts.normalized)=="BMP10"),], xaxt="n", xlab="", ylab="BMP10 Normalized Expression", ylim=c(3,17))
axis(1, at=1:nrow(metadata), labels=metadata$Sample.ID, las=2)

## UNC50 is the least variable gene, meaning it's expression rarely changes sample-to-sample
index.least_var_gene=which.min(expression.variabality_dataframe$expression.variabality)
expression.variabality_dataframe[index.least_var_gene,]
plot(counts.normalized[which(rownames(counts.normalized)=="UNC50"),], xaxt="n", xlab="", ylab="UNC50 Normalized Expression", ylim=c(3,17))
axis(1, at=1:nrow(metadata), labels=metadata$Sample.ID, las=2)

## get the indexes of the 500 genes with the most variance
index.highVaGenes=head(order(gene.vars, decreasing=T), 500)

## keep the top 500 most variable genes and remove the others, now only 500 genes are remaining
counts.normalized=counts.normalized[index.highVaGenes,]
dim(counts.normalized)





##################
#### PCA #########
##################
##**## principal component analysis (PCA) is a common way to cluster our expression data.  the theory behind it is somewhat complicated but basically we're trying to summarize the variability across all genes/samples into a subset of principal components or PCs. We'll end up with a number of PCs (PC1, PC2, PC3, etc) and we can plot them to better understand patterns in our data ##**##

## generate our PCs
pca <- prcomp(t(counts.normalized), scale=F)

## now instead of a table of samples/genes, we have a table of samples/PCs
pca$x

## we can plot these PCs on a scatter plot. example, here's PC1 vs PC2
plot(pca$x[,1], pca$x[,2], xlab="PC1",ylab="PC2")

## PC1/PC2 are the most important PCs because they account for the greatest amount of variability in the dataset, but we can try other combinations of PCs as well
plot(pca$x[,1], pca$x[,3], xlab="PC1",ylab="PC3")  ##PC1 vs PC3
plot(pca$x[,4], pca$x[,5], xlab="PC4",ylab="PC5")  ##PC4 vs PC5

## sticking with PC1 vs PC2, we clearly see two clusters but what do they represent? to figure this out, we'll need to color/label our points to see if there's any column in the metadata that explains our clusters.  we'll use a very popular visualization package called ggplot 
install.packages("ggplot2")
library(ggplot2)

## create a data frame for ggplot, including the 2 PCs we want to plot (PC1 and PC2) as well as the metadata information that we want to color our points by. Here we'll first see if the two clusters on PC1 and PC2 correspond to the sexes (male v female)
ggplot.dataframe=data.frame(PC_on_xaxis=pca$x[,1], PC_on_yaxis=pca$x[,2],group=metadata$Sex, id=metadata$Sample.ID) ##note group=metadata$sex
ggplot(data=ggplot.dataframe, aes(x=PC_on_xaxis, y=PC_on_yaxis, colour=group)) +
  theme_bw()+
  geom_point(size=7) + 
  xlab("PC1") + 
  ylab("PC2")

## there's clearly a split on the y-axis (PC2) with respect to sex. all females have a high PC2 value while all males have a low PC2 value. in other words, sex is on PC2. However, this doesn't explain the two larger clusters we see on PC1 ie the left and right sides of the graph. trying a different variable from our metadata, it seems that these clusters are a result of the different heart tissues: atrial vs ventricle. in other words, 'tissue' is on PC1
ggplot.dataframe=data.frame(PC_on_xaxis=pca$x[,1], PC_on_yaxis=pca$x[,2],group=metadata$Tissue, id=metadata$Sample.ID)  ##note now group=metadata$tissue
ggplot(data=ggplot.dataframe, aes(x=PC_on_xaxis, y=PC_on_yaxis, colour=group)) +
  theme_bw()+
  geom_point(size=7) + 
  xlab("PC1") + 
  ylab("PC2")





#############################################
## HEIRARCHIAL CLUSTERING (DENDROGRAMS) #####
#############################################
##**## heirarchial clustering is another very common way to cluster gene expression data.  here, instead of transforming expression into principle components, we're going to measure the distance (or dissimilarity) between the 99 input samples and create a dendogram or tree ##**##


## what's the euclidean distance between the first two samples?
counts.normalized[,1:2]  #expression for first two samples
dist(t(counts.normalized[,1:2]))  #distance or dissimilarity between the first two samples

## now the euclidean distance between all 99 samples
dist.out=dist(t(counts.normalized))

## use the hclust() function to perform heirarchial clustering on these distances
hcluster.out <- hclust(dist.out)

## we can plot this; each 'line' is called a branch and the samples at the bottom are called 'leaves'
plot(hcluster.out)

## this dendrogram is hart to understand, let's use ggplot to create a better figure. we'll need a new package installed called ggdendro 
install.packages("ggdendro")
library(ggdendro)

## create the ggdendro object using our hcluster.out from the prev step
ggplot_dendrogram <- dendro_data(hcluster.out, type="rectangle")

## we again want to color our points based on a column in our metadata. to do this, we first need to see which samples in the dendrogram correspond to our metadata. since we saw in the PCA that the main source of variation (PC1) is tissue, let's set that as the variable of interest. we should see a clear split with all ventricle samples on one side and atrial samples on the other
indexes.group=match(as.character(ggplot_dendrogram$labels$label), metadata$Sample.ID)
ggplot_dendrogram$labels$group=metadata$Tissue[indexes.group]  ##note using metadata$tissue
ggplot() + 
  geom_segment(data=ggplot_dendrogram$segments, aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=ggplot_dendrogram$labels, aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
  geom_point(data=ggplot_dendrogram$labels, aes(x, y, color=group), size=4) + 
  ylim(-max(ggplot_dendrogram$segments$y)*(1/5), max(ggplot_dendrogram$segments$y))+
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.y=element_blank(), axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())


## now we'll annotate by sex. since we saw in our PCA analysis that sex was on PC2, we should see some sub-clustering
ggplot_dendrogram$labels$group=metadata$Sex[indexes.group]  ##note now using metadata$sex
ggplot() + 
  geom_segment(data=ggplot_dendrogram$segments, aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=ggplot_dendrogram$labels, aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
  geom_point(data=ggplot_dendrogram$labels, aes(x, y, color=group), size=4) + 
  ylim(-max(ggplot_dendrogram$segments$y)*(1/5), max(ggplot_dendrogram$segments$y))+
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.y=element_blank(), axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())





####################################################
## DIFFERENTIAL EXPRESSION ANALYSIS WITH DESEQ #####
####################################################
## instead of just visualizing expression via PCAs and dendrograms, we can perform linear modeling to determine genes which are significantly different between two conditions of interest. here we'll use DESEq2 to identify genes which are significantly different between left ventricle and atrial appendage tissue

## provide deseq2 the raw counts (not normalized), the metadata, and the variable of interest (in this case, Tissue)
de1.dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~Tissue)
de1.deseq <- DESeq(de1.dds)

## now extract the results with the argument contrasts(<variable of interest>, <condition1>, <condidion2>)
de1.resultsTable=as.data.frame(results(de1.deseq, contrast=c("Tissue", "Heart - Atrial Appendage", "Heart - Left Ventricle")))
View(de1.resultsTable)

## plot comparing log2fold change with pvalue (also called a volcano plot)
ggplot(de1.resultsTable, aes(x=log2FoldChange, y=(-log10(padj)))) +
  geom_point()+
  theme_bw()+
  ylab("p-value (-log10)")+xlab("log2 Fold Expression")+
  geom_hline(yintercept=1.3, lty=2)+
  theme(legend.position="none")

## top most up- and down-regulated genes
de1.resultsTable[tail(order(de1.resultsTable$log2FoldChange), 5),]  ##genes increased expression in left ventricle
de1.resultsTable[head(order(de1.resultsTable$log2FoldChange), 5),]  ##genes decreased expression in left ventricle

## how many total genes are significantly differentially expressed?
length(which(de1.resultsTable$padj<0.05))  #12,301 DE genes

## 12,301 genes is alot, too many for pathway enrichment. There are a few ways to get fewer genes
sig_DEgenes=rownames(de1.resultsTable)[head(order(de1.resultsTable$padj, decreasing=F), 250)]                      #filter for the 250 most significant
sig_DEgenes=rownames(de1.resultsTable)[which(de1.resultsTable$padj<=0.05&abs(de1.resultsTable$log2FoldChange)>2)]  #filter for those with > 2 or < -2 fold change



#####################################################
## PATHWAY ENRICHMENT ANALYSIS WITH ENRICHR ########
#####################################################
library(enrichR)
View(enrichr(sig_DEgenes, "KEGG_2021_Human")[[1]])
View(enrichr(sig_DEgenes, "GO_Biological_Process_2021")[[1]])
##list of possible pathway databases other than KEGG and GO: https://maayanlab.cloud/Enrichr/#libraries




