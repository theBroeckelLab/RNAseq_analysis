########################################################################
#### Co-expression network analysis tutorial                         ###
#### 09.26.2023                                                      ###
####                                                                 ###
#### Three Parts                                                     ###
#### -1. Simple co-expression network analysis with WGCNA            ###
#### -2. Differential co-expression network analysis with DGCA       ###
#### -3. Single-sample co-expression network analysis with lionessR  ###
########################################################################







##############################################################
####  1. Simple co-expression network analysis with WGCNA ####
##############################################################

## lit to review
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/


## install/load packages
#BiocManager::install("WGCNA")
library(WGCNA)
library(DESeq2)


## read in AA U01 expression data
metadata= read.table("Z:/Projects/Project Management/iPSC-CM U01 RNA sequencing/AA CEU Data to Nik/MCW-U01_RNAseq_counts/AA/U01_RNAseq_AA_GRCh38_metadata.txt", row.names=1, header=T, sep="\t")
 counts=read.table("Z:/Projects/Project Management/iPSC-CM U01 RNA sequencing/AA CEU Data to Nik/MCW-U01_RNAseq_counts/AA/U01_RNAseq_AA_GRch38_rawCounts.txt", sep="\t", row.names = 1, header=T, stringsAsFactors=F)
colnames(counts) <- rownames(metadata)


## keep only the Unstim samples
idx.unstim=which(metadata$Condition=="Unstim")
metadata=metadata[idx.unstim,]
counts=counts[,idx.unstim]


## confirm there are at least 15 samples to run co-expression network analysis
if(nrow(metadata)<15) {stop("too few samples for WGCNA")}


## keep only genes with median read count >= 10 reads
counts=counts[apply(counts,1,median)>=10,]

## normalize the data using vst from DESEq2 package
counts.vst=varianceStabilizingTransformation(as.matrix(counts))

## remove batch effect in data (optional)
#counts.vst=removeBatchEffect(counts.vst, batch=metadata$sex)

## keep only the top N most variable genes (optional)
N=500
counts.vst=counts.vst[head(order(rowVars(counts.vst), decreasing=T), N),]

## set parameters for WGCNA - see publication
merge_thresh=0.2
minModSize=30
network_type="unsigned"
beta=NA



## find optimal value for beta parameter
if (network_type=="signed") {powers = 6:30}
if (network_type=="unsigned") {powers = 3:15}
sft=pickSoftThreshold(t(counts.vst), powerVector=powers, networkType=network_type)
sft_idx=(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
beta=min(sft$fitIndices$Power[(which(sft_idx>=0.80))])
ggplot(data.frame(power=sft$fitIndices$Power, corr=sft_idx), aes(x=power, y=corr, label=power))+
  geom_text(size=5, color=ifelse(sft$fitIndices$Power==beta, 'red','black'))+
  xlab("Soft Threshold (power)")+
  ylab("Scale Free Topology Model Fit")+
  geom_hline(yintercept=0.80, color="steelblue",lty=2)



## create co-expression adjacency matrix (NxN matrix with every gene-gene correlation)
wgcna.adj.out=adjacency(t(counts.vst), power = beta, type=network_type)
wgcna.adj.out[1:5,1:5]

## transform adjacency matrix to a toplogical overlap measure (TOM) matrix
wgcna.tom.out=TOMsimilarity(wgcna.adj.out, TOMType=network_type)

## Use TOM matrix to identify co-expression modules
dissTOM=(1-wgcna.tom.out)
geneTree=hclust(as.dist(dissTOM), method = "average")
dynamicMods=cutreeDynamic(dendro=geneTree, distM = dissTOM, deepSplit=2, pamRespectsDendro=FALSE, minClusterSize=minModSize)
merge = mergeCloseModules(t(counts.vst), dynamicMods, cutHeight=merge_thresh)
moduleLabels=merge$colors
modules=data.frame(gene=rownames(counts.vst), module=paste0("M", moduleLabels))
modules[1:5,]
#there are 3 modules of co-expressed genes (along with M0, which contains genes that aren't co-expressed)
table(modules$module)

## calculate module eigengenes (the 1st principle component of a module's expression)
wgcna.eigengene.out=merge$newMEs
wgcna.eigengene.out


## calculate connectivity (all of a gene's gene-gene correlations summed)
k=intramodularConnectivity(wgcna.adj.out, moduleLabels)[,1:3]
k[1:5,]
#connectivity should follow a power-law distribution
hist(k$kOut, breaks=25) ##total connectivity of a gene
hist(k$kOut, breaks=25) ##connectivity of a gene only considering other genes in it's module
hist(k$kOut, breaks=25) ##connectivity of a gene only considering genes not in it's module


## compile modules and connectivity
wgcna.gene.out=data.frame(Modules=paste0("M", moduleLabels), k)
wgcna.gene.out[1:5,]



## save data
save(file="Z:/Projects/Project Management/test.RData", wgcna.gene.out, counts, counts.vst, metadata, wgcna.eigengene.out)


## clear environment and free memory
rm(list=ls())
gc()






################################################################
####  2. Differential expression network analysis with DGCA ####
################################################################

## lit to review
#https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-016-0349-1
#https://htmlpreview.github.io/?https://github.com/andymckenzie/DGCA/blob/master/inst/doc/DGCA.html


## install/load packages
install.packages("DGCA")
library(DGCA)
library(DESeq2)


## read in AA U01 expression data
metadata= read.table("Z:/Projects/Project Management/iPSC-CM U01 RNA sequencing/AA CEU Data to Nik/MCW-U01_RNAseq_counts/AA/U01_RNAseq_AA_GRCh38_metadata.txt", row.names=1, header=T, sep="\t")
counts=read.table("Z:/Projects/Project Management/iPSC-CM U01 RNA sequencing/AA CEU Data to Nik/MCW-U01_RNAseq_counts/AA/U01_RNAseq_AA_GRch38_rawCounts.txt", sep="\t", row.names = 1, header=T, stringsAsFactors=F)
colnames(counts) <- rownames(metadata)

## define a variable of interest to run comparison (ie for Unstim vs ET1, variable of interest is 'Condition')
variable_ofInterest="Condition"

## confirm there are at least 15 samples IN EACH GROUP to run co-expression network analysis
table(metadata[[variable_ofInterest]])
if(any(as.numeric(table(metadata[[variable_ofInterest]]))<15)) {stop("too few samples for WGCNA")}

## keep only genes with median read count >= 10 reads
counts=counts[apply(counts,1,median)>=10,]

## normalize the data using vst from DESEq2 package
counts.vst=varianceStabilizingTransformation(as.matrix(counts))

## remove batch effect in data (optional)
#counts.vst=removeBatchEffect(counts.vst, batch=metadata$sex)

## keep only the top N most variable genes (optional)
N=500
counts.vst=counts.vst[head(order(rowVars(counts.vst), decreasing=T), N),]



## create a meta matrix for DGCA
dgca_meta <- matrix(0, ncol=2, nrow=ncol(counts.vst))
rownames(dgca_meta)=colnames(dgca_counts)
colnames(dgca_meta)=c("Unstim","ET1")
dgca_meta[which(metadata$Condition=="Unstim"),"Unstim"]=1
dgca_meta[which(metadata$Condition=="ET1"),"ET1"]=1



## Run DGCA
genesToProcess=rownames(counts.vst)
dgca.dcls=list()
dgca.gene.out=data.frame(Gene=genesToProcess)
#loop through each gene and process with DGCA
for (i in 1:length(genesToProcess)) {
  if(i%%100==0) {print(paste0("processing gene...", i, " of ", length(genesToProcess)))}
  #run main differential function ddcorAll()
  d.out=suppressMessages(ddcorAll(inputMat=counts.vst, design=dgca_meta,
                                  compare=c("Unstim", "ET1"), nPairs="all", nPerm=0,
                                  adjust="BH", getDCorAvg=F, classify=FALSE, corrType="pearson",
                                  dCorAvgType="gene_average", dCorAvgMethod="mean",
                                  splitSet=genesToProcess[i], verbose=F))
  #save gene N's full results as dataframe in list
  dgca.dcls[[genesToProcess[i]]]=d.out
  #summarize gene N's results in dataframe
  dgca.gene.out$Sig_No.Total[i]=length(which(d.out$pValDiff_adj<=0.05))
  dgca.gene.out$All_Median.ABS.ZDiff[i]=median(abs(d.out$zScoreDiff))
}


## Review results - finding differentially co-expressed links (DCLs)
dgca.dcls$USP9Y  ## DGCA results for USP9Y - no DCLs between Unstim and ET1
dgca.dcls$NPPB  ## DGCA results for NPPB - 258 DCLs (pValDiff_adj<0.05) between Unstim and ET1
hist(dgca.gene.out$Sig_No.Total, breaks=20) ## ~150 genes have 0 DCLs
View(dgca.gene.out)


## Plot one gene X gene expression to understand DCLs
#higly sig diff correlation between NPPB and NRK
dgca.dcls$NPPB[1,]
#pull expression from counts.vst and plot
geneXgene=data.frame(Condition=rep(metadata$Condition,2), Gene=rep(c("NPPB","NRK"), each=nrow(metadata)),
                     NPPB.expression=as.numeric(counts.vst[which(rownames(counts.vst)=="NPPB"),]),
                     NRK.expression=as.numeric(counts.vst[which(rownames(counts.vst)=="NRK"),]))
#no correlation expression between NPPB and NRK in Unstim samples
ggplot(geneXgene[which(geneXgene$Condition=="Unstim"),], aes(x=NPPB.expression, y=NRK.expression, color=Condition))+
  geom_point(size=5)+theme_bw()
#Strong negitave expression correlation between NPPB and NRK in ET1 samples
ggplot(geneXgene[which(geneXgene$Condition=="ET1"),], aes(x=NPPB.expression, y=NRK.expression, color=Condition))+
  geom_point(size=5)+theme_bw()



## Idenfity differentially co-expressed genes (DCGS) which genes that have more DCLs than expected by chance
#NRK has 211 DCLs, use binomial model to determine if this is significant
prob=(sum(dgca.gene.out$Sig_No.Total)/2)/(((nrow(counts.vst)^2)/2)-(nrow(counts.vst)/2))
for (i in 1:nrow(dgca.gene.out)) {dgca.gene.out$Binom.Pval[i]=pbinom(dgca.gene.out$Sig_No.Total[i], nrow(counts.vst)-1, prob, lower.tail=FALSE) }
dgca.gene.out$Binom.Qval=p.adjust(dgca.gene.out$Binom.Pval, method="bonferroni")
#NRK's 211 DCLs is significant, NRK is a DCG
dgca.gene.out[which(dgca.gene.out$Gene=="NRK"),]
#in total, 236 genes have more DCLs than expected by chance (aka 236 DCGs)
length(which(dgca.gene.out$Binom.Qval<=0.05))



## Save output
dgca.out=list(dgca.dcls=dgca.dcls, dgca.dcgs=dgca.gene.out)
save(dgca.dcls, dgca.gene.out, metadata, counts, counts.vst, dgca_meta, file="Z:/Projects/Project Management/test.RData")


## clear environment and free memory
rm(list=ls())
gc()







#####################################################################
####  3. Single-sample expression network analysis with lionessR ####
#####################################################################

## review lit
#https://www.sciencedirect.com/science/article/pii/S2589004219300872
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7870822/
#https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-6235-7


## install packages
#BiocManager::install("lionessR")
library(lionessR)

## read in synthego data
meta=read.table("Z:/Projects/Project Management/Synthego/Data/pilot 20 targets - runs1-2/metadata.txt", sep="\t", header=T, row.names=NULL)
meta$TargetGene[which(meta$TargetGene=="Trac")]="TRAC"
counts=read.table("Z:/Projects/Project Management/Synthego/Data/pilot 20 targets - runs1-2/counts.txt", sep="\t", header=T, row.names = 1)

## normalize counts
counts.vst=varianceStabilizingTransformation(as.matrix(counts))

## get names of KO target genes
ko.genes=setdiff(unique(meta$TargetGene), c("TRAC","Mock"))

## isolate the N most variable genes 
N=500
topVar.genes=rownames(counts.vst)[head(order(as.numeric(apply(counts.vst, 1, var)), decreasing=T), N)]

## combine most variable genes with the KO target genes and filter counts
analysisGenes=unique(c(ko.genes, topVar.genes))
counts.vst=counts.vst[which(rownames(counts.vst)%in%analysisGenes),]

## run lionessR for one replicate of one KO
idx.target=which(meta$TargetGene=="MYH7"&meta$ReplicateID=="replicate_1"&meta$Plate=="plate_1")
idx.controls=which(meta$TargetGene=="Mock")
lioness_oneSample=lioness(counts.vst[,c(idx.target, idx.controls)])

## analyze results (MYH7-rep1-plate1 in first column)
dat.out_oneSample=as.data.frame(assay(lioness_oneSample))
View(dat.out_oneSample)
#three gene-gene pairs with the strongest positive relationship in MYH7-rep1-plate1 sample
dat.out_oneSample[head(order(dat.out_oneSample$MYH7_P1B1R1_S125, decreasing=T), 3),]
#three gene-gene pairs with the strongest negitave relationship in MYH7-rep1-plate1 sample
dat.out_oneSample[head(order(dat.out_oneSample$MYH7_P1B1R1_S125, decreasing=F), 3),]



## loop through each KO gene- instead of using 1 replicate, get z-scores for all 6 replicates then average together
z.compiled=data.frame(rep(NA, nrow(counts.vst)^2))
lioness.out=list()
for (j in 1:length(ko.genes)) {
  print(paste0("processing ", ko.genes[j], " (", j, " of 20)"))
  idx.koGene=which(meta$TargetGene==ko.genes[j])
  for (i in 1:length(idx.koGene)) {
    idx.inSamples=c(which(meta$TargetGene=="Mock"), idx.koGene[i])
    dat.out=lioness(counts.vst[,idx.inSamples])
    if (i==1) {lioness.out[[j]]=data.frame(assay(dat.out)[,ncol(assay(dat.out))]); next}
    lioness.out[[j]]=data.frame(cbind(lioness.out[[j]], data.frame(assay(dat.out)[,ncol(assay(dat.out))])))
  }
  colnames(lioness.out[[j]])=meta$sample[idx.koGene]
  z.compiled[,j]=rowMeans(lioness.out[[j]])
}
rownames(z.compiled)=rownames(lioness.out[[1]])
names(lioness.out)=ko.genes
colnames(z.compiled)=ko.genes


## compile all z-scores from every lionessR run
all.zs=as.numeric(as.matrix(z.compiled))

## confirm that they are normally distrubuted
hist(all.zs, breaks=50, xlab="Delta Correlation")
ks.test(all.zs, "pnorm") # Kolmogorov-Smirnov test will indicate whether data is normal (p<=0.05) or not

## if data is normal than any z-score > mean+(2*sd) or < mean-(2*sd) is significant
z.upper= (mean(all.zs, na.rm=T)) + (2*sd(all.zs, na.rm=T))
z.lower= (mean(all.zs, na.rm=T)) - (2*sd(all.zs, na.rm=T))

## identify sig gene-gene pairs with the single MYH7-rep1-plate1 sample from above
sig.pos=rownames(dat.out_oneSample)[which(dat.out_oneSample$MYH7_P1B1R1_S125>z.upper)]  ##sig positive gene-gene changes
sig.neg=rownames(dat.out_oneSample)[which(dat.out_oneSample$MYH7_P1B1R1_S125<z.lower)]  ##sig negitave gene-gene changes


## 44,800 of 260,100 gene-gene pairs (17.2%) are correlated in MYH7-rep1-plate1 sample
sum(length(sig.pos), length(sig.neg))
nrow(dat.out_oneSample)
100*(sum(length(sig.pos), length(sig.neg))) / nrow(dat.out_oneSample)


## save data
save(counts, counts.vst, meta, z.compiled, all.zs, file="Z:/Projects/Project Management/test.RData")



## clear environment and free memory
rm(list=ls())
gc()
