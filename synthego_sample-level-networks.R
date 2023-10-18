################################################################################
#### Sample-level network analysis with Synthego KOs                         ###
#### 10.16.2023                                                              ###
####                                                                         ###
################################################################################

## review lit
#https://www.sciencedirect.com/science/article/pii/S2589004219300872
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7870822/
#https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-6235-7

## install packages
#BiocManager::install("lionessR")
library(lionessR)
library(igraph)
library(enrichR)
library(DESeq2)





#######################################################
## run lionessR for one single replicate of one KO ####
#######################################################


## read in synthego data
meta=read.table("Z:/Projects/Project Management/Synthego/Data/pilot 20 targets - runs1-2/metadata.txt", sep="\t", header=T, row.names=NULL)
counts=read.table("Z:/Projects/Project Management/Synthego/Data/pilot 20 targets - runs1-2/counts.txt", sep="\t", header=T, row.names = 1)

## 6 replicates for 20 KO genes, 1 negative mock control, and 1 positive (TRAC) control
table(meta$TargetGene, meta$SeqRun)


## get names of KO target genes
ko.genes=setdiff(unique(meta$TargetGene), c("TRAC","Mock"))
ko.genes

## normalize counts
counts.vst=varianceStabilizingTransformation(as.matrix(counts))


## isolate the N most variable genes for analysis
N=500
topVar.genes=rownames(counts.vst)[head(order(as.numeric(apply(counts.vst, 1, var)), decreasing=T), N)]


## filter counts data for only the N most variable genes (make sure KO target genes aren't filtered out)
analysisGenes=unique(c(ko.genes, topVar.genes))
counts.vst=counts.vst[which(rownames(counts.vst)%in%analysisGenes),]


## pull the indexes for 6 mock controls
idx.controls=which(meta$TargetGene=="Mock")
meta[idx.controls,]

## pull the indexes for replicate_1 plate_1 of MYH7
idx.target=which(meta$TargetGene=="MYH7"&meta$ReplicateID=="replicate_1"&meta$Plate=="plate_1") 
meta[idx.target,]

## Filter normalized counts for 6 controls and 1 MYH7 KO replicate
counts.vst_oneSample=counts.vst[,c(idx.target, idx.controls)]

## Run lioness
lioness_oneSample=as.data.frame(assay(lioness(counts.vst_oneSample)))
##**## you'll likely see a number of warning messages, this is 
##*because we typically want 15+ samples to run lionessR, however
##*we're only running 7 with the current analysis. this is a 
##*major drawback/consideration when interpreting results 

## View z-score results - each column is a sample-level network and each row is a gene-gene pair, data represent z-scores transformed from correlations
View(lioness_oneSample)

## although we get a network for each input sample, we're mostly interested in the MYH7 sample (not the 6 control sample networks)
View(data.frame(genes=rownames(lioness_oneSample), MYH7_P1B1R1_S125=lioness_oneSample$MYH7_P1B1R1_S125))


## in the MYH7 KO sample, we can review only the MYH7 z-scores 
hist(lioness_oneSample$MYH7_P1B1R1_S125)                                               #histogram of z-scores for MYH7 gene-gene correlations
lioness_oneSample[head(order(lioness_oneSample$MYH7_P1B1R1_S125, decreasing=T), 3),]   #three gene-gene pairs with the strongest positive MYH7 relationship
lioness_oneSample[head(order(lioness_oneSample$MYH7_P1B1R1_S125, decreasing=F), 3),]   #three gene-gene pairs with the strongest negative MYH7 relationship



## or we can use all gene-gene pairs in the MYH7 sample for network visualization, first reformat data into long format
rownames(lioness_oneSample)=gsub("NM_", "NM", rownames(lioness_oneSample))
rownames(lioness_oneSample)=gsub("NR_", "NR", rownames(lioness_oneSample))
df.net1=data.frame(gene1=sapply(strsplit(rownames(lioness_oneSample), "_"), `[`, 1),
                   gene2=sapply(strsplit(rownames(lioness_oneSample), "_"), `[`, 2),
                   corr1=lioness_oneSample$MYH7_P1B1R1_S125, 
                   network.partner=0)
df.net1[1:5,]


## when do we consider two genes as partners in the network?  Here, any with z-score 2 standard deviations from the average z-score
sd.factor=2
limit.upper=mean(df.net1$corr1, na.rm=T)+(sd.factor*sd(df.net1$corr1, na.rm=T))     #z-scores greater than this are considered partners
df.net1$network.partner[which(df.net1$corr1>=limit.upper)]=1                        #set df.net1$network.partner=1 for any gene-gene z-scores > upper limit
df.net1[147394:147398,]                                                             #ACTA1 and MYH7 are partners with pos z-score
limit.lower=mean(df.net1$corr1, na.rm=T)-(sd.factor*sd(df.net1$corr1, na.rm=T))     #z-scores less than this are also considered partners
df.net1$network.partner[which(df.net1$corr1<=limit.lower)]=1                        #set df.net1$network.partner=1 for any gene-gene z-scores < lower limit
df.net1[552:558,]                                                                   #C1S and ACE2 are partners with neg z-score


## create the graph/network object using package igraph
me.utm=graph_from_data_frame(df.net1[which(df.net1$network.partner!=0),1:2], directed=F)
me.utm

## by default- gene-gene partners are counted twice (ie C1S-ACE2 and ACE2-C1S)
df.net1[c(555,22442),]

## removes these redundancies with simplify()
me.utm=simplify(me.utm)

## calculate network genes, edges, and density 
length(V(me.utm))                    #network has 373 genes (we started with ~500 most var genes, those with 0 partners were filtered out)
length(E(me.utm))                    #network has 5,927 gene-gene partnerings
length(E(me.utm))/length(V(me.utm))  #average gene has ~16 partners


## genes with only 1 partner can make network visualization difficult, often best to remove
sum(degree(me.utm)<=1)  #22 genes with only 1 partner
me.utm=delete.vertices(me.utm, degree(me.utm)<=1)


## re-calculate network genes, edges, and density 
length(V(me.utm))                    #network has 351 genes (genes with 0 partners have been filtered out)
length(E(me.utm))                    #network has 5,905 gene-gene partnerings
length(E(me.utm))/length(V(me.utm))  #average gene has ~17 partners


## choose colors for network visualization
cols=c("blue3","darkgreen","orangered","grey40")

## initially set all genes to grey
V(me.utm)$color=rep(cols[4], length(V(me.utm)))

## change the KO gene to blue
V(me.utm)$color[which(V(me.utm)$name%in%"MYH7")]=cols[1]

## identify MYH7 partner genes with a positive relationship, change to green
limit.upper_genes=df.net1$gene2[which(df.net1$gene1=="MYH7")][which(df.net1$corr1[which(df.net1$gene1=="MYH7")]>limit.upper)]
V(me.utm)$color[which(V(me.utm)$name%in%limit.upper_genes)]=cols[2]

## identify MYH7 partner genes with a negative relationship, change to red
limit.lower_genes=df.net1$gene2[which(df.net1$gene1=="MYH7")][which(df.net1$corr1[which(df.net1$gene1=="MYH7")]<limit.lower)]
V(me.utm)$color[which(V(me.utm)$name%in%limit.lower_genes)]=cols[3]


## by default- genes are plotted from their order in network thus the first few genes will be 'buried' under all the others in the network vis
V(me.utm)[1:9]

## reorder the network genes so that grey genes are plotted first, followed by red/green, then blue KO gene on top
v.reo=c(); for (k in 1:length(cols)) {v.reo=c(v.reo, which(V(me.utm)$color==cols[k]))}
me.utm=permute(me.utm, match(V(me.utm)$name, V(me.utm)$name[rev(v.reo)]))

## render plot - igraph has tons of parameters for visualization (see https://r.igraph.org/)
print(plot.igraph(me.utm,
                  main="MYH7",
                  vertex.size=4,
                  vertex.color=V(me.utm)$color,
                  vertex.frame.color=V(me.utm)$color,
                  vertex.label="",
                  xlim=c(-1,1),ylim=c(-1,1),
                  edge.width=0.1,
                  layout=layout.fruchterman.reingold(me.utm)))  #there are many types of network layouts, fruchterman-reingold usually best for biological data


## final network summary
length(V(me.utm))                                                          #number of network genes
length(E(me.utm))                                                          #number of network edges (ie partners)
round(length(E(me.utm)) / length(V(me.utm)))                               #average partners per gene
length(which(c(limit.upper_genes, limit.lower_genes)%in%V(me.utm)$name))   #number of partners for MYH7
length(which(limit.upper_genes%in%V(me.utm)$name))                         #number of positive MYH7 partners
length(which(limit.lower_genes%in%V(me.utm)$name))                         #number of negative MYH7 partners

## run enrichment with the MYH7 partners
myh7.partners=c(limit.upper_genes, limit.lower_genes)
enrich.out=enrichr(myh7.partners, "KEGG_2021_Human")[[1]]
View(enrich.out)

## pull all gene-gene partners in network
all.network.partners=as_data_frame(me.utm, what="edges")
all.network.partners[1:10,]

##number of partners per gene
network.partner.counts=as.data.frame(table(as.character(as.matrix(all.network.partners))))
network.partner.counts[1:5,]                                                     #number of partners for the first 5 genes
network.partner.counts[order(network.partner.counts$Freq, decreasing=T)[1:3],]   #TNNI3 (a known hypertrophic gene) has the most partners, ie a network hub
network.partner.counts[order(network.partner.counts$Freq, decreasing=T)[17],]    #MYH7 has the 17th most partners

## network should have a power-law distribution characterized by a minority of hub genes and a majority of peripheral genes 
hist(table(as.character(as.matrix(as_data_frame(me.utm, what="edges")))), breaks=20, ylab="number of genes",xlab="number of edges", main="")


## clear environment and free memory
rm(list=ls())
gc()













##########################################################
## loop lionessR to use all 6 replicates of all 20 KO ####
##########################################################

## read in synthego data
meta=read.table("Z:/Projects/Project Management/Synthego/Data/pilot 20 targets - runs1-2/metadata.txt", sep="\t", header=T, row.names=NULL)
counts=read.table("Z:/Projects/Project Management/Synthego/Data/pilot 20 targets - runs1-2/counts.txt", sep="\t", header=T, row.names = 1)

## 6 replicates for 20 KO genes, 1 negative mock control, and 1 positive (TRAC) control
table(meta$TargetGene, meta$SeqRun)


## get names of KO target genes
ko.genes=setdiff(unique(meta$TargetGene), c("TRAC","Mock"))
ko.genes

## normalize counts
counts.vst=varianceStabilizingTransformation(as.matrix(counts))


## isolate the N most variable genes 
N=500
topVar.genes=rownames(counts.vst)[head(order(as.numeric(apply(counts.vst, 1, var)), decreasing=T), N)]


## filter for N most variable genes (make sure KO target genes aren't filtered out)
analysisGenes=unique(c(ko.genes, topVar.genes))
counts.vst=counts.vst[which(rownames(counts.vst)%in%analysisGenes),]


## create empty data frame for z-scores, columns are KO genes and rows will be all possible gene-gene pairs
z.compiled=data.frame(matrix(nrow=nrow(counts.vst)^2, ncol=length(ko.genes), dimnames=list(1:nrow(counts.vst)^2, ko.genes)))

## create empty list; it will have 20 data frames (one for each KO) each with 6 columns (one for each KO replicate)
lioness.out=list()

## loop through each KO (1 through 20)
for (j in 1:length(ko.genes)) {
  #pull index for KO gene replicates
  idx.koGene=which(meta$TargetGene==ko.genes[j])
  ## loop through each replicate (1 through 6)
  for (i in 1:length(idx.koGene)) {
    print(paste0("processing gene ", j, " of ", length(ko.genes), ": ", ko.genes[j], " (replicate ", i, " of 6)"))
    #pull 6 control samples + replicate i of KO gene j
    idx.inSamples=c(which(meta$TargetGene=="Mock"), idx.koGene[i])
    #run lionessR
    dat.out=lioness(counts.vst[,idx.inSamples])
    #add lionessR results for replicate i to the lioness.out list
    if (i==1) {lioness.out[[j]]=data.frame(assay(dat.out)[,ncol(assay(dat.out))]); next}
    lioness.out[[j]]=data.frame(cbind(lioness.out[[j]], data.frame(assay(dat.out)[,ncol(assay(dat.out))])))
  }
  #each dataframe in lioness.out list contains the z-scores for all 6 replicates of KO gene j
  colnames(lioness.out[[j]])=meta$SampleName[idx.koGene]
  lioness.out[[j]][1:5,]
  #take the mean z-score across all 6 replicates, add to z.compiled data frame
  z.compiled[,j]=rowMeans(lioness.out[[j]])
}
rownames(z.compiled)=rownames(lioness.out[[1]])
names(lioness.out)=ko.genes


## review z.compiled- this includes the mean z-score across 6 replicates for each KO gene
z.compiled[1:5,]


## compile all z-scores from every lionessR run
all.zs=as.numeric(as.matrix(z.compiled))
hist(all.zs, breaks=50, xlab="Delta Correlation")



## create empty data frame to store network stats
network.stats=data.frame(matrix(nrow=length(ko.genes), ncol=5, dimnames=list(ko.genes, c("net.nodes","net.density","ko.partners_tot","ko.partners_up","ko.partners_down"))))
network.stats

## loop through each KO and generate network using igraph (see above for more detailed notes)
for (i in 1:length(ko.genes)) {
  #igraph visualizations
  goi=ko.genes[i]
  rownames(z.compiled)=gsub("NM_", "NM", rownames(z.compiled))
  rownames(z.compiled)=gsub("NR_", "NR", rownames(z.compiled))
  df.net1=data.frame(gene1=sapply(strsplit(rownames(z.compiled), "_"), `[`, 1),
                     gene2=sapply(strsplit(rownames(z.compiled), "_"), `[`, 2),
                     corr1=z.compiled[[goi]], corr1.thresholded=0)
  sd.factor=2
  limit.upper=mean(all.zs, na.rm=T)+(sd.factor*sd(all.zs, na.rm=T))                  
  limit.lower=mean(all.zs, na.rm=T)-(sd.factor*sd(all.zs, na.rm=T))                   
  df.net1$corr1.thresholded[which(df.net1$corr1>=limit.upper)]=1
  df.net1$corr1.thresholded[which(df.net1$corr1<=limit.lower)]=1
  me.utm=graph_from_data_frame(df.net1[which(df.net1$corr1.thresholded!=0),1:2], directed=F)
  me.utm=simplify(me.utm)
  me.utm=delete.vertices(me.utm, degree(me.utm)<=1)
  ##set colors
  cols=c("blue3","darkgreen","orangered","grey40")
  V(me.utm)$color=rep(cols[4], length(V(me.utm)))
  V(me.utm)$color[which(V(me.utm)$name%in%goi)]=cols[1]
  limit.upper_genes=df.net1$gene2[which(df.net1$gene1==goi)][which(df.net1$corr1[which(df.net1$gene1==goi)]>limit.upper)]
  limit.lower_genes=df.net1$gene2[which(df.net1$gene1==goi)][which(df.net1$corr1[which(df.net1$gene1==goi)]<limit.lower)]
  V(me.utm)$color[which(V(me.utm)$name%in%limit.upper_genes)]=cols[2]
  V(me.utm)$color[which(V(me.utm)$name%in%limit.lower_genes)]=cols[3]
  v.reo=c(); for (k in 1:length(cols)) {v.reo=c(v.reo, which(V(me.utm)$color==cols[k]))}
  me.utm=permute(me.utm, match(V(me.utm)$name, V(me.utm)$name[rev(v.reo)]))
  l=layout.fruchterman.reingold(me.utm)
  print(plot.igraph(me.utm,
                    main=ko.genes[i],
                    vertex.size=4,
                    vertex.color=V(me.utm)$color,
                    vertex.frame.color=V(me.utm)$color,
                    vertex.label="",
                    xlim=c(-1,1),ylim=c(-1,1),
                    edge.width=0.1,
                    layout=l))
  #network summary
  network.stats$net.nodes[i]=length(V(me.utm))
  network.stats$net.density[i]=round(length(E(me.utm)) / length(V(me.utm)))
  network.stats$ko.partners_tot[i]=length(which(c(limit.upper_genes, limit.lower_genes)%in%V(me.utm)$name)) 
  network.stats$ko.partners_up[i]=length(which(limit.upper_genes%in%V(me.utm)$name))
  network.stats$ko.partners_down[i]=length(which(limit.lower_genes%in%V(me.utm)$name))
}


## view/write network stats
network.stats
write.xlsx(network.stats, "Z:/Projects/Project Management/Synthego/networkStats.xlsx")


## OPTIONAL: save network data
save(counts, counts.vst, meta, z.compiled, all.zs, file="Z:/Projects/Project Management/Synthego/test.RData")


## clear environment and free memory
rm(list=ls())
gc()








