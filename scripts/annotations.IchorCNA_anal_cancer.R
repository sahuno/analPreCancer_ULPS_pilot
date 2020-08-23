#### Anal cancer ananlysis

library(biomaRt)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationHub)
library(data.table)
library(GenomeInfoDb)
library(pheatmap)
library(plot.matrix)
library(data.table)
library(tidyverse)
library(reshape2)

##comfig, set the following parameters


#path of CNA files
cna.anal.paths <- "../analPreCancer_ULPS_pilot/data/raw_data/cna_seg/"
source("../analPreCancer_ULPS_pilot/scripts/ichorCNA_annotations_ULPSBroad.R") ###load sccript to extract copy number from *cna.seg

dataAnal <- data.frame(stringsAsFactors=FALSE)

### use this to annotate all samples together
for (dat in list.files(path=cna.anal.paths,full.names=TRUE,pattern = "cna.seg", recursive = TRUE)) {
  df.data <- AnnotateEnsembleGenes.anal(geneList=gr.ensemblGenes.biomart.hg19,filePath=dat)
  #data1$v1 <- NA
  #data1$v1 <- df.data[1]
  dataAnal <- as.data.frame(c(dataAnal,df.data))
}

##get gene names
getGeneName<- dataAnal[,grepl("Gene",names(dataAnal))][1]
getCytobands <- dataAnal[,grepl("cytoband",names(dataAnal))][1]

removeMatch <- c("Gene","cytoband") ### patterns to remove
##get cna nunbers
df.annotations.genes <- cbind(getGeneName,getCytobands,dataAnal[,!grepl(paste(removeMatch, collapse = "|"),names(dataAnal))])

###insert NA's into empty cells and get ridd of all empty cases get
df.annotations.genes[df.annotations.genes==""]<-NA
df.annotations.genes <- df.annotations.genes[complete.cases(df.annotations.genes), ]

##get unique  cells 
df.annotations.genes <- distinct(df.annotations.genes, Gene_Names, .keep_all = TRUE)

##add band to gene names
df.annotations.genes.loci  <- df.annotations.genes 
df.annotations.genes.loci$gene_loci <- NA
df.annotations.genes.loci$geneBand <- paste0(df.annotations.genes.loci$Gene_Names," (",df.annotations.genes.loci$cytoband,")")
df.saved.gene.loci <- df.annotations.genes.loci[,c(1,2,ncol(df.annotations.genes.loci))] ### save gene name and locis, use later to interect
#rownames(df.annotations.genes.loci) <- df.annotations.genes.loci[,1] 
#df.annotations.genes.loci <- df.annotations.genes.loci[-1]


rownames(df.annotations.genes) <- df.annotations.genes[,1] 
df.annotations.genes.final <- df.annotations.genes[-1]

#names(df.annotations.genes.final)



##comut ready data - for plotting
#melt.df.annotations.genes <-  melt(head(df.annotations.genes[-2]),id="Gene_Names",id.name="Gene",variable.name = "sample",value.name ="corrected_copy_number") ##test with header
melt.df.annotations.genes.anal <-  melt(df.annotations.genes[-2],id="Gene_Names",id.name="Gene",variable.name = "sample",value.name ="corrected_copy_number")
#melt.df.annotations.genes.anal.v2 <-  melt(df.annotations.genes.loci[-2],id="Gene_Names",id.name="Gene",variable.name = "sample",value.name ="corrected_copy_number")

melt.df.annotations.genes.anal %>% filter(Gene_Names=="ASC")

# top.cna.genes.anal <- c("FGFR3", "PDGFRA", "CSF1R","TP53",
#                         "FLT3","EGFR","APC","KDR","HRAS","IDH1",
#                         "RET","SMARCB1","STK11","CDH1","NOTCH1","ERBB2",
#                         "PTEN","PIK3CA","CDKN2A","HNF1A","MPL","ALK","BRAF","KMT2D")

##subset selected genes for comut
anal_cancer_specific_genes <- str_sort(c("CCL22","STK11","DNMT3B","TERC","PIK3CA","TP63","TGFBR2","TERT","DUSP4","ATM","BAP1","TRAF3"), numeric = TRUE)
sub.melt.df.annotations.cna.genes.anal = melt.df.annotations.genes.anal[which(melt.df.annotations.genes.anal$Gene_Names %in% anal_cancer_specific_genes),]

melt.df.annotations.genes.anal[which(melt.df.annotations.genes.anal$Gene_Names %in% "TGFBR2"),]
sub.melt.df.annotations.cna.genes.anal[which(sub.melt.df.annotations.cna.genes.anal$Gene_Names %in% "TGFBR2"),]
melt.df.annotations.genes.anal[which(melt.df.annotations.genes.anal$Gene_Names %in% "ATM"),]

###sumamry stats on 
sub.melt.df.annotations.cna.genes.summary.anal <- sub.melt.df.annotations.cna.genes.anal %>% group_by(Gene_Names) %>% dplyr::summarise(count=n(),minCN=min(corrected_copy_number),
                                                                                                                             medianCN=median(corrected_copy_number),maxCN=max(corrected_copy_number),
                                                                                                                             CNgain=sum(corrected_copy_number>2 & corrected_copy_number<4),
                                                                                                                             CNgain_Percent=((sum(corrected_copy_number>2 & corrected_copy_number<4))/n())*100,
                                                                                                                             CNamp_Percent=((sum(corrected_copy_number==4))/n())*100,
                                                                                                                             CNhamp_Percent=((sum(corrected_copy_number>4))/n())*100,
                                                                                                                             CNneutr_Percent=((sum(corrected_copy_number==2))/n())*100,
                                                                                                                             CNdel_Perc=((sum(corrected_copy_number<2))/n())*100,
                                                                                                                             CN_anyCN_above1_Percent=((sum(corrected_copy_number>1))/n())*100,
                                                                                                                             CN_atLeast_gain=(sum(corrected_copy_number>2)),
                                                                                                                             CN_atLeast_gain_Percent=((sum(corrected_copy_number>2))/n())*100)

###plot percentage of individuals with any cn higher tan neutral in the cohort
#ggplot(sub.melt.df.annotations.cna.genes.summary.anal, aes(x=as.factor(Gene_Names),y=CN_anyCN_above1_Percent)) + geom_col(color="blue")



sub.melt.df.annotations.cna.genes.anal$cnaEvent <- NA
test.cna.event.anal <- sub.melt.df.annotations.cna.genes.anal %>% mutate(cnaEvent = ifelse(corrected_copy_number==2,"Neutral",ifelse(corrected_copy_number == 1, "Homodeletion",ifelse(corrected_copy_number == 4, "Amplification",ifelse(corrected_copy_number > 4, "High amplification",
                                                                                                                                                                                                                                    ifelse(corrected_copy_number == 3, "Gain", cnaEvent))))))
### errm when i add one more doesn't work anymore
test.cna.event.anal %>% filter(cnaEvent =="High amplification")
test.cna.event.anal.cmut <- test.cna.event.anal[,c(2,1,4)]

test.cna.event.anal.cmut$sample <- gsub("\\.","-",test.cna.event.anal.cmut$sample) #standardize names 


#remove low tfx samples
sampleOrder_minusLow4 <-c("RP-2258_BRP23736_v1_WGS_OnPrem","RP-2258_BRP23748_v1_WGS_OnPrem","RP-2258_BRP23876_v1_WGS_OnPrem","RP-2258_BRP23878_v1_WGS_OnPrem")
test.cna.event.anal.cmut.minus4samples <- test.cna.event.anal.cmut %>% filter(!sample %in% c("RP-2258_BRP23736_v1_WGS_OnPrem","RP-2258_BRP23748_v1_WGS_OnPrem","RP-2258_BRP23876_v1_WGS_OnPrem","RP-2258_BRP23878_v1_WGS_OnPrem"))

test.cna.event.anal.cmut.minus4samples <- merge(test.cna.event.anal.cmut.minus4samples,df.saved.gene.loci[,c(1,3)],by.x="Gene_Names",by.y="Gene_Names")

###change sample names to p01-16
test.cna.event.anal.cmut.minus4samples$sample <- plyr::mapvalues(test.cna.event.anal.cmut.minus4samples$sample, from=c("RP-2258_BRP23739_v1_WGS_OnPrem", "RP-2258_BRP23740_v1_WGS_OnPrem", "RP-2258_BRP23743_v1_WGS_OnPrem", "RP-2258_BRP23844_v1_WGS_OnPrem", "RP-2258_BRP23846_v1_WGS_OnPrem", "RP-2258_BRP23737_v1_WGS_OnPrem", "RP-2258_BRP23845_v1_WGS_OnPrem", "RP-2258_BRP23854_v1_WGS_OnPrem", "RP-2258_BRP23734_v1_WGS_OnPrem", "RP-2258_BRP23746_v1_WGS_OnPrem", "RP-2258_BRP23747_v1_WGS_OnPrem", "RP-2258_BRP23751_v1_WGS_OnPrem", "RP-2258_BRP23863_v1_WGS_OnPrem", "RP-2258_BRP23865_v1_WGS_OnPrem", "RP-2258_BRP23872_v1_WGS_OnPrem", "RP-2258_BRP23873_v1_WGS_OnPrem"),
                                                                      to=c("P01_HSIL", "P02_HSIL", "P03_HSIL", "P04_HSIL", "P05_HSIL", "P06_HSIL", "P07_HSIL", "P08_HSIL", "P09_LSIL", "P10_LSIL", "P11_LSIL", "P12_LSIL", "P13_LSIL", "P14_LSIL", "P15_LSIL", "P16_LSIL"))


write.table(test.cna.event.anal.cmut.minus4samples[,c(2,4,3)], file="../analPreCancer_ULPS_pilot/data/output/coMutPlot_ready_data_allSamples.anal.tsv",quote=FALSE, col.names = c("sample", "category","value"),sep='\t')


###find deleted genes  #####  not  needed
names(df.annotations.genes.final)

del.RP.2258_BRP23736_v1_WGS_OnPrem <-  df.annotations.genes.final %>% filter(RP.2258_BRP23736_v1_WGS_OnPrem == 1)
del.RP.2258_BRP23736_v1_WGS_OnPrem %>% filter(del.RP.2258_BRP23736_v1_WGS_OnPrem, grepl("13",del.RP.2258_BRP23736_v1_WGS_OnPrem$cytoband))

del.RP.2258_BRP23736_v1_WGS_OnPrem[which(del.RP.2258_BRP23736_v1_WGS_OnPrem$cytoband == "^13"),]

grepl()
test.cna.event.anal.cmut %>% filter(sample== "23736")
unique(test.cna.event.anal.cmut.minus4samples$Gene_Names)

df.saved.gene.loci %>% filter(Gene_Names %in% c("CCL22","STK11","DNMT3B","TERC","PIK3CA","TP63","TGFBR2","TERT","DUSP4","ATM","BAP1","TRAF3")) %>% select(geneBand)
names(df.saved.gene.loci)
geneSort <- str_sort(c("CCL22","STK11","DNMT3B","TERC","PIK3CA","TP63","TGFBR2","TERT","DUSP4","ATM","BAP1","TRAF3"), numeric = TRUE)






