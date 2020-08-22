

###### Copy number annoations for Seg files comign from broad, ichor cna run 


##########  run the following codes functions  #### below ##########
####use hg19 ensemsbl dataset to annotate   ###############################
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl.hg19 <- getBM(c("chromosome_name", "start_position", "end_position","strand","hgnc_symbol","band","ensembl_gene_id","gene_biotype"),
                      filters = c("chromosome_name"),
                      values = list(chromosome_name=c(1:22, "X", "Y", "M")),
                      mart = grch37)

##rename cytobands, ie `chr+band``
ensembl.hg19$band <- paste0(ensembl.hg19$chromosome_name,ensembl.hg19$band)
###correct strand information for 
ensembl.hg19.strandChnaged <- ensembl.hg19
ensembl.hg19.strandChnaged$strand <- gsub("-1","-",ensembl.hg19.strandChnaged$strand)
ensembl.hg19.strandChnaged$strand <- gsub("1","+",ensembl.hg19.strandChnaged$strand)

##final Granges Object of ensembl genes hg19
ensembl.hg19.strandChnaged <- makeGRangesFromDataFrame(ensembl.hg19.strandChnaged,start.field= "start_position", end.field= "end_position",keep.extra.columns = T)

##function to annotate cna ###########
AnnotateEnsembleGenes.anal <- function(geneList=ensembl.hg19.strandChnaged,filePath=""){
  df.seg <- read.delim(filePath, header=TRUE, sep = "\t", strip.white = TRUE)
  dr.seg <- makeGRangesFromDataFrame(df.seg,start.field=  "start", end.field= "end", seqnames.field="chr", keep.extra.columns = T)
  sampleName <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(filePath)))
  ##find overlaps between cna.seg and ensemble genes
  Ov.refGenes.in.CNA.dr.seg <- findOverlaps(query=ensembl.hg19.strandChnaged, subject=dr.seg, type = "within") ##complete overlaps
  unique(queryHits(Ov.refGenes.in.CNA.dr.seg))
  queryLength(Ov.refGenes.in.CNA.dr.seg)  ##sanity checks
  
  ###assigning copy number values to gene hits ##same as previous line
  ensembl.hg19.strandChnaged$copy_number <- NA
  ##Get Subject hits as dataFrame
  sub_hit <- as.data.frame(dr.seg[subjectHits(Ov.refGenes.in.CNA.dr.seg)])
  ensembl.hg19.strandChnaged[queryHits(Ov.refGenes.in.CNA.dr.seg)]$copy_number <- sub_hit[,grepl("copy.number",names(sub_hit))]
  annotated_Seg_EnsembleGenes <- ensembl.hg19.strandChnaged
  
  #print(tools::file_path_sans_ext(tools::file_path_sans_ext(basename(filePath))))
  #print(annotated_Seg_EnsembleGenes[seqnames(annotated_Seg_EnsembleGenes)=="10", ][85:95])
  #print(annotated_Seg_EnsembleGenes[which(annotated_Seg_EnsembleGenes$Gene.name=="PTEN")])
  #print(annotated_Seg_EnsembleGenes[which(annotated_Seg_EnsembleGenes$Gene.name=="CCND1")])
  write.table(annotated_Seg_EnsembleGenes, file=paste0("../analPreCancer_ULPS_pilot/data/output/",sampleName,".annotated.","ensemblGenes.tsv"), quote=FALSE, sep='\t', col.names = NA)
  
  ####write 
  df.metadata <- as.data.frame(elementMetadata(annotated_Seg_EnsembleGenes))
  df.sample <- as.data.frame(df.metadata[,grepl("copy_number",names(df.metadata))])
  geneNames <- as.data.frame(df.metadata[,grepl("hgnc_symbol",names(df.metadata))])
  cytobands <- as.data.frame(df.metadata[,grepl("band",names(df.metadata))])
  
  
  df.output <- data.frame(geneNames,df.sample,cytobands)
  names(df.output)[1] <- "Gene_Names"
  names(df.output)[2] <- sampleName
  names(df.output)[3] <- "cytoband"
  return(df.output)
}