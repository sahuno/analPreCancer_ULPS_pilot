
#### Date; August 19th 2020
### fuction to eextract ichorCNA paramters 
##author, Samuel Ahuno

library(reshape2)
library(hablar)
library(lubridate)
library(tidyverse)
library(stringr)
library(data.table)


#### get params from google bucket
#gsutil cp 'gs://fc-d305704f-d91b-41da-a3bf-1103dc9fee0b/FFPE_ichor/*BRP*.params.txt' .

##### paramters
#filePath <- "/Users/samuelahuno/Polak_lab10082019/Ghana_cfDNA/IchorCNA_vr2_RepTime_Results/results_ichor_minSegBin20_AltFractThresh_0.01/ichorCNA/GBR06256_1000kb/GBR06256_1000kb.params.txt"
# fread("/Users/samuelahuno/Polak_lab10082019/Ghana_cfDNA/IchorCNA_vr2_RepTime_Results/results_ichor_minSegBin20_AltFractThresh_0.01/ichorCNA/GBR06256_1000kb/GBR06256_1000kb.params.txt"
#       , header=FALSE,  col.names = c("parameter","value"),sep = "\t", strip.white = TRUE,blank.lines.skip = TRUE, nrows = 9,skip = 5)


IchorCNA_param_extract <- function(param_text_file=""){
  ###read param files
  read.ichor.param <- fread(param_text_file, header=FALSE,  col.names = c("parameter","value"),sep = "\t", strip.white = TRUE,blank.lines.skip = TRUE, nrows = 9,skip = 5)
  ###get file names, use later
  sampleName <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(param_text_file)))
  
  ##remove any non-unifrm characters
  read.ichor.param$sample <- sampleName
  read.ichor.param$parameter <- gsub(" ","_",gsub(":","",read.ichor.param$parameter))
  ###make data wide
  ichor.params <- reshape2::dcast(read.ichor.param, sample ~ parameter, value.var = "value")
  ichor.params.final <- ichor.params %>% select("sample","Tumor_Fraction","Ploidy","Subclone_Fraction","Fraction_CNA_Subclonal","Fraction_Genome_Subclonal","Gender","ChrX_median_log_ratio","ChrY_coverage_fraction","Coverage")
  return(ichor.params.final)
}

###test
#test.param <- IchorCNA_param_extract(filePath)

#spacec <- "Space with me:"
#  gsub(" ","",gsub(":","",spacec))
##initiazialize empty vector
ichor.param.dataStore <- data.frame(stringsAsFactors=FALSE)
### for loop to run all samples
for (params_txt in list.files(path="/Users/samuelahuno/Polak_lab10082019/Collabs_internal/analPreCancer_ULPS_pilot/data/raw_data/dir_params_txt",full.names=TRUE,pattern = "params.txt", recursive = TRUE)) {
  df.data <- IchorCNA_param_extract(params_txt)
  #data1$v1 <- NA
  #data1$v1 <- df.data[1]
  ichor.param.dataStore <- as.data.frame(rbind(ichor.param.dataStore,df.data))
}


#kofi <- data.frame("ll" =1:5)
#kofi %>% select("ll")

ichor.param.dataStore.copy <- ichor.param.dataStore

##coerce data to numeric 
ichor.param.dataStore.copy <- ichor.param.dataStore.copy %>% hablar::convert(num(Tumor_Fraction, Ploidy,Subclone_Fraction,Fraction_CNA_Subclonal,
                                              Fraction_Genome_Subclonal,ChrX_median_log_ratio,
                                              ChrY_coverage_fraction,Coverage),chr(Gender))


##function to convert to perentage
cal_percent <- function(x, na.rm=FALSE) (x * 100)
ichor.param.dataStore.copy <- ichor.param.dataStore.copy %>% mutate_at(vars(c("Tumor_Fraction","Subclone_Fraction","Fraction_CNA_Subclonal","Fraction_Genome_Subclonal","ChrY_coverage_fraction")),cal_percent)


##read other informatio
read.pid.infor <- fread("/Users/samuelahuno/Polak_lab10082019/Collabs_internal/analPreCancer_ULPS_pilot/data/raw_data/Anal_ichor_Filenames_v2.txt", header=TRUE,sep = "\t",na.strings = "NA")
remove_space_colNames <- gsub(" ","_",names(read.pid.infor)) ##remove spaces from names
names(read.pid.infor) <- remove_space_colNames

###some irregular file names/cahracters should be removed
read.pid.infor$sample <- gsub("\\.","-",gsub("\xca","",read.pid.infor$sample))
read.pid.infor$First_seen <- gsub("\xca","",read.pid.infor$First_seen)
read.pid.infor$Biopsy_Result <- gsub(" ","",read.pid.infor$Biopsy_Result)


##put NA everywhere  there is no data
read.pid.infor[read.pid.infor==""]<-NA 


#remove 19* from year so as to be systematic
#df.DOB <- "19/09/2919"
#str_remove(df.DOB, "19")
#str_replace(df.DOB, "19", "")
#gsub(".*\\/19","",read.pid.infor$Birthdate)  
#gsub("\\/19.*","",read.pid.infor$Birthdate) 

##no need if Tm has recoded dates with full year
# read.pid.infor$Birthdate <- sub('19(?=..$)', '', read.pid.infor$Birthdate, perl = TRUE)
# read.pid.infor$Birthdate <- sub('(?=..$)', '19', read.pid.infor$Birthdate, perl = TRUE)
# read.pid.infor$`First seen`  <- sub('(?=..$)', '20', read.pid.infor$`First seen`, perl = TRUE) #do same for first seen

##convert dates to class date
read.pid.infor$Birthdate <- mdy(read.pid.infor$Birthdate)
read.pid.infor$`First_seen` <- mdy(read.pid.infor$`First_seen`)

#tets_date <- as_date("2019/08/14")
#date_today <- today(tzone = "") ##set todays date
#interval(tets_date,date_today) %/% months(1)

##convert date of births to time
read.pid.infor <- read.pid.infor %>% mutate(Age = interval(Birthdate,First_seen) %/% years(1), monthsSinceFirstSeen = interval(First_seen,First_seen) %/% months(1))
####merge ichor results and participant information
merged.pid.infor.ichor.param.dataStore.copy <- merge(read.pid.infor,ichor.param.dataStore.copy,by.x = "sample", by.y = "sample")


###remove the 4 samples lower than 3% tumor fraction
sampleOrder_minusLow4 <- merged.pid.infor.ichor.param.dataStore.copy %>% filter(Tumor_Fraction <3.0) %>% select(sample)
sampleOrder_minusLow4 <-c("RP-2258_BRP23736_v1_WGS_OnPrem","RP-2258_BRP23748_v1_WGS_OnPrem","RP-2258_BRP23876_v1_WGS_OnPrem","RP-2258_BRP23878_v1_WGS_OnPrem")

merged.pid.infor.ichor.param.dataStore.copy <- merged.pid.infor.ichor.param.dataStore.copy[!which(merged.pid.infor.ichor.param.dataStore.copy$sample %in% sampleOrder_minusLow4),]

(merged.pid.infor.ichor.param.dataStore.copy$sample)

merged.pid.infor.ichor.param.dataStore.copy$sample <- plyr::mapvalues(merged.pid.infor.ichor.param.dataStore.copy$sample, from=c("RP-2258_BRP23739_v1_WGS_OnPrem", "RP-2258_BRP23740_v1_WGS_OnPrem", "RP-2258_BRP23743_v1_WGS_OnPrem", "RP-2258_BRP23844_v1_WGS_OnPrem", "RP-2258_BRP23846_v1_WGS_OnPrem", "RP-2258_BRP23737_v1_WGS_OnPrem", "RP-2258_BRP23845_v1_WGS_OnPrem", "RP-2258_BRP23854_v1_WGS_OnPrem", "RP-2258_BRP23734_v1_WGS_OnPrem", "RP-2258_BRP23746_v1_WGS_OnPrem", "RP-2258_BRP23747_v1_WGS_OnPrem", "RP-2258_BRP23751_v1_WGS_OnPrem", "RP-2258_BRP23863_v1_WGS_OnPrem", "RP-2258_BRP23865_v1_WGS_OnPrem", "RP-2258_BRP23872_v1_WGS_OnPrem", "RP-2258_BRP23873_v1_WGS_OnPrem"),
          to=c("P01_HSIL", "P02_HSIL", "P03_HSIL", "P04_HSIL", "P05_HSIL", "P06_HSIL", "P07_HSIL", "P08_HSIL", "P09_LSIL", "P10_LSIL", "P11_LSIL", "P12_LSIL", "P13_LSIL", "P14_LSIL", "P15_LSIL", "P16_LSIL"))


##select data for commuts
comut.anal.FISH <- merged.pid.infor.ichor.param.dataStore.copy %>% select (c("sample","FISH_3q")) %>% mutate(category="FISH")
comut.anal.lesionGrade<-merged.pid.infor.ichor.param.dataStore.copy %>% select (c("sample","Lesion_Grade")) %>% mutate(category="Lesion_Grade")

comut.anal.Gender <-merged.pid.infor.ichor.param.dataStore.copy %>% select (c("sample","Gender.x")) %>% mutate(category="Gender")
comut.anal.Biopsy <-merged.pid.infor.ichor.param.dataStore.copy %>% select (c("sample","Biopsy_Result")) %>% mutate(category="Biopsy")
comut.anal.Race <-merged.pid.infor.ichor.param.dataStore.copy %>% select (c("sample","Race")) %>% mutate(category="Race")
comut.anal.Age <- merged.pid.infor.ichor.param.dataStore.copy %>% select (c("sample","Age")) %>% mutate(category="Age")

comut.anal.Age$AgeDiagnosis <- NA
comut.anal.Age <- comut.anal.Age %>% mutate(AgeDiagnosis = ifelse(Age < 50,"below 50 years", ifelse(Age >= 50, "50 years and above", AgeDiagnosis)))


comut.anal.Tfx <- merged.pid.infor.ichor.param.dataStore.copy %>% select (c("sample","Tumor_Fraction"))

##arange by 
sampleOrder <- merged.pid.infor.ichor.param.dataStore.copy %>% arrange(Lesion_Grade,desc(Tumor_Fraction)) %>% select(sample)
sampleOrder_minusLowTFx <- merged.pid.infor.ichor.param.dataStore.copy %>% arrange(Lesion_Grade,match(FISH_3q, c("3q", "normal", "pending"))) %>% select(sample)




####create tables for HIV status
comut.anal.hiv <- merged.pid.infor.ichor.param.dataStore.copy %>% select (c("sample")) %>% mutate(category="HIV_status",values="Positive")


##commut plot table formats, use this
#write.table(ihc.only,file = "../data/ihc.only.4pandas_edited.tsv",col.names = c("sample","category","value"),quote=FALSE, sep='\t')
write.table(merged.pid.infor.ichor.param.dataStore.copy,file = "../analPreCancer_ULPS_pilot/data/output/merged.ichorResults.PID.infor.4pandas.tsv",quote=FALSE, sep='\t')
##
write.table(comut.anal.FISH[,c(1,3,2)],file = "../analPreCancer_ULPS_pilot/data/output/comut.anal.FISH.4pandas.tsv",col.names = c("sample","category","value"),quote=FALSE, sep='\t')
write.table(comut.anal.lesionGrade[,c(1,3,2)],file = "../analPreCancer_ULPS_pilot/data/output/comut.anal.lesionGrade.4pandas.tsv",col.names = c("sample","category","value"),quote=FALSE, sep='\t')
write.table(comut.anal.Biopsy[,c(1,3,2)],file = "../analPreCancer_ULPS_pilot/data/output/comut.anal.Biopsy.4pandas.tsv",col.names = c("sample","category","value"),quote=FALSE, sep='\t')
write.table(comut.anal.Gender[,c(1,3,2)],file = "../analPreCancer_ULPS_pilot/data/output/comut.anal.Gender.4pandas.tsv",col.names = c("sample","category","value"),quote=FALSE, sep='\t')
write.table(comut.anal.Race[,c(1,3,2)],file = "../analPreCancer_ULPS_pilot/data/output/comut.anal.Race.4pandas.tsv",col.names = c("sample","category","value"),quote=FALSE, sep='\t')
write.table(comut.anal.Age[,c(1,3,4)],file = "../analPreCancer_ULPS_pilot/data/output/comut.anal.Age.4pandas.tsv",col.names = c("sample","category","value"),quote=FALSE, sep='\t')
write.table(comut.anal.Tfx,file = "../analPreCancer_ULPS_pilot/data/output/comut.anal.Tfx.4pandas.tsv",col.names = c("sample","Tumor_fraction"),quote=FALSE, sep='\t')
write.table(comut.anal.hiv,file = "../analPreCancer_ULPS_pilot/data/output/omut.anal.hiv.4pandas.tsv",col.names = c("sample","category","value"),quote=FALSE, sep='\t')


merged.pid.infor.ichor.param.dataStore.copy %>% select (c("sample","Race")) %>% mutate(category="Race")


names(merged.pid.infor.ichor.param.dataStore.copy)
