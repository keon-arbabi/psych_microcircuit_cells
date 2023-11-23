
# prepare-data-cmc-rna-seq.R
# Common Mind Consortium merged counts and metadata. Code adapted from: https://github.com/CommonMindConsortium/covarr-de
suppressPackageStartupMessages({
  
  # install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
  library(synapser)
  
  library(tidyverse)
  library(data.table)
  library(compositions)
}) 

# Read raw counts from MPP and HBCC
# counts.MPP = fread(synGlet('syn21867938')$path, data.table=FALSE)
ALL_USED_IDs = c()

# Just DLPFC counts
ALL_USED_IDs = c(ALL_USED_IDs, 'syn24172729')
counts.MPP = fread(synGet('syn24172729')$path, data.table=FALSE)
rownames(counts.MPP) = counts.MPP$Geneid
counts.MPP = counts.MPP[,-c(1:6)]

# HBCC - DLFPC
ALL_USED_IDs = c(ALL_USED_IDs, 'syn21886235')
counts.HBCC = fread(synGet('syn21886235')$path, data.table=FALSE)
rownames(counts.HBCC) = counts.HBCC$Geneid
counts.HBCC = counts.HBCC[,-c(1:6)]
geneCountsMerged = cbind(counts.MPP, counts.HBCC)

# data download
downloadFile = function(id, version=NULL){
  fread(synGet(id, version = version)$path, data.table = F)
}

# Get ancestry vector calculated using gemtools
ANCESTRY_ID = 'syn17894713'
ALL_USED_IDs = c(ALL_USED_IDs, ANCESTRY_ID)
ANCESTRY.HBCC = downloadFile(ANCESTRY_ID) %>% 
  dplyr::rename(ID = 'Genotyping_Sample_ID')
ANCESTRY_ID = 'syn17346100'
ALL_USED_IDs = c(ALL_USED_IDs, ANCESTRY_ID)
ANCESTRY.MPP = downloadFile(ANCESTRY_ID) %>% 
  dplyr::rename('DNA_report..Genotyping.Sample_ID' = 'Genotyping_Sample_ID')
ANCESTRY = rbind(
  ANCESTRY.HBCC[,colnames(ANCESTRY.MPP)[-2]], 
  ANCESTRY.MPP[,colnames(ANCESTRY.MPP)[-2]]
)
# Get genotype ids from synapse. 
GENOTYPE_ID = "syn18358480"
ALL_USED_IDs = c(ALL_USED_IDs, GENOTYPE_ID)
GENOTYPE = downloadFile(GENOTYPE_ID) %>% 
  dplyr::select(Individual_ID, `Genotyping_Sample_ID`, `Exclude`) %>% 
  dplyr::inner_join(ANCESTRY) %>% 
  dplyr::filter(is.na(`Exclude`))
# Get RNASeq QCmetadata
CLINICAL_ID = "syn3354385"
clinical = downloadFile(CLINICAL_ID)
ASSAY_ID = "syn24173489"
rnaseq = downloadFile(ASSAY_ID)
metadata = right_join(
  clinical, 
  rnaseq,
)
ALL_USED_IDs = c(ALL_USED_IDs,CLINICAL_ID, ASSAY_ID)


metadata = metadata %>%
  dplyr::select(
    `Individual ID`, Institution, Cohort, `Reported Gender`, Sex, Ethnicity, 
    `Age of Death`, `PMI (in hours)`, Dx, pH, SampleID, one_of(
      'RIN','Brain_Region', 'Ribozero_Batch', 'Library_Batch', 
      'Flowcell_Batch', 'Mapped_Reads', 'Intragenic_Rate','Intronic_Rate', 
      'Intergenic_Rate','Genes_Detected', 'Expression_Profiling_Efficiency', 
      'rRNA_Rate', 'Total_Reads', 'Percent_Aligned', 'Transcripts_Detected',
      'Exclude?')
  ) %>% 
  dplyr::filter(`Brain_Region` %in% c("DLPFC"))
metadata$`Age of Death`[metadata$`Age of Death`=="90+"] = 90
metadata$`Age of Death` = as.numeric(metadata$`Age of Death`)

# Merge metadata 
METADATA = metadata %>%
  dplyr::left_join(GENOTYPE, by = c("Individual ID" = "Individual_ID")) %>%  
  dplyr::rename(
    Region = `Brain_Region`,
    PMI = `PMI (in hours)`,
    RIN = `RIN`,
    ReportExclude = `Exclude?`,
    GenotypeExclude = Exclude,
    LibraryBatch = `Library_Batch`,
    FlowcellBatch = `Flowcell_Batch`,
    RibozeroBatch = `Ribozero_Batch`,
    MappedReads = `Mapped_Reads`,
    IntragenicRate = `Intragenic_Rate`, 
    IntronicRate = `Intronic_Rate`, 
    IntergenicRate = `Intergenic_Rate`, 
    GenesDetected = `Genes_Detected`,
    ExpProfEfficiency = `Expression_Profiling_Efficiency`, 
    rRNARate = `rRNA_Rate`,
    TotalReads = `Total_Reads`, 
    AlignmentRate = `Percent_Aligned`, 
    TranscriptsDetected = `Transcripts_Detected`,    
    TranscriptsDetected = `Transcripts_Detected`,    
    Reported_Gender = `Reported Gender`,
    ageOfDeath = `Age of Death`,
    IndividualID = `Individual ID`) %>%
  dplyr::select(
    SampleID, IndividualID, Institution, Cohort,Reported_Gender, Sex, 
    Ethnicity, ageOfDeath, PMI, Dx, RIN, EV.1, EV.2, EV.3, EV.4, EV.5, 
    LibraryBatch, FlowcellBatch, RibozeroBatch, MappedReads, TotalReads, 
    GenesDetected, AlignmentRate, IntragenicRate, IntergenicRate, IntronicRate, 
    ExpProfEfficiency, rRNARate, TranscriptsDetected, ReportExclude, 
    GenotypeExclude, pH)

# data preprocessing
METADATA = METADATA %>%
  dplyr::filter(SampleID %in% colnames(geneCountsMerged), !is.na(SampleID)) %>%
  dplyr::mutate(Cohort = forcats::fct_recode(Cohort, `MSSM-Penn-Pitt` = 'MSSM-Penn-Pitt', `NIMH-HBCC`='NIMH-HBCC'))  %>%
  dplyr::filter(!Dx %in% c("undetermined")) 
#dplyr::filter(ageOfDeath > 17) 
unique(METADATA$Dx)

paste("BP =", nrow(METADATA %>% filter(Dx == "BP")))
paste("AFF =", nrow(METADATA %>% filter(Dx == "AFF")))
paste("SCZ =", nrow(METADATA %>% filter(Dx == "SCZ")))
paste("Control =", nrow(METADATA %>% filter(Dx == "Control")))

# set Control to baseline
METADATA$Dx = factor(METADATA$Dx, c('Control','SCZ','BP'))

ind = METADATA$SampleID[which(METADATA$ReportExclude == 1 | METADATA$GenotypeExclude)]
writeLines(paste('Following',length(ind),'samples are marked exclude'))
writeLines(paste(ind, collapse = ', '))
METADATA = METADATA  %>% dplyr::filter(!(SampleID %in% ind)) 

ind = METADATA$SampleID [is.na(METADATA$Ethnicity) | is.na(METADATA$Institution) | is.na(METADATA$Dx)]
writeLines(paste('Following', length(ind), 'counts are missing any metadata'))
writeLines(paste(ind, collapse = ', '))
METADATA = METADATA  %>% dplyr::filter(!(SampleID %in% ind)) 

ind = METADATA$SampleID [is.na(METADATA$PMI)]
writeLines(paste('Following', length(ind), 'counts are missing PMI'))
writeLines(paste(ind, collapse = ', '))
METADATA = METADATA  %>% dplyr::filter(!(SampleID %in% ind)) 

ind = METADATA$SampleID [is.na(METADATA$Reported_Gender)]
writeLines(paste('Following', length(ind), 'counts are missing gender'))
writeLines(paste(ind, collapse = ', '))
METADATA = METADATA  %>% dplyr::filter(!(SampleID %in% ind)) 

ind = METADATA$SampleID [is.na(METADATA$ageOfDeath)]
writeLines(paste('Following', length(ind), 'counts are missing age of death'))
writeLines(paste(ind, collapse = ', '))
METADATA = METADATA  %>% dplyr::filter(!(SampleID %in% ind))

ind = METADATA$SampleID [is.na(METADATA$EV.1)]
writeLines(paste('Following', length(ind), 'counts are missing ancestry information'))
writeLines(paste(ind, collapse = ', '))
METADATA = METADATA  %>% dplyr::filter(!(SampleID %in% ind))

# drop individuals who where sequenced twice here
tab = table(METADATA$Individual_ID)

dropSamples = sapply( names(tab[tab > 1]), function(id){
  
  idx = which(METADATA$Individual_ID == id)
  df = METADATA[idx,]
  
  # remove individuals that have less then max reads
  i = which(df$MappedReads < max(df$MappedReads))
  rownames(METADATA)[idx[i]]
} )
METADATA = METADATA[!(rownames(METADATA) %in% dropSamples),]

# include estimated cell fractions in METADATA
df_cellFractions = read.table(synGet('syn22333694')$path, row.names=1)

celFrac_ilr = ilr(df_cellFractions)
colnames(celFrac_ilr) = paste0("cellFrac_ilr_", 1:3)

METADATA = merge(METADATA, df_cellFractions, by.x = 'SampleID', by.y='row.names')
METADATA = merge(METADATA, celFrac_ilr, by.x = 'SampleID', by.y='row.names')

# Match covariates to expression data
indToRetain = intersect(METADATA$SampleID, colnames(geneCountsMerged))

geneCountsMerged = geneCountsMerged[,indToRetain]
saveRDS(geneCountsMerged, file = "./data/commonmind/geneCountsMerged.rds")

rownames(METADATA) = METADATA$SampleID
saveRDS(METADATA, file = "./data/commonmind/METADATA.rds")

