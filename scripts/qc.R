
# differential-expression-analysis.R 
# load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(edgeR)
  library(RUVSeq)
  library(psych)
  library(stats)
  library(scales)
  library(CovariateAnalysis) #devtools::install_github('th1vairam/CovariateAnalysis@dev')
  library(mvIC) #remotes::install_github("GabrielHoffman/mvIC")
  library(variancePartition)
  library(svglite)
  library(gridExtra)
  library(ggpubr)
  library(ggsci)
  library(ggrepel)
})
cols_1 = c("PyrL2n3" = "#D62728FF", "PyrL5n6" = "#FF7F0EFF", "PVALB" = "#2CA02CFF", "SST" = "#9467BDFF", "VIP" = "#1F77B4FF")

# load data ----
load(file = "./data/lcm-seq/data.Rdata")
# clean-up medications 
metadata = metadata %>%
  mutate(Tobacco = factor(ifelse(grepl("Y", Tob.ATOD),"Y","N"))) %>%
  mutate(Antidepressants = factor(ifelse(grepl("D", Meds.ATOD),"Y","N"))) %>%
  mutate(Benzos_Anticonvulsants = factor(ifelse(grepl(paste0(c("B","C"), collapse = "|"), Meds.ATOD),"Y","N"))) %>%
  mutate(Antipsychotics = factor(ifelse(grepl("P", Meds.ATOD),"Y","N"))) %>%
  mutate(Lithium = factor(ifelse(grepl("L", Meds.ATOD),"Y","N"))) 
# add star qc 
metadata = metadata %>% left_join(., star_qc, by = 'ID')
rownames(metadata) = metadata$ID

# covariates clustering ----
# determine relationship between covariates
factor_covars = c("HU","DX","CT","MOD","Sex","Race","Pool")
cont_covars = c("Age","PMI","pH","RIN","Intronic_rate","Intragenic_rate","Intergenic_rate","Library_size","Detected_genes",
                "Tobacco","Antidepressants","Benzos_Anticonvulsants","Antipsychotics","Lithium"
                #"Uniquely_mapped_pct","Multi_mapped_pct","Unmapped_pct"
                )

# Find inter relation between covariates
covariates = metadata[,c(factor_covars, cont_covars), drop = F] %>%
  # clean up names
  dplyr::rename(Donor_ID = HU, Diagnosis = DX, Cell_type = CT, Tissue_RIN = RIN, Tissue_pH = pH, Tissue_PMI = PMI)
rownames(covariates) = metadata$ID  
# correlation/association between covariates at an FDR <= 0.1
covariates_correlation = getAssociationStatistics(covariates, PVAL = 0.05)
tmp = covariates_correlation$ESTIMATE
#tmp[covariates_correlation$PVAL > 0.05] = 0
h = Heatmap(tmp, col = colorRamp2(c(-1,0,1), c('blue','white','red')), name = 'AssocEstimate')

#svglite(file = "./figures/draft/covariate_associations_heatmap.svg", width = 7, height = 5.8)
ComplexHeatmap::draw(h, heatmap_legend_side = 'left')
#dev.off()

# explore metadata ----
my.theme = theme_bw() %+replace% theme(legend.position = 'top', plot.title = element_text(hjust = 0.5, face = "bold"))
p = list()
# Age of death
p[[1]] = ggplot(covariates, aes(x = Cell_type, y = Age)) + geom_boxplot() + stat_compare_means(vjust = 1, hjust = 0)
p[[1]] = p[[1]] + ggtitle('Age of Death') + my.theme
# RIN
p[[2]] = ggplot(covariates, aes(x = Cell_type, y = RIN)) + geom_boxplot() + stat_compare_means(vjust = 1, hjust = 0)
p[[2]] = p[[2]] + ggtitle('RIN') + my.theme
# PMI
p[[3]] = ggplot(covariates, aes(x = Cell_type, y = PMI)) + geom_boxplot() + stat_compare_means(vjust = 1, hjust = 0)
p[[3]] = p[[3]] + ggtitle('PMI (in hours)') + my.theme
# Intronic Rate
p[[4]] = ggplot(covariates, aes(x = Cell_type, y = Intronic_rate)) + geom_boxplot() + stat_compare_means(vjust = 1, hjust = 0)
p[[4]] = p[[4]] + ggtitle('Intronic Rate') + my.theme
# Intergenic Rate
p[[5]] = ggplot(covariates, aes(x = Cell_type, y = Intergenic_rate)) + geom_boxplot() + stat_compare_means(vjust = 1, hjust = 0)
p[[5]] = p[[5]] + ggtitle('Intergenic Rate') + my.theme
# Genes Detected
p[[6]] = ggplot(covariates, aes(x = Cell_type, y = Detected_genes)) + geom_boxplot() + stat_compare_means(vjust = 1, hjust = 0)
p[[6]] = p[[6]] + ggtitle('Genes Detected') + my.theme
# Library Size
p[[7]] = ggplot(covariates, aes(x = Cell_type, y = Library_size)) + geom_boxplot() + stat_compare_means(vjust = 1, hjust = 0)
p[[7]] = p[[7]] + ggtitle('Library Size') + my.theme
# Percent Unmapped Reads
p[[8]] = ggplot(covariates, aes(x = Cell_type, y = Unmapped_pct)) + geom_boxplot() + stat_compare_means(vjust = 1, hjust = 0)
p[[8]] = p[[8]] + ggtitle('Unmapped Reads (percent)') + my.theme
# Percent Multi-mapped Reads
p[[9]] = ggplot(covariates, aes(x = Cell_type, y = Multi_mapped_pct)) + geom_boxplot() + stat_compare_means(vjust = 1, hjust = 0)
p[[9]] = p[[9]] + ggtitle('Multimapped Reads (percent)') + my.theme

#svglite(file = "./figures/draft/covariate_compare_ct.svg", width = 9, height = 10)
multiplot(plotlist = p, cols = 3)
#dev.off()

# outlier analysis ----
# keep genes that are expressed in at least one combination of CT and DX
dge0 = DGEList(counts_matrix)
dge0 = calcNormFactors(dge0, method = "TMM")

filt_lst = by(metadata, list(metadata$CT, metadata$DX), function(x){
  rowSums(cpm(counts_matrix[,x$ID]) > 1) >= 15
})
keep = colSums(do.call(rbind, filt_lst)) > 0 
# filter expression 
table(keep)
dge = dge0[keep,]

# find principal components of expression to plot
PC = prcomp(voom(dge, plot = TRUE)$E, scale. = T, center = T)
plotdata = data.frame(ID = rownames(PC$rotation), 
                      PC1 = PC$rotation[,1], 
                      PC2 = PC$rotation[,2])
# percentage from each PC
summary(PC)$importance[2, 1:5, drop = FALSE]
eigen = PC$sdev^2
pc1 = eigen[1]/sum(eigen)
pc2 = eigen[2]/sum(eigen)
# identify outliers - samples 4 SDs from the mean
outliers = as.character(plotdata$ID[c(which(plotdata$PC1 < mean(plotdata$PC1) - 4*sd(plotdata$PC1)), which(plotdata$PC1 > mean(plotdata$PC1) + 2*sd(plotdata$PC1))), drop = T])
outliers = c(outliers, as.character(plotdata$ID[c(which(plotdata$PC2 < mean(plotdata$PC2) - 4*sd(plotdata$PC2)), which(plotdata$PC2 > mean(plotdata$PC2) + 2*sd(plotdata$PC2))), drop = T] ))

# plot first 2 PCs
plotdata = left_join(plotdata, metadata, by = "ID") %>%
  dplyr::mutate(label = ID) %>% 
  dplyr::mutate(label = ifelse((label %in% outliers), label, NA)) 
ggplot(plotdata, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = CT, shape = DX, size = scale(Intergenic_rate))) +
  geom_text_repel(aes(label = label), size = 4, hjust = 0) +
  scale_color_manual(values = cols_1) + 
  theme_classic2() + 
  theme(legend.position = "right") #+ facet_grid(CT~.)
# save
#ggsave(file = "./figures/draft/sample_outliers_pca.svg", height = 6, width = 7)

# Plot aberrant distribution of logcpm counts
voom(dge, plot = FALSE)$E %>%
  rownameToFirstColumn('gene_ensembl') %>%
  tidyr::gather(ID, logCPM, -gene_ensembl) %>%
  left_join(., metadata, by = "ID") %>% 
  filter(ID %in% outliers) %>%
ggplot(aes(x = logCPM, color = ID)) + 
  geom_density() + 
  theme_classic2() + 
  theme(legend.position = 'top') + 
  facet_grid("DX", scale = 'free')
#ggsave(filename = "./figures/draft/outlier_expression.svg", height = 8, width = 10)

# significant covariates ----
# find correlation between PC's of gene expression with covariates across CTs
covariates_select = covariates %>% dplyr::select(-c(Library_size, Detected_genes))
preAdjustedSigCovars = suppressWarnings(
  runPCAandPlotCorrelations(genesBySamples = voom(dge, plot = FALSE)$E, 
                            samplesByCovariates =  covariates_select, 
                            CORRELATION_TYPE = "spearman",
                            dataName = "",
                            isKeyPlot = TRUE, 
                            MIN_PVE_PCT_PC = 2))
preAdjustedSigCovars$significantCovars
svglite(file = "./figures/draft/covariate_pca_corr.svg", width = 7, height = 4)
preAdjustedSigCovars[["PC_res"]][[2]]$plotData
dev.off()

# model identification ----
# forward stepwise regression with the mvIC criterion enables automated variable selection for high dimensional datasets
# https://github.com/GabrielHoffman/mvIC

# keep genes that are expressed in at least one combination of CT and DX
dge0 = DGEList(counts_matrix)
dge0 = calcNormFactors(dge0, method = "TMM")

model_lst = by(metadata, list(metadata$CT), function(x){
  # filter within each CT
  counts = dge0$counts[, x$ID]
  keep = rowSums(cpm(counts) > 1) >= 15
  dge = DGEList(counts[keep,], group = x$DX)
  
  # get unwanted variation using RUVseq
  design = model.matrix(~0 + DX, x)
  y = estimateGLMCommonDisp(dge, design) 
  y = estimateGLMTagwiseDisp(y, design) 
  fit = glmFit(y, design) 
  res = residuals(fit, type = "deviance")
  
  w_vars = RUVr(dge$counts, rownames(dge$counts), k=5, res)
  covars = cbind(covariates[x$ID,], w_vars$W)
  
  #print(cor(covars$W_1, covars$Intergenic_rate))

  # set base formula
  base_formula = glue::glue("~ ", glue::glue_collapse("Diagnosis", sep = " + "))
  # step-wise regression
  phase1 = mvForwardStepwise(exprObj = voom(dge, plot = FALSE)$E,
                             baseFormula = base_formula,
                             data = covars,
                             criterion = "BIC",
                             variables = array(c("(1|MOD)","(1|Sex)","(1|Race)","(1|Tobacco)","scale(Age)",
                                                 "(1|Antidepressants)","(1|Benzos_Anticonvulsants)","(1|Antipsychotics)","(1|Lithium)"
                                                 )))
  
  # fit a base model with the following covariates
  phase2 = mvForwardStepwise(exprObj = voom(dge, plot = FALSE)$E,
                             baseFormula = phase1$formula, 
                             data = covars,
                             criterion = "BIC",
                             variables = array(c("scale(Intronic_rate)","scale(Intergenic_rate)","(1|Pool)","scale(Tissue_PMI)","scale(Tissue_pH)","scale(Tissue_RIN)",
                                                 "W_1","W_2","W_3","W_4","W_5",
                                                 "scale(Library_size)","scale(Detected_genes)"
                                                 )))
  knitr::kable(phase2$settings)
  knitr::kable(phase2$trace)
  as.character(phase2$formula)[2]
})

#write_rds(model_lst, file = "./output/model_lst.rds")
model_lst = readRDS(file = "./output/model_lst.rds")

  
