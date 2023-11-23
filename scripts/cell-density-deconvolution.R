
# cell-density-deconvolution.r
# load libraries ----
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(markerGeneProfile)
  library(lme4)
  library(tidyverse)
  library(Seurat)
  library(ggbeeswarm)
  library(ggsignif)
  library(cowplot)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  library(broom)
  library(stats)
  library(WGCNA)
})

scales::show_col(pal_d3("category20c")(9))
scales::show_col(pal_d3("category20b")(9))

cols_1 = c("CTRL" = "#3182BDFF", "MDD" = "#E6550DFF", "BD" = "#31A354FF", "SCZ" = "#756BB1FF")

# load lcm-seq data ---
load(file = "./data/lcm-seq/data.Rdata")

# load CommonMind data ----
cmc_metadata = readRDS(file = "./data/cmc/METADATA.rds") %>% 
  # filter to Pitt institution only
  filter(Institution == "Pitt") %>% 
  # remove empty
  dplyr::select(-c("ReportExclude","GenotypeExclude")) %>%
  mutate(DX = factor(Dx, levels = c("Control","BP","SCZ"), labels = c("CTRL","BD","SCZ"))) %>% 
  mutate(Reported_Gender = factor(Reported_Gender, levels = c("Male","Female")))
# matching cmc counts matrix 
cmc_counts_matrix = readRDS(file = "./data/cmc/geneCountsMerged.rds") %>% dplyr::select(row.names(cmc_metadata))
# check matching rows and columns 
all(rownames(cmc_metadata) == colnames(cmc_counts_matrix)) 

# for ENSEMBL id ENSG00000279457.4, return ENSG00000279457
trimEnsembl = function(x){
  gsub("(.*)\\.(.*)", "\\1", x) 
}
# get non-unique gene names 
cmc_counts_matrix = cmc_counts_matrix %>% 
  rownames_to_column(var = "ensembl_id_unique") %>% 
  mutate(ensembl_id = trimEnsembl(ensembl_id_unique)) %>%
  dplyr::relocate(ensembl_id, .after = ensembl_id_unique)

# remove duplicate genes by selecting the gene with highest expression
rownames(cmc_counts_matrix) = cmc_counts_matrix$ensembl_id_unique
select_rows = collapseRows(datET = cmc_counts_matrix[,-c(1:2)],
                           rowGroup = cmc_counts_matrix$ensembl_id,
                           rowID = cmc_counts_matrix$ensembl_id_unique,
                           method = "MaxMean")
cmc_counts_matrix = cmc_counts_matrix[select_rows$selectedRow,]
rownames(cmc_counts_matrix) = cmc_counts_matrix$ensembl_id
cmc_counts_matrix = cmc_counts_matrix[,-c(1:2)]

# filter genes 
keep = rowSums(cpm(cmc_counts_matrix) > 0.5) >= .1*ncol(cmc_counts_matrix)
table(keep)
cmc_counts_matrix_filt = cmc_counts_matrix[keep,]

# compute RNA-seq residuals ----
# estimate effective library size using TMM
dge = DGEList(cmc_counts_matrix_filt)
dge = calcNormFactors(dge)
# get normalized matrix 
cmc_norm_matrix = as.data.frame(voom(dge)$E)
cmc_norm_matrix = rownames_to_column(cmc_norm_matrix, var = "ensembl_id")

# hold out sex, age, and PMI
design = model.matrix(~DX + Reported_Gender + scale(ageOfDeath)+ scale(PMI), cmc_metadata)
vm = voom(dge, design, plot = TRUE)
fit = lmFit(vm, design)
fit = eBayes(fit)
# plot final fit
plotSA(fit, main = "Final model: Mean-variance trend")

# get residuals 
cmc_residuals_matrix = as.data.frame(residuals(fit, cmc_norm_matrix[,-1]))
# add back DX
add_vars = grep("DX", colnames(fit$design), value = T)
cmc_residuals_matrix = cmc_residuals_matrix +
  fit$coefficients[,add_vars] %*% t(fit$design[,add_vars])
rownames(cmc_residuals_matrix) = cmc_norm_matrix$ensembl_id
cmc_residuals_matrix = rownames_to_column(cmc_residuals_matrix, var = "ensembl_id")

# check residuals 
hist(unlist(cmc_residuals_matrix[-1]), breaks = 100)
# qqnorm(unlist(cmc_residuals_matrix[-1]), pch = 1, frame = FALSE)
# qqline(unlist(cmc_residuals_matrix[-1]), col = "steelblue", lwd = 2)

# get cell type markers ----
# load markers 
markers = read.csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'))[,-1]
colnames(markers) = colnames(markers) %>% make.names() %>% tolower()
markers = markers %>% dplyr::rename(gene_entrez = entrez.gene.id)
# update ensembl ids mapped to entrez genes 
markers = left_join(markers %>% distinct(), ensembl, by = "gene_entrez") %>%
  dplyr::select(-c(gene, ensembl.gene.id)) %>%
  dplyr::select(c(gene_symbol, gene_ensembl, gene_entrez, everything()))
#markers = markers %>% filter(used.in.mgp == "TRUE")

# format as list
cell_types = markers %>% filter(!is.na(subclass)) %>% pull(subclass) %>% unique()
marker_list = lapply(cell_types, function(cell_type){
  return(markers %>% filter(subclass == cell_type, gene_ensembl %in% cmc_residuals_matrix$ensembl_id,
  ) %>% pull(gene_ensembl))
})
names(marker_list) = cell_types

# MGP deconvolution ----
# Run marker gene profile (MGP) analysis on normalized data using subclass markers 
estimations_norm =  mgpEstimate(
  exprData = cmc_norm_matrix,
  genes = marker_list,
  geneColName = 'ensembl_id',
  outlierSampleRemove = FALSE, 
  geneTransform = NULL, 
  groups = NULL, 
  seekConsensus = FALSE, 
  removeMinority = FALSE)
# get as data frame
estimations_norm_df = as.data.frame(estimations_norm$estimates) %>%
  rownames_to_column(var = "SampleID")

# Run marker gene profile (MGP) analysis on residualized data using subclass markers 
estimations_resid =  mgpEstimate(
  exprData = cmc_residuals_matrix,
  genes = marker_list,
  geneColName = 'ensembl_id',
  outlierSampleRemove = FALSE, 
  geneTransform = NULL, 
  groups = NULL, 
  seekConsensus = FALSE, 
  removeMinority = FALSE)
# get as data frame
estimations_resid_df = as.data.frame(estimations_resid$estimates) %>%
  rownames_to_column(var = "SampleID")

# both are highly correlated 
cor(unlist(estimations_norm_df[,-1]), unlist(estimations_resid_df[,-1]))

# save non-residualized mgp estimates for later
tmp = estimations_norm_df %>%
  right_join(cmc_metadata %>% dplyr::select(SampleID, IndividualID, DX), ., by = "SampleID") %>%
  pivot_longer(cols = colnames(estimations_norm_df)[-1], names_to = "CT", values_to = "MGP") 
write.csv(tmp, file = "./output/cmc_deconvolution.csv")

# use residualized mgps for following plots 
mgp_resid_long  = cmc_metadata %>%
  dplyr::select(SampleID, IndividualID, DX, Reported_Gender, ageOfDeath, PMI) %>%
  left_join(., estimations_resid_df, by = "SampleID") %>%
  pivot_longer(cols = colnames(estimations_resid_df)[-1], names_to = "CT", values_to = "mgp") 
  
# MGP qc ----
# # set desired estimations obj
# estimations = estimations_resid
# 
# # collect various qc-related measures from MGP estimates object
# for(i in 1:length(cell_types)){ # switch too
#   cells_df = estimations$usedMarkerExpression[i] %>% as.data.frame()
#   masterlist = paste0(rownames(cells_df), collapse=', ')
#   num_markers = length(rownames(cells_df))
#   rm_marker_ratios = estimations$removedMarkerRatios[i]
#   if(!is.null(estimations$trimmedPCAs[[i]])){
#     percent_variance = ((summary(estimations$trimmedPCAs[[i]]))[6]) %>% as.data.frame()
#     percent_variance_PC1 = percent_variance[2,1]
#   } else{
#     percent_variance_PC1 = NA
#   } 
#   if(i==1){
#     master_df = data.frame("markers_used" = masterlist, "removed_marker_ratios" = rm_marker_ratios,
#                            "percent_variance_PC1" = percent_variance_PC1,"num_markers" = num_markers)  
#   } else{
#     df = data.frame("markers_used" = masterlist, "removed_marker_ratios" = rm_marker_ratios,
#                     "percent_variance_PC1" = percent_variance_PC1, "num_markers" = num_markers)
#     master_df = rbind(master_df, df)
#   }
# }
# qc_metrics = rownames_to_column(master_df, var = "celltype") 
# qc_metrics = qc_metrics[complete.cases(qc_metrics),]
# 
# # plot number of markers used 
# qc_metrics %>%  ggplot(aes(x = celltype, y = num_markers)) +
#   geom_bar(stat = "identity", fill = "#e0abf5") +
#   geom_text(aes(label = num_markers), hjust = 0) +
#   geom_hline(yintercept = 4) + 
#   labs(title = "Number of Markers Per Celltype", 
#        x="Cell Type", y = "Markers Used") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
#   coord_flip()
# #ggsave(filename = "./figures/density_deconvolution/n_markers.jpeg", width = 9, height = 5.5)
# 
# # plot variance explained by PC1 
# qc_metrics %>%  ggplot(aes(x = celltype, y = percent_variance_PC1))+
#   geom_bar(stat = "identity", fill = ifelse(qc_metrics$percent_variance_PC1 > 0.35, "#AFEEEE", "#808080")) +
#   geom_text(aes(label = round(percent_variance_PC1,2)), hjust = 0) +
#   geom_hline(yintercept = 0.35) +
#   labs(title = "Percent Variance Explained by Each MGP",
#        x="MGPs", y = "Percent Variance Explained by PC 1")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90)) + coord_flip()
# #ggsave(filename = "./figures/density_deconvolution/pct_markers.jpeg", width = 9, height = 5.5)

# FISH cell densities ----
# load densities (cells per 0.1 mm^2)
cell_dens = read.csv(file = "./data/lcm-seq/cell_densities_full.csv") %>% 
  dplyr::rename(HU = Subject, PVALB = PV, DX = Subject.Group) %>%
  mutate(DX = factor(DX, levels = c("Control","MDD","Bipolar","SCHIZ"), labels = c("CTRL","MDD","BD","SCZ"))) %>%
  mutate(HU = as.character(HU)) 
# tidy   
cell_dens_long = metadata %>%
  dplyr::select(HU, DX, Age, Sex, PMI) %>%
  right_join(., cell_dens %>% dplyr::select(-DX), by = "HU") %>%
  dplyr::rename(PyrL2n3 = PYR_23, PyrL5n6 = PYR_56) %>%
  pivot_longer(cols = c(PyrL2n3, PyrL5n6, PVALB, VIP, SST), names_to = "CT", values_to = "Density") %>%
  # make sure there are no duplicates 
  distinct(HU, DX, CT, .keep_all = TRUE) %>%
  # log transformation for computing residuals 
  mutate(Log_density = log(Density + 1)) %>%
  # drop NAs and outliers +/- 1.5 IQR   
  group_by(CT) %>%
  filter(!Log_density %in% boxplot.stats(Log_density)$out) %>% 
  ungroup() %>%
  drop_na()

# compute density residuals ----
# full model 
fit = lm(scale(Log_density) ~ DX + Sex + scale(Age) + scale(PMI), data = cell_dens_long)

# residualize and add DX back 
design = model.matrix(fit)
add_vars = grep("DX", colnames(design), value = T)
resid = residuals(fit) + fit$coefficients[add_vars] %*% t(design[,add_vars])
cell_dens_long$res = as.vector(resid)
# save for later 
#write.csv(cell_dens_long, file = "./output/FISH_densities.csv")

# # check residuals (CT must be accounted for)
# plot(fit)
# shapiro.test(residuals(fit))
# hist(residuals(fit), breaks = 50)
# qqnorm(cell_dens_long$res, pch = 1, frame = FALSE)
# qqline(cell_dens_long$res, col = "steelblue", lwd = 2)
# 
# # box-cox variance-stabilizing transformation
# boxcox(fit, plotit = TRUE, lambda = seq(-0.25, 0.25, by = 0.05))
# fit = lm((((dens ^ -0.05) - 1) / -0.05)~ DX + Reported_Gender + scale(ageOfDeath) + scale(PMI), data = cell_dens_long)
# 
# # remove outliers according to cook's distance
# cd = cooks.distance(fit)
# sum(cd > 2 / length(cd))
# cd_subset = cd < 2 / length(cd)
# cell_dens_long = cell_dens_long[cd_subset,]
# fit_2 = lm(dens ~ DX + Reported_Gender + ageOfDeath + PMI, data = cell_dens_long)
# 
# # check model assumptions
# plot(fit)
# shapiro.test(residuals(fit))

# visualize mpg and cell density across disorders ----

# plot cell type densities across disorders 
p1 = cell_dens_long %>% 
  filter(CT %in% c("PyrL2n3","PyrL5n6","PVALB","VIP","SST")) %>%
  mutate(CT = factor(CT, levels = c("PyrL2n3","PyrL5n6","PVALB","VIP","SST"))) %>%
ggplot(aes(x = DX, y = res, color = DX)) +
  geom_violin(aes(fill = DX), alpha = 0.4) +
  geom_quasirandom(na.rm = T, shape = 16, alpha = 0.9) + 
  stat_summary(fun = mean, geom = "point", aes(fill = DX), size = 4, shape = 23, col = "black") +
  labs(y = 'Cellular density residuals\n(AU, ACC)', x = "", title = "", fill = "Diagnostic group", color = "Diagnostic group") +  
  scale_fill_manual(values = cols_1) +
  scale_color_manual(values = cols_1) + 
  theme_classic2() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(),
        strip.text = element_text(face = "bold"),
        plot.margin = margin(0,0,-1,0,"cm")) +
  stat_compare_means(comparisons = list(c('CTRL','MDD'), c('CTRL','BD'), c('CTRL','SCZ')), method = "wilcox.test", size = 3) +
  facet_wrap(~CT, ncol = 5, scales = "free_y")
p1

# plot cell type proportions across disorders 
p2 = mgp_resid_long  %>%
  filter(CT %in% c("IT","L6b","PVALB","VIP","SST")) %>%
  mutate(CT = factor(CT, levels = c("IT","L6b","PVALB","VIP","SST"))) %>%
ggplot(aes(x = DX, y = mgp, color = DX)) +
  geom_violin(aes(fill = DX), alpha = 0.4) +
  geom_quasirandom(na.rm = T, shape = 16, alpha = 0.9) + 
  stat_summary(fun = mean, geom = "point", aes(fill = DX), size = 4, shape = 23, col = "black") +
  scale_x_discrete(limits = c("CTRL", "MDD", "BD", "SCZ")) +
  labs(y = 'rCTP residuals\n(AU, PFC)', x = "", color = "Diagnosis", fill = "Diagnosis", title = "") +  
  scale_fill_manual(values = cols_1) + 
  scale_color_manual(values = cols_1) +
  theme_classic2() +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(),
        strip.text = element_text(face = "bold"),
        plot.margin = margin(-1,0,-1,0,"cm")) +
  stat_compare_means(comparisons = list(c('CTRL','BD'), c('CTRL','SCZ')), method = "wilcox.test", size = 3) +
  facet_wrap(~CT, ncol = 7, scales = "free_y")
p2

# linear models ---- 
# linear model for cell densities 
cell_dens_lm = cell_dens_long %>%
  mutate(CT = factor(CT, levels = c("PyrL2n3","PyrL5n6","PVALB","VIP","SST"), labels = c("PyrL2n3/IT","PyrL5n6/L6b","PVALB","VIP","SST"))) %>%
  group_by(CT) %>%
  do(tidy(lm(scale(Log_density) ~ DX + Sex + scale(Age) + scale(PMI), .), conf.int = T)) %>%  
  ungroup() %>%
  mutate(padj = p.adjust(`p.value`, method = 'BH')) %>%
  mutate(assay = "dens") %>%
  mutate(term = recode(term, 
                       `(Intercept)` = "Intercept", 
                       `DXMDD` = "MDD",
                       `DXBD` = "BD",
                       `DXSCZ` = "SCZ",
                       `scale(ageOfDeath)` = "Age",
                       `Reported_GenderFemale` = "Reported_Gender"))

# linear model for cell proportion estimates 
cell_mgp_lm = mgp_resid_long  %>%
  filter(CT %in% c("IT","L6b","PVALB","VIP","SST")) %>%
  mutate(CT = factor(CT, levels = c("IT","L6b","PVALB","VIP","SST"), labels = c("PyrL2n3/IT","PyrL5n6/L6b","PVALB","VIP","SST"))) %>%
  group_by(CT) %>%
  do(tidy(lm(scale(mgp) ~ DX + Reported_Gender + scale(ageOfDeath) + scale(PMI), data = .), conf.int = T)) %>%
  ungroup() %>%
  mutate(padj = p.adjust(`p.value`, method = 'BH')) %>%
  mutate(assay = "mgp") %>%
  mutate(term = recode(term, 
                       `(Intercept)` = "Intercept", 
                       `DXMDD` = "MDD",
                       `DXBD` = "BD",
                       `DXSCZ` = "SCZ",
                       `scale(ageOfDeath)` = "Age",
                       `Reported_GenderFemale` = "Reported_Gender"))
# combine
lm_df = rbind(cell_dens_lm, cell_mgp_lm) %>%
  mutate(sig = case_when(padj < 0.05 ~ TRUE, TRUE ~ FALSE))

# plot standard beta coefficients 
p3 = lm_df %>%
  filter(term %in% c("MDD", "BD", "SCZ")) %>% 
  mutate(term = factor(term, levels = c("MDD", "BD", "SCZ"))) %>% 
  #mutate(CT = factor(CT, levels = c("PYR/GLU","VIP","SST","PVALB"))) %>%
ggplot(aes(x = term, y = estimate, fill = assay)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = "dodge") +
  geom_text(aes(label = ifelse(sig, "*", ""), y = conf.low - 0.3, group = assay), color = "red",
            position = position_dodge(width = .9), size = 20 / .pt) +
  scale_y_continuous(limits = c(-1.8,1.5)) +
  labs(y = 'Std. Beta coeff.', x = 'Diagnostic group', title = "") + 
  guides(fill = guide_legend(title = "Data type")) +
  scale_fill_d3() +
  theme_classic2() +
  theme(plot.title = element_blank(),
        axis.line = element_line(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text = element_text(colour = "black"),
        strip.text = element_text(face = "bold"),
        plot.margin = margin(-1,0,0,0,"cm")) +
  facet_wrap(~CT, drop = T, scale = "free", nrow = 1)

# combine figure
p1 + p2 + p3 + 
  plot_layout(ncol = 1)
ggsave(filename = "./figures/draft/dx_mgp_dens.svg", width = 10, height = 9)

