# differential-expression-analysis.R 
# load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(edgeR)
  library(DESeq2)
  library(FactoMineR)
  library(umap)
  library(Rtsne)
  library(transcripTools)
  library(variancePartition)
  library(ComplexHeatmap)
  library(ggpubr)
  library(ggVennDiagram)
  library(patchwork)
  library(pheatmap)
  library(svglite)
})
cols_1 = c("PyrL2n3" = "#D62728FF", "PyrL5n6" = "#FF7F0EFF", "PVALB" = "#2CA02CFF", "SST" = "#9467BDFF", "VIP" = "#1F77B4FF")

# load data ----
load(file = "./data/lcm-seq/data.Rdata")
# add star qc 
metadata = metadata %>% left_join(., star_qc, by = 'ID')
rownames(metadata) = metadata$ID
# clean-up medications 
metadata = metadata %>%
  mutate(Tobacco = factor(ifelse(grepl("Y", Tob.ATOD),"Y","N"))) %>%
  mutate(Antidepressants = factor(ifelse(grepl("D", Meds.ATOD),"Y","N"))) %>%
  mutate(Benzos_Anticonvulsants = factor(ifelse(grepl(paste0(c("B","C"), collapse = "|"), Meds.ATOD),"Y","N"))) %>%
  mutate(Antipsychotics = factor(ifelse(grepl("P", Meds.ATOD),"Y","N"))) %>%
  mutate(Lithium = factor(ifelse(grepl("L", Meds.ATOD),"Y","N"))) 

# # remove outliers detected during qc steps
# metadata = metadata %>% filter(!ID %in% outliers)
# match and check
counts_matrix = counts_matrix[, match(rownames(metadata), colnames(counts_matrix))]
all(rownames(metadata) == colnames(counts_matrix))

# calculate scaling factors to convert raw library sizes into effective library sizes
dge0 = DGEList(counts_matrix)
dge0 = calcNormFactors(dge0, method = "TMM")

# stratified limma-voom ----
# run for each cell type separately 
lcm_lst_limma = by(metadata, list(metadata$CT), simplify = FALSE, function(x){
  # keep genes expressed in at least 15 samples (~80% of the smallest number of replicates (n = 19))  
  keep = rowSums(cpm(dge0[,x$ID]) > 1) >= 15
  dge = dge0[keep, x$ID]
  # specify the model to be fitted (no intercept, means model)
  # we do this before using voom since voom uses variances of the model residuals (observed - fitted)
  design = model.matrix(~0 + DX + Sex + Pool + MOD + scale(Age) + scale(PMI) + scale(RIN) + scale(Intergenic_rate), x) 
  colnames(design) = gsub("\\(|\\)","", colnames(design))
  
  contrs = makeContrasts(DXMDD = DXMDD-DXCTRL,
                         DXBD = DXBD-DXCTRL,
                         DXSCZ = DXSCZ-DXCTRL, levels = colnames(design))
  # voom to obtain log-counts-per-million values and observational level weights
  vm = voomWithQualityWeights(dge, design, plot = FALSE)
  # fit a linear model using weighted least squares for each gene
  vfit = lmFit(vm, design)
  # estimate contrast for each gene
  vfit = contrasts.fit(vfit, contrs)
  # empirical Bayes smoothing of standard errors 
  # shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error
  # `robust = TRUE` protects the procedure against hyper-variable or hypo-variable genes, especially when analyzing RNA-seq data
  efit = eBayes(vfit, robust = TRUE)
  
  # plot residual standard deviation versus average log expression for a fitted linear model
  # plotSA(efit, main = as.character(unique(x$CT)))
  
  # get results per disorder contrast
  contrasts = colnames(contrs)
  tmp_lst = lapply(contrasts, function(y){
    res = topTable(efit, coef = y, sort = "none", n = Inf, adjust.method = "BH") %>% rownames_to_column(var = "gene_ensembl")
    res$CT = unique(as.character(x$CT))
    res$DX = y
    return(res)
  })
  res = do.call(rbind, tmp_lst)
})
# combine and clean 
lcm_df_limma = do.call(rbind, lcm_lst_limma) %>%
  relocate(DX, .after = "CT") %>% 
  relocate(gene_ensembl, .before = "CT") %>%
  mutate(DX = substring(DX, 3, nchar(DX)),
         DE = adj.P.Val < 0.2,
         CT = factor(CT, levels =  c("PVALB","SST","VIP","PyrL2n3","PyrL5n6")),
         DX = factor(DX, levels = c("MDD","BD","SCZ"), labels = c("MDD","BD","SCZ"))) %>%
  # add gene symbols 
  left_join(., ensembl, by = "gene_ensembl") %>%
  dplyr::relocate(gene_symbol, .after = gene_ensembl)
#write_rds(lcm_df_limma, file = "./output/lcm_df_limma.rds")

# plots ----
# keep genes that are expressed in at least one combination of CT and DX
filt_lst = by(metadata, list(metadata$CT, metadata$DX), function(x){
  rowSums(cpm(counts_matrix[,x$ID]) > 1) >= 15
})
keep = colSums(do.call(rbind, filt_lst)) > 0 
counts_matrix_filt = counts_matrix[keep,]
lcpm_matrix_filt = lcpm_matrix[keep,]
# variance stabilizing transformation 
dds = DESeqDataSetFromMatrix(counts_matrix_filt, metadata, design = ~ CT + DX)
vst = vst(dds, blind = FALSE)
vst_matrix = vst@assays@data@listData[[1]]

## sample dimensionality reduction ----
## sample PCA
tmp = mostVar(vst_matrix, n = 2000)
pca_data = PCA(t(tmp), ncp = 4, graph = FALSE)
percent_var = pca_data$eig[,2]
pca_data = pca_data$ind$coord %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  left_join(., metadata, by = "ID")

p1 = ggplot(pca_data, aes(x = Dim.1, y = Dim.4, color = CT, shape = DX)) +
  geom_point(size = 3) +
  scale_color_manual(values = cols_1) +
  xlab(paste0("PC1: ", round(percent_var[1],0), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[4],0), "% variance")) +
  coord_fixed() +
  guides(color = guide_legend(label.position = "left"),
         shape = guide_legend(label.position = "left")) +
  labs(color = "", shape = "") +
  theme_classic2() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "left", legend.margin = margin())

## sample UMAP
umap_setting = umap.defaults
umap_setting$spread = 8
umap_setting$n_neighbors = 30
umap_setting$min_dist = 5

tmp = mostVar(cpm(counts_matrix, log = F), n = 2000)
umap_out = umap(t(tmp), umap_setting)
umap_data = umap_out$layout %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  left_join(., metadata, by = "ID")

p1 = ggplot(umap_data, aes(x = V1, y = V2, color = CT, shape = DX)) +
  geom_point(size = 3) +
  scale_color_manual(values = cols_1) +
  labs(x = "UMAP 1", y = "UMAP 2") + 
  guides(color = guide_legend(label.position = "left"),
         shape = guide_legend(label.position = "left")) +
  labs(color = "", shape = "") +
  theme_classic2() +
  theme(legend.position = "none", legend.margin = margin(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, face = "bold"))
p1

kmeans_out = kmeans(umap_out$layout, centers = 5)
kmeans_out = kmeans_out$cluster
true_clust_unique = unique(metadata$CT)
mapping = match(metadata$CT, true_clust_unique)
ari = mclust::adjustedRandIndex(kmeans_out, mapping)
((ari + 1) / 2) * 100

# cell type marker expression 
markers = ensembl %>% filter(gene_symbol %in% c("PVALB","SST","VIP","CUX2","FEZF2","STMN2")) %>% pull(gene_ensembl)
p2 = cpm(counts_matrix, log = TRUE) %>%
  as.data.frame() %>%
  filter(row.names(.) %in% markers) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  merge(metadata %>% dplyr::select(ID, CT, DX), ., by = "ID") %>%
  pivot_longer(cols = -c(ID, CT, DX), values_to = "exp", names_to = "gene_ensembl") %>%
  left_join(., ensembl %>% dplyr::select(gene_ensembl, gene_symbol), by = "gene_ensembl") %>%
  mutate(gene_symbol = factor(gene_symbol, levels = c("PVALB","SST","VIP","CUX2","FEZF2","STMN2"))) %>%
  mutate(exp = na_if(exp, 0)) %>%
ggplot(., aes(x = CT, y = exp, fill = CT)) + 
  geom_violin(width = 0.75, alpha = 1) +
  #scale_y_continuous(trans = "log2") + 
  scale_fill_manual(values = cols_1) +
  guides(fill = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2)) +
  labs(y = "Expression (logCPM)", x = "Cell type") + 
  theme_classic2() + 
  theme(legend.position = "none", 
        legend.box = "vertical",
        legend.margin = margin(),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 16, face = "bold", color = "black")) +
  stat_compare_means(method = "wilcox.test",                        label = "p.signif",
                     label.y = c(30, 35), # Adjust the vertical positioning of the labels
                     label.x = c(1.5, 2.5), # Adjust the horizontal positioning of the labels
                     vjust = 0.2) +
  facet_wrap(~gene_symbol, drop = TRUE, scales = "free_y", nrow = 2) 
p2

svglite(filename = "./figures/draft/tmp.svg", height = 5.5, width = 10)
p1 + p2 + plot_layout(widths = c(1,1.2))
dev.off()

## number of DEGs ----
lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds") 
lcm_df_limma %>%  
  distinct(CT, DX, gene_ensembl, .keep_all = TRUE) %>%
  mutate(Dir = factor(sign(logFC))) %>%
  mutate(Dir = recode_factor(Dir, `-1` = "Down", `1` = "Up")) %>%
  group_by(CT, DX, Dir) %>% 
  dplyr::summarize("pt05" = sum(adj.P.Val < 0.05),
            "pt10" = sum(adj.P.Val < 0.10),
            "pt20" = sum(adj.P.Val < 0.20)) %>%
  pivot_longer(cols = c(pt05, pt10, pt20),
               values_to = "Freq",
               names_to = "FDR_thresh") %>%
  mutate(Freq = if_else(Dir == "Down", -Freq, Freq)) %>%
ggplot(., aes(x = DX, y = Freq, fill = Dir)) +
  geom_bar(data = . %>% filter(FDR_thresh == "pt05"), stat = "identity", width = 0.8, alpha = 1) +
  geom_bar(data = . %>% filter(FDR_thresh == "pt10"), stat = "identity", width = 0.8, alpha = 0.7) +
  geom_bar(data = . %>% filter(FDR_thresh == "pt20"), stat = "identity", width = 0.8, alpha = 0.5) +
  geom_tile(aes(y = NA_integer_, alpha = factor(FDR_thresh))) + 
  scale_fill_manual(values = c("#1F77B4FF","#D62728FF")) +
  scale_alpha_manual("FDR Threshold", values = c(1,0.6,0.3), labels = c("0.05","0.10","0.20")) + 
  scale_y_continuous(limits = c(-90, 90), breaks = seq(-90, 90, 20), labels = abs(seq(-90, 90, 20))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "top")) +
  labs(y = "Number of DE genes", x = "", fill = "Diagnostic group") +
  coord_flip() +
  facet_wrap(~CT, ncol = 1) +
  theme_classic2() +
  theme(strip.text = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.y = element_text(face = "bold"), 
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom")

#ggsave(filename = "./figures/draft/de_n.svg", height = 8, width = 6, dpi = 600)

lcm_df_limma %>% distinct(DX, CT, gene_ensembl, .keep_all = T) %>% filter(adj.P.Val < 0.20) %>% pull(gene_ensembl) %>% length()

# ggplot() +
#   geom_point(aes(x = DX, y = logFC), data = lcm_df_limma, color = "grey", alpha = 0.05) +
#   ggbeeswarm::geom_quasirandom(aes(x = DX, y = logFC), data = lcm_df_limma %>% filter(DE == TRUE, logFC > 0), color = "#D62728FF", width = 0.2, size = 2.2, shape = 16) +
#   ggbeeswarm::geom_quasirandom(aes(x = DX, y = logFC), data = lcm_df_limma %>% filter(DE == TRUE, logFC < 0), color = "#1F77B4FF", width = 0.2, size = 2.2, shape = 16) + 
#   geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") + 
#   scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 2), labels = abs(seq(-6, 6, 2))) +
#   labs(y = "Log fold-change", x = "Diagnostic group", fill = "Logfold-change") +
#   coord_flip() +
#   facet_wrap(~CT, ncol = 1) +
#   theme_classic2() + 
#     theme(strip.text = element_text(face = "bold", size = 16),
#           axis.text = element_text(size = 12, color = "black"),
#           axis.text.y = element_text(face = "bold"), 
#           axis.title = element_text(size = 16, face = "bold"))
# #ggsave(filename = "./figures/draft/de_logfc.svg", height = 8, width = 6, dpi = 600)

# compare to bulk tables ----

lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds") 

# bulk microarray (MDD, BD and SCZ)
bulk_microarray = read.csv(file = "./data/gandal/microarray_meta_analyses.csv", row.names = 1) %>% 
  dplyr::rename(gene_ensembl = ensembl_gene_id) %>%
  mutate(logFDR_microarray = sign(SCZ.beta_log2FC)*-log10(SCZ.FDR), DX = "SCZ") %>%
  dplyr::select(DX, gene_ensembl, logFDR_microarray)
bulk_microarray = rbind(bulk_microarray, 
                       read.csv(file = "./data/gandal/microarray_meta_analyses.csv", row.names = 1) %>% 
                         dplyr::rename(gene_ensembl = ensembl_gene_id) %>%
                         mutate(logFDR_microarray = sign(BD.beta_log2FC)*-log10(BD.FDR), DX = "BD") %>%
                         dplyr::select(DX, gene_ensembl, logFDR_microarray))
bulk_microarray = rbind(bulk_microarray, 
                       read.csv(file = "./data/gandal/microarray_meta_analyses.csv", row.names = 1) %>% 
                         dplyr::rename(gene_ensembl = ensembl_gene_id) %>%
                         mutate(logFDR_microarray = sign(MDD.beta_log2FC)*-log10(MDD.FDR), DX = "MDD") %>%
                         dplyr::select(DX, gene_ensembl, logFDR_microarray))

bulk_rnaseq = fread(file = "./data/gandal/RNAseq_CMC.csv", data.table = F) %>%
  dplyr::rename(gene_ensembl = V1) %>%
  mutate(logFDR_rnaseq = sign(SCZ.logFC)*-log10(SCZ.adj.P.Val), DX = "SCZ") %>%
  dplyr::select(DX, gene_ensembl, logFDR_rnaseq)
bulk_rnaseq = rbind(bulk_rnaseq, 
                    fread(file = "./data/gandal/RNAseq_CMC.csv", data.table = F) %>%
                      dplyr::rename(gene_ensembl = V1) %>%
                      mutate(logFDR_rnaseq = sign(BD.logFC)*-log10(BD.adj.P.Val), DX = "BD") %>%
                      dplyr::select(DX, gene_ensembl, logFDR_rnaseq))

lcm_rnaseq = lcm_df_limma %>%
  distinct(CT, DX, gene_ensembl, .keep_all = TRUE) %>%
  mutate(logFDR_lcmseq = sign(logFC)*-log10(adj.P.Val)) %>%
  dplyr::select(CT, DX, gene_ensembl, logFDR_lcmseq) 

bulk_overlaps = lcm_rnaseq %>%
  filter(abs(logFDR_lcmseq) > -log10(0.2)) %>%
  left_join(bulk_microarray %>% filter(abs(logFDR_microarray) > -log10(0.2)), by = c("gene_ensembl", "DX")) %>%
  left_join(bulk_rnaseq %>% filter(abs(logFDR_rnaseq) > -log10(0.2)), by = c("gene_ensembl", "DX")) %>%
  mutate(matching_status_bulk_microarray = case_when(
    sign(logFDR_lcmseq) == sign(logFDR_microarray) ~ "Concordant",
    sign(logFDR_lcmseq) != sign(logFDR_microarray) ~ "Discordant",
    TRUE ~ NA_character_
  ),
  matching_status_bulk_rnaseq = case_when(
    sign(logFDR_lcmseq) == sign(logFDR_rnaseq) ~ "Concordant",
    sign(logFDR_lcmseq) != sign(logFDR_rnaseq) ~ "Discordant",
    TRUE ~ NA_character_
  )) %>%
  pivot_longer(cols = c(matching_status_bulk_microarray, matching_status_bulk_rnaseq),
               names_to = "dataset",
               values_to = "status") %>%
  mutate(dataset = case_when(
    dataset == "matching_status_bulk_microarray" ~ "bulk_microarray",
    dataset == "matching_status_bulk_rnaseq" ~ "bulk_rnaseq",
    TRUE ~ NA_character_
  )) %>%
  left_join(., lcm_df_limma, by = c("CT","DX","gene_ensembl"))
  
bulk_overlaps_n = bulk_overlaps %>%
  group_by(CT, DX, dataset, status) %>%
  summarize(num_overlaps = n_distinct(gene_ensembl),
            overlapping_genes = paste(unique(gene_ensembl), collapse = ", ")) %>%
  drop_na() 

bulk_overlaps_n_sum = bulk_overlaps_n %>%
  mutate(overlapping_genes = strsplit(overlapping_genes, ", ")) %>%
  group_by(status) %>%
  summarize(num_overlaps = length(unique(unlist(overlapping_genes))),
            overlapping_genes = paste(unique(unlist(overlapping_genes)), collapse = ", ")) 

phyper(
  q = 43 - 1,
  m = lcm_df_limma %>% filter(adj.P.Val <= 0.2) %>% pull(gene_ensembl) %>% unique() %>% length(),
  n = bulk_microarray %>% filter(abs(logFDR_microarray) >= -log10(0.2)) %>% pull(gene_ensembl) %>% unique() %>% length(),
  k = 43, 
  lower.tail = FALSE
)

# compare overlaps tabels ----

## between DX
N = as.numeric(length(unique(lcm_df_limma$gene_ensembl)))
K = as.numeric(length(unique(lcm_df_limma %>% filter(DE == TRUE) %>% pull(gene_ensembl))))
n = K
ls = list(
  lcm_df_limma %>% filter(DX == "MDD", DE == TRUE) %>% pull(gene_ensembl) %>% unique(),
  lcm_df_limma %>% filter(DX == "BD", DE == TRUE) %>% pull(gene_ensembl) %>% unique(),
  lcm_df_limma %>% filter(DX == "SCZ", DE == TRUE) %>% pull(gene_ensembl) %>% unique()
)
tb = table(do.call(c, ls))
k = as.numeric(length(names(tb[tb >= 2])))

phyper(k-1, K, N-K, n, lower.tail = FALSE)

## between CT
N = as.numeric(length(unique(lcm_df_limma$gene_ensembl)))
K = as.numeric(length(unique(lcm_df_limma %>% filter(DE == TRUE) %>% pull(gene_ensembl))))
n = K
ls = list(
  lcm_df_limma %>% filter(CT == "PVALB", DE == TRUE) %>% pull(gene_ensembl) %>% unique(),
  lcm_df_limma %>% filter(CT == "SST", DE == TRUE) %>% pull(gene_ensembl) %>% unique(),
  lcm_df_limma %>% filter(CT == "VIP", DE == TRUE) %>% pull(gene_ensembl) %>% unique(),
  lcm_df_limma %>% filter(CT == "PyrL2n3", DE == TRUE) %>% pull(gene_ensembl) %>% unique(),
  lcm_df_limma %>% filter(CT == "PyrL5n6", DE == TRUE) %>% pull(gene_ensembl) %>% unique()
)
tb = table(do.call(c, ls))
k = as.numeric(length(names(tb[tb >= 2])))

phyper(k-1, K, N-K, n, lower.tail = FALSE)

## venn diagrams ----
# number of DE genes across DX
de_dx = list(
  MDD = lcm_df_limma %>% filter(DE == TRUE, logFC > 0, DX == "MDD") %>% pull(gene_ensembl) %>% unique(),
  BD = lcm_df_limma %>% filter(DE == TRUE, logFC > 0, DX == "BD") %>% pull(gene_ensembl) %>% unique(),
  SCZ = lcm_df_limma %>% filter(DE == TRUE, logFC > 0, DX == "SCZ") %>% pull(gene_ensembl) %>% unique()
)
p1 = ggVennDiagram(de_dx, label = "count", label_alpha = 0, label_size = 10) + 
  scale_fill_gradient(low = "white", high = "#D62728FF") +
  scale_color_manual(values = rep("#D62728FF",4)) +
  theme(legend.position = "none")

# number of DE genes across DX
de_dx = list(
  MDD = lcm_df_limma %>% filter(DE == TRUE, logFC < 0, DX == "MDD") %>% pull(gene_ensembl),
  BD = lcm_df_limma %>% filter(DE == TRUE, logFC < 0, DX == "BD") %>% pull(gene_ensembl),
  SCZ = lcm_df_limma %>% filter(DE == TRUE, logFC < 0, DX == "SCZ") %>% pull(gene_ensembl)
)
p2 = ggVennDiagram(de_dx, label = "count", label_alpha = 0, label_size = 10) + 
  scale_fill_gradient(low = "white", high = "#1F77B4FF") +
  scale_color_manual(values = rep("#1F77B4FF",4)) +
  theme(legend.position = "none")

# number of upregulated DE genes across CT
de_ct = list(
  PVALB = lcm_df_limma %>% filter(DE == TRUE, logFC > 0, CT == "PVALB") %>% pull(gene_ensembl),
  SST = lcm_df_limma %>% filter(DE == TRUE, logFC > 0, CT == "SST") %>% pull(gene_ensembl),
  VIP = lcm_df_limma %>% filter(DE == TRUE, logFC > 0, CT == "VIP") %>% pull(gene_ensembl),
  PYR = lcm_df_limma %>% filter(DE == TRUE, logFC > 0, CT %in% c("PyrL2n3","PyrL5n6")) %>% pull(gene_ensembl)
)
p3 = ggVennDiagram(de_ct, label = "count", label_alpha = 0, label_size = 10) + 
  scale_fill_gradient(low = "white", high = "#D62728FF") +
  scale_color_manual(values = rep("#D62728FF",4)) +
  theme(legend.position = "none")

# number of DE genes across CT
de_ct = list(
  PVALB = lcm_df_limma %>% filter(DE == TRUE, logFC < 0, CT == "PVALB") %>% pull(gene_ensembl),
  SST = lcm_df_limma %>% filter(DE == TRUE, logFC < 0, CT == "SST") %>% pull(gene_ensembl),
  VIP = lcm_df_limma %>% filter(DE == TRUE, logFC < 0, CT == "VIP") %>% pull(gene_ensembl),
  PYR = lcm_df_limma %>% filter(DE == TRUE, logFC < 0, CT %in% c("PyrL2n3","PyrL5n6")) %>% pull(gene_ensembl)
)
p4 = ggVennDiagram(de_ct, label = "count", label_alpha = 0, label_size = 10) + 
  scale_fill_gradient(low = "white", high = "#1F77B4FF") +
  scale_color_manual(values = rep("#1F77B4FF",4)) +
  theme(legend.position = "none")

p1/p2/p3/p4
ggsave(filename = "./figures/draft/venns.svg", width = 5.5, height = 18)

## DE heatmaps ----

lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds") 

dx = "MDD"

# Define separate scaling factor for the number of genes included across all CT groups
scaling_factor = 0.20 # You can adjust this value to scale the number of included top genes
scaling_factor_highest = 0.14  # You can adjust this value for the CT group with the highest proportion

# MDD
# 0.22
# 0.14
# BD
# 0.30
# 0.14
# SCZ
# 0.10
# 0.06

tmp = lcm_df_limma %>%
  filter(
    #!CT %in% c("PyrL2n3", "PyrL5n6"), 
    DX == dx, !gene_symbol == "", gene_biotype == "protein_coding") %>%
  mutate(signed_logFDR = sign(logFC) * -log10(adj.P.Val))

# Calculate the number of DE genes for each CT and DX combination
num_genes_DE_per_CT_DX = tmp %>%
  group_by(CT, DX) %>%
  filter(DE == TRUE) %>%
  summarise(num_DE_genes = n()) %>%
  ungroup()

# Merge the proportions with the main dataframe
tmp = tmp %>%
  left_join(num_genes_DE_per_CT_DX, by = c("CT", "DX"))

# Calculate the proportion of genes to be selected for each CT group based on the DE gene counts
proportions_per_CT = tmp %>%
  distinct(CT, num_DE_genes) %>%
  mutate(proportion = num_DE_genes / max(num_genes_DE_per_CT_DX$num_DE_genes)) %>%
  select(CT, proportion)

# Identify the CT group with the highest proportion of DE genes
highest_proportion_CT = proportions_per_CT %>%
  drop_na() %>%
  filter(proportion == max(proportion)) %>%
  pull(CT) %>%
  as.character()

# Calculate the number of top genes to select for each CT group based on the proportions and scaling factors
num_top_genes_per_CT = proportions_per_CT %>%
  mutate(num_top_genes = ifelse(CT == highest_proportion_CT, 
                                round(proportion * sum(tmp$DE) * scaling_factor_highest),
                                round(proportion * sum(tmp$DE) * scaling_factor))) %>%
  select(CT, num_top_genes) %>%
  drop_na()

# Function to slice data for each CT group
slice_data_for_CT = function(data) {
  num_top_genes = num_top_genes_per_CT %>%
    filter(CT == unique(data$CT)) %>%
    pull(num_top_genes)
  
  if (nrow(data) < num_top_genes) {
    return(data)
  } else {
    return(data %>%
             slice_min(order_by = signed_logFDR, n = num_top_genes) %>%
             bind_rows(data %>%
                         slice_max(order_by = signed_logFDR, n = num_top_genes)))
  }
}

# Apply the slicing function for each CT group
top_genes_full = tmp %>%
  inner_join(num_top_genes_per_CT, by = "CT") %>%
  group_split(CT) %>%
  lapply(slice_data_for_CT) %>%
  bind_rows() %>%
  arrange(match(CT, levels(lcm_df_limma$CT))) %>%
  pull(gene_symbol, CT)

top_genes = unique(top_genes_full) 

rnaseq_res = read.csv(file = "./data/gandal/RNAseq_CMC.csv", row.names = 1) %>%
  rownames_to_column(var = "gene_ensembl") %>%
  left_join(., ensembl %>% dplyr::select(gene_ensembl, gene_symbol), by = "gene_ensembl") %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  filter(gene_symbol %in% top_genes) %>%
  mutate(Bulk_RNAseq = sign(SCZ.logFC)*-log10(SCZ.adj.P.Val)) %>%
  dplyr::select(gene_symbol, Bulk_RNAseq)

microarray_res = read.csv(file = "./data/gandal/microarray_meta_analyses.csv")[-1] %>% 
  dplyr::select(hgnc_symbol, SCZ.beta_log2FC, SCZ.FDR) %>%
  dplyr::rename(gene_symbol = hgnc_symbol) %>%
  filter(gene_symbol %in% top_genes) %>%
  mutate(Bulk_microarray = sign(SCZ.beta_log2FC)*-log10(SCZ.FDR)) %>%
  dplyr::select(gene_symbol, Bulk_microarray)

mat = tmp %>% 
  filter(gene_symbol %in% top_genes) %>%
  dplyr::select(gene_symbol, CT, DX, signed_logFDR) %>%
  pivot_wider(id_cols = gene_symbol,
              names_from = CT, 
              values_from = signed_logFDR,
              values_fn = mean) %>%
  arrange(match(gene_symbol, top_genes)) %>%
  left_join(., microarray_res, by = "gene_symbol") %>%
  #left_join(., rnaseq_res, by = "gene_symbol") %>%
  column_to_rownames(var = "gene_symbol") %>%
  #na.omit() %>%
  as.matrix()

colnames(mat)[4:5] = c("L2/3 PYR","L5/6 PYR")
mat[abs(mat) < 0.05] = 0
mat[is.na(mat)] = 0

sig_mat = mat
sig_mat[abs(sig_mat) > -log10(0.2)] = "*"
sig_mat[sig_mat != "*"] = ""

my_breaks = seq(-2,2,0.15)
colors = c(rev(colorRampPalette(c("white", "#1F77B4FF"))(length(my_breaks)/2)), colorRampPalette(c("white", "#D62728FF"))(length(my_breaks)/2)) 
# colors = center.palette(mat, palette_length = 100, color1 = "#1F77B4FF", color2 = "#D62728FF")$colors
# wpos = which(colors %in% "#FFFFFF")
# colors = append(colors, rep("#FFFFFF", 9), wpos[1]) # more white
# n = 1
# colors = c(rep(colors[1], (wpos[1])*n - 5), colors, rep(colors[length(colors)], (length(colors) - wpos[2])*n)) # more red and blue

shared_genes = lcm_df_limma %>%
  distinct(CT, DX, gene_ensembl, .keep_all = TRUE) %>%
  group_by(CT, gene_symbol) %>%
  filter(sum(DE == TRUE) >= 2) %>%
  group_by(CT, gene_symbol, logFC_sign = sign(logFC)) %>%
  filter(n() >= 2) %>%
  mutate(DX_shared = paste(unique(DX[DE == TRUE]), collapse = ",")) %>%
  ungroup() %>%
  group_by(CT, gene_symbol, logFC_sign) %>%
  filter(n() >= 2) %>%
  mutate(DX_shared = paste(unique(DX[DE == TRUE]), collapse = ",")) %>%
  ungroup() %>%
  select(gene_symbol, CT, logFC_sign, DX_shared) %>%
  distinct() %>%
  filter(!gene_symbol == "") %>%
  mutate(DX_shared_N = str_count(DX_shared, ",") + 1) %>%
  left_join(., ensembl, "gene_symbol")

shared_genes %>% filter(DX_shared_N == 2) %>% pull(gene_symbol) %>% unique() %>% length()
shared_genes %>% filter(DX_shared_N == 3) %>% pull(gene_symbol) %>% unique() %>% length()
shared_genes %>% filter(DX_shared_N > 0) %>% pull(gene_symbol) %>% unique() %>% length()


# 1  2  3 
# 11 38  6 

black_genes = shared_genes %>% filter(DX_shared_N == 1) %>% pull(gene_symbol) %>% unique()
bold_genes = shared_genes %>% filter(DX_shared_N %in% c(2,3)) %>% filter(!gene_symbol %in% black_genes) %>% pull(gene_symbol) %>% unique()
bold_underline_genes = shared_genes %>% filter(DX_shared_N == 3) %>% filter(!gene_symbol %in% black_genes) %>% pull(gene_symbol) %>% unique()

make_face_names = function(mat, rc_fun, rc_names_b = NA, 
                           rc_names_i = NA) {
  f_names = rc_fun(mat)
  ids_b = rc_names_b %>% match(rc_fun(mat))
  ids_i = rc_names_i %>% match(rc_fun(mat))
  
  ids_b %>%
    walk(
      function(i)
        f_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  ids_i %>%
    walk(
      function(i)
        f_names[i] <<-
        bquote(italic(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  
  f_names
}

svglite(filename = paste("./figures/draft/degs_heatmap_",dx,".svg"), width = 3, height = 10)
pheatmap(mat,
         display_numbers = sig_mat,
         number_color = "white",
         color = colors,
         breaks = my_breaks,
         border_color = "black",
         cluster_rows = F, cluster_cols = F,
         angle_col = "45",
         cellheight = 10, cellwidth = 12,
         gaps_col = 5,
         labels_row = make_face_names(mat, 
                                      rownames, rc_names_b = bold_genes, rc_names_i = black_genes)
)
dev.off()

genes_select = shared_genes %>% 
  filter(DX_shared_N == 1) %>% 
  pull(gene_symbol) %>%
  unique()
shared_table = lcm_df_limma %>%
  filter(gene_symbol %in% genes_select) %>%
  mutate(group = paste(CT,DX,sep = "_"),
         logFDR = sign(logFC)*-log10(adj.P.Val)) %>%
  dplyr::select(gene_symbol, gene_description, group, logFDR) %>%
  pivot_wider(values_from = logFDR, names_from = group) %>%
  dplyr::select(gene_symbol, gene_description, starts_with("PVALB"), starts_with("SST"), starts_with("VIP"), 
                starts_with('PyrL2n3'), starts_with('PyrL5n6')) 

tmp = shared_table[, -c(1, 2)]
tmp[is.na(tmp)] = 0
ordering = seriation::seriate(tmp, margin = 1)
ordering = order.dendrogram(as.dendrogram(ordering[[1]]))
shared_table = shared_table[ordering, ]
  
write.csv(shared_table, file = "./output/shared_degs.csv", row.names = F)
  
  
  
  
  
  
  
  














