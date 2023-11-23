# correspondence-analysis.R
# load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(edgeR)
  library(limma)
  library(Seurat)
  library(RRHO2)
  library(svglite)
  library(ggsvg)
  library(ggpubr)
  library(gridExtra)
})
#setwd("/nethome/kcni/karbabi/r_projects/cortical_microcircuits_across_psychiatric_disorders")
#options(device = "RstudioGD")
#options(scipen=999)

# concordance w. sinlge-nucleus ctrl ----
# load and format human single-cell reference (https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq)
sc_ref_metadata = fread(file = "./data/allen-brain/Human Multiple Cortical Areas SMART-seq/metadata.csv", data.table = F) %>% 
  filter(!(subclass_label == "" | is.na(subclass_label))) %>%
  filter(str_detect(cluster_label, "Exc L2-3 LINC00507|Exc L6 FEZF2") | subclass_label %in% c("VIP","SST","PVALB")) %>%
  mutate(CT = case_when(str_detect(cluster_label, "Exc L2-3 LINC00507") ~ "PyrL2n3", 
                        str_detect(cluster_label, "Exc L6 FEZF2") ~ "PyrL5n6",
                        TRUE ~ subclass_label)) %>%
  mutate(CT = factor(CT, levels =  c("PVALB","SST","VIP","PyrL2n3","PyrL5n6")),
         external_donor_name_label = factor(external_donor_name_label),
         pseudobulk_group = interaction(external_donor_name_label, CT))
rownames(sc_ref_metadata) = sc_ref_metadata$sample_name
table(sc_ref_metadata$CT)
# counts matrix 
sc_ref_matrix = fread(file = "./data/allen-brain/Human Multiple Cortical Areas SMART-seq/matrix.csv", data.table = F) %>%
  column_to_rownames(var = "sample_name")
sc_ref_matrix = sc_ref_matrix[match(sc_ref_metadata$sample_name, rownames(sc_ref_matrix)),]
# check
all(rownames(sc_ref_metadata) == rownames(sc_ref_matrix))

# get pseudobulk
mm = model.matrix(~ 0 + external_donor_name_label:CT, sc_ref_metadata)
dim(mm)
pseudobulk_matrix = t(sc_ref_matrix) %*% mm
colnames(pseudobulk_matrix) = gsub("external_donor_name_label", "", colnames(pseudobulk_matrix))
colnames(pseudobulk_matrix) = gsub(":CT", ".", colnames(pseudobulk_matrix))
# normalize 
pseudobulk_lcpm = edgeR::cpm(pseudobulk_matrix, prior.count = 2, log = TRUE)
pseudobulk_lcpm[1:20,1:10]

# get average expression of each gene for each cell type
celltypes = as.character(unique(sc_ref_metadata$CT))
sc_ref_df = bind_rows(
  lapply(celltypes, function(ct){
    res = data.frame(mean = apply(pseudobulk_lcpm[,grepl(ct, colnames(pseudobulk_lcpm))], 1, mean),
                     CT = ct) %>%
      rownames_to_column(var = "gene_symbol")
  })) %>%
  left_join(., ensembl, by = "gene_symbol") %>%
  distinct(gene_symbol, CT, .keep_all = TRUE)
# save
#write.csv(sc_ref_df, file = "./output/sc_ref_avg_exp.csv")

#load 
sc_ref_df = read.csv(file = "./output/sc_ref_avg_exp.csv", row.names = 1)
metadata_ctrl = metadata %>% filter(DX == "CTRL")
# get average expression of each gene for each cell type in controls 
lcm_exp_df = do.call(rbind,
                     by(metadata_ctrl, simplify = FALSE, metadata_ctrl$CT, function(x){
                       data.frame(mean = apply(lcpm_matrix[,x$ID], 1, mean),
                                  CT = as.character(unique(x$CT))) %>%
                         rownames_to_column(var = "gene_ensembl")
                     })) %>%
  left_join(., ensembl, by = "gene_ensembl") %>%
  distinct(CT, gene_symbol, .keep_all = TRUE)

# highlight DE genes
plot_df = lcm_exp_df %>%
  dplyr::select(CT, gene_symbol, mean) %>%
  left_join(., lcm_df_limma %>% dplyr::select(gene_symbol, CT, DE), by = c("CT","gene_symbol")) %>%
  left_join(., sc_ref_df %>% dplyr::select(CT, gene_symbol, mean), by = c("CT","gene_symbol")) %>%
  mutate(CT = factor(CT, levels = c("PVALB","SST","VIP","PyrL2n3","PyrL5n6"), 
                     labels = c("PVALB","SST","VIP","PyrL2n3","PyrL5n6"))) 
# visualize
p3 = ggplot(data = plot_df, aes(x = mean.x, y = mean.y)) +
  geom_point(data = plot_df %>% filter(DE == FALSE), aes(mean.x, mean.y), shape = 16, color = "darkgrey", alpha = 0.05) +
  geom_point(data = plot_df %>% filter(DE == TRUE), aes(mean.x, mean.y), shape = 16, color = "#1F77B4FF", alpha = 0.8) + 
  geom_smooth(data = plot_df, aes(mean.x, mean.y), method = "lm", se = FALSE, color = "#FF7F0EFF", lwd = 1.1) +
  stat_cor(method = "pearson") +
  labs(x = "LCM-seq gene-wise average expression\n(control only, logCPM)", y = "snRNA-seq gene-wise average expression\n(ref. atlas, logCPM)") +
  facet_wrap(~CT, nrow = 2, ncol = 3, scales = "free") +
  theme_classic2() +
  theme(strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),) 
#ggsave(filename = "./figures/draft/human_sc_comparisons_c.svg", height = 7.5, width = 9)

# get percent of cells that have non-zero expression of each gene 
sc_ref_df = do.call(rbind, 
                    by(sc_ref_metadata, sc_ref_metadata$CT, simplify = FALSE, function(x){
                      data.frame(non_zero_sum = apply(sc_ref_matrix[x$sample_name,], 2, function(x) sum(x!=0)),
                                 non_zero_pct = apply(sc_ref_matrix[x$sample_name,], 2, function(x) sum(x!=0)/length(x)),
                                 CT = as.character(unique(x$CT))) %>%
                        rownames_to_column(var = "gene_symbol")
                    })) %>%
  left_join(., ensembl %>% dplyr::select(gene_symbol, gene_ensembl), by = "gene_symbol") %>%
  mutate(EXP = non_zero_pct > 0.10) %>%
  distinct(gene_symbol, CT, .keep_all = TRUE)
#write.csv(sc_ref_df, file = "./output/sc_ref_non_zero.csv")

sc_ref_df = read.csv(file = "./output/sc_ref_non_zero.csv", row.names = 1)
lcm_df_limma = lcm_df_limma %>% 
  left_join(., sc_ref_df %>% dplyr::select(CT, gene_ensembl, non_zero_pct, EXP), by = c("CT","gene_ensembl")) %>% 
  mutate(DX = factor(DX, levels = c("MDD","BD","SCZ"), labels = c("MDD","BD","SCZ"))) %>%
  mutate(CT = factor(CT, levels = c("PVALB","SST","VIP","PyrL2n3","PyrL5n6"))) 
#write_rds(lcm_df_limma, file = "./output/lcm_df_limma.rds")

cols = pal_d3("category10")(10)[-c(1:2)]

p1 = lcm_df_limma %>%
  filter(DE == TRUE) %>%
  mutate(gene_gene_biotype = gsub("_", " ", gene_biotype)) %>%
  group_by(CT, gene_biotype) %>% 
  dplyr::summarize(Freq = n()) %>%
  drop_na() %>%
  filter(Freq > 2) %>%
  ggplot(., aes(x = CT, y = Freq, fill = reorder(gene_biotype, -Freq))) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.8, preserve = "single")) + 
  geom_text(aes(label = Freq), position = position_dodge2(width = 0.8, preserve = "single"), vjust = -0.25) +
  scale_y_continuous(expand = c(0,0), limits = c(0,180)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(ncol = 1, byrow = FALSE)) +
  labs(x = "", y = "", fill = "Transcript class") + 
  theme_classic2() +
  theme(strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        legend.position = "right", legend.justification = "left", legend.margin = margin(0,0,0,0))

p2 = lcm_df_limma %>%
  filter(DE == TRUE) %>%
  group_by(CT, EXP) %>% 
  dplyr::summarize(Freq = n()) %>%
  drop_na() %>%
  ggplot(., aes(x = CT, y = Freq, fill = reorder(EXP, -Freq))) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.8, preserve = "single")) + 
  geom_text(aes(label = Freq), position = position_dodge2(width = 0.8, preserve = "single"), vjust = -0.25) +
  scale_y_continuous(expand = c(0,0), limits = c(0,180)) +
  scale_fill_d3() +
  labs(x = "Cell type", y = "Number of DE genes (FDR < 0.20)", fill = "Detected expression in snRNA-seq\nref. atlas (>10% of cells)") + 
  theme_classic2() +
  theme(strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        legend.position = "right", legend.justification = "left", legend.margin = margin(0,0,0,0))

#svglite(filename = paste("./figures/draft/human_sc_comparisons.svg"), height = 7.5, width = 16)
((p1 / p2) | p3) +
  plot_layout(widths = c(1,1.3))
#dev.off()

# correspondence w. single-nucleus scz ----
load(file = "./data/lcm-seq/data.Rdata")
lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds")

## load pseudobulk data ----
pb_metadata = fread(file = "./data/mtsinai-mclean/pseudobulk_fine_grained_metadata.tsv", data.table = F) %>%
  mutate(Bulk = paste(unique_donor_ID, Celltype, sep = "_")) %>%
  mutate(Phenotype = factor(Phenotype, levels = c("CON","SZ"))) 
rownames(pb_metadata) = pb_metadata$Bulk

pb_counts = fread(file = "./data/mtsinai-mclean/pseudobulk_fine_grained_counts.tsv", data.table = F) %>%
  unite(unique_pb_ID, c(unique_donor_ID, Celltype), sep = "_", remove = TRUE) %>%
  column_to_rownames(var = "unique_pb_ID")
pb_counts = t(pb_counts)

pb_counts = pb_counts[, match(rownames(pb_metadata), colnames(pb_counts))]
all(rownames(pb_metadata) == colnames(pb_counts))

# calculate scaling factors to convert raw library sizes into effective library sizes
dge0 = DGEList(pb_counts)
dge0 = calcNormFactors(dge0, method = "TMM")
# add to metadata
pb_metadata = pb_metadata %>%
  left_join(., dge0$samples %>% rownames_to_column(var = "Bulk") %>% dplyr::select(Bulk, lib.size), by = "Bulk") %>%
  dplyr::rename(num_reads = lib.size)

## celltype stratified limma-voom ----
# run for each cell type separately 
pb_limma = do.call(rbind,
                   by(pb_metadata, list(pb_metadata$Celltype), function(x){
                     # only include subjects with a minimum of 5 cells in the cluster
                     x = x %>% filter(num_cells >= 5)
                     if(nrow(x) > 10) {
                       # keep genes expressed in at least x% of samples
                       dge = dge0[,x$Bulk]
                       dge = dge[rowSums(dge$counts >= 1) >= ncol(dge)*0.4, ]
                       
                       design = model.matrix(~ Phenotype + Gender + Age + PMI + Cohort + Batch + log(num_cells) + log(num_reads), x) 
                       colnames(design) = gsub("\\(|\\)","", colnames(design))
                       
                       vm = voom(dge, design = design, plot = TRUE)
                       title(main = unique(x$Celltype), cex.main = 2, line = -2)
                       fit = lmFit(vm, design)
                       fit = eBayes(fit)
                       
                       print(unique(x$Celltype))
                       print(summary(decideTests(fit, p.value = 0.05, adjust.method = "BH"))[,2])
                       
                       tt = topTable(fit, coef = "PhenotypeSZ",  n = Inf, sort = "none", adjust.method = "BH")
                       tt$gene_symbol = rownames(tt)
                       tt$Celltype = as.character(unique(x$Celltype))
                       return(tt)
                       }
                     })
)
rownames(pb_limma) = NULL
#write.csv(pb_limma, file = "./data/mtsinai-mclean/de_limma.csv", row.names = FALSE)

## find best cell type matches ----
combs = as.matrix(expand.grid(as.character(unique(lcm_df_limma$CT)), as.character(unique(pb_limma$Celltype))))
combs = split(combs, seq(nrow(combs)))

res = lapply(combs, function(x){
  v1 = lcm_df_limma %>% 
    filter(DX == "SCZ", CT == x[1]) %>% 
    pull(t, gene_symbol)
  v2 = pb_limma %>% 
    filter(Celltype == x[2]) %>% 
    pull(t, gene_symbol)
  v1 = v1[intersect(names(v1), names(v2))]
  v2 = v2[intersect(names(v1), names(v2))]
  data.frame(ct1 = x[1], ct2 = x[2], 
             n = length(v1),
             cor = cor(v1,v2, method = "spearman"), 
             pval = cor.test(v1,v2)$p.value)
})
res = data.frame(do.call(rbind, res))

## RRHO plots ----
# http://htmlpreview.github.io/?https://github.com/RRHO2/RRHO2/blob/master/vignettes/RRHO2.html
pb_limma = read.csv(file = "./data/mtsinai-mclean/de_limma.csv")
pb_limma_filt = pb_limma %>%
  mutate(Celltype = case_match(Celltype,
                              'Ex-L23'~'PyrL2n3', 
                              'Ex-L56'~'PyrL5n6', 
                              'In-PV_Chandelier'~'PVALB', 
                              'In-SST'~"SST", 
                              'In-VIP'~'VIP')) %>%
  filter(Celltype %in% c("PVALB","SST","VIP","PyrL2n3","PyrL5n6")) %>%
  mutate(Celltype = factor(Celltype, levels = c("PVALB","SST","VIP","PyrL2n3","PyrL5n6"))) 

plot_lst = NULL
for(ct in unique(lcm_df_limma$CT)){
  # prepare ranked gene lists
  tmp_lcm = lcm_df_limma %>%
    filter(CT == ct, DX == "SCZ") %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    mutate(rank = -log10(P.Value)*sign(logFC),
           logFDR = -log10(adj.P.Val)*sign(logFC)) %>%
    arrange(dplyr::desc(rank)) 
  rank_lcm = tmp_lcm %>% dplyr::select(gene_symbol, rank)
  
  tmp_sc = pb_limma_filt %>%
    filter(Celltype == ct) %>%
    mutate(rank = -log10(P.Value)*sign(logFC),
           logFDR = -log10(adj.P.Val)*sign(logFC)) %>%
    arrange(dplyr::desc(rank)) 
  rank_sc = tmp_sc %>% dplyr::select(gene_symbol, rank)
  
  # same genes only
  keep_genes = intersect(rank_lcm$gene_symbol, rank_sc$gene_symbol)
  rank_lcm = rank_lcm %>% filter(gene_symbol %in% keep_genes)
  rank_sc = rank_sc %>% filter(gene_symbol %in% keep_genes)
  
  rrho_res = RRHO2_initialize(rank_lcm, rank_sc, labels = c(paste(ct,"LCM-seq"), paste(ct,"single-nucleus RNA-seq")), log10.ind = TRUE, method = "hyper")
  svglite(paste0("./figures/draft/rrho/sc_scz_rrho",ct,".svg"), scaling = 1.6)
  if(ct %in% c("PyrL2n3","PyrL5n6")){
    threshold = 20
  } else {
    threshold = 8
  }
  RRHO2_heatmap(rrho_res, maximum = threshold, minimum = 0, useRaster = TRUE)
  dev.off()
  plot_lst[[ct]] = svg_to_rasterGrob(paste(readLines(paste0("./figures/draft/rrho/sc_scz_rrho",ct,".svg")), collapse = "\n"))
  
  plot(
    inner_join(tmp_lcm %>% dplyr::select(gene_symbol, logFDR),
               tmp_sc %>% dplyr::select(gene_symbol, logFDR), 
               by = "gene_symbol", suffix = c("_lcm","_sc")) %>%
      mutate(col = case_when(logFDR_lcm > -log10(0.2) & logFDR_sc > -log10(0.2) ~ "#D62728FF",
                             logFDR_lcm < log10(0.2) & logFDR_sc < log10(0.2) ~ "#1F77B4FF",
                             .default = "grey"),
             lab = case_when(abs(logFDR_lcm) > -log10(0.4) | abs(logFDR_sc) > -log10(0.4) ~ gene_symbol)) %>%
      ggplot(., aes(x = logFDR_lcm, y = logFDR_sc, color = col, label = lab)) +
      geom_point() +
      ggrepel::geom_text_repel() + 
      scale_color_identity() +
      scale_x_continuous(trans = "log10") +
      scale_y_continuous(trans = "log10") +
      labs(title = paste(ct)) +
      theme_classic()
  )
  ggsave(filename = paste0("./figures/draft/rrho/scatter_scz_",ct,".png"), width = 5, height = 5)
  
}
svglite(filename = "./figures/draft/sc_scz_rrho.svg", width = 20, height = 5)
grid.arrange(grobs = plot_lst, ncol = 5, nrow = 1)
dev.off()

# correspondence w. bulk scz ----

load(file = "./data/lcm-seq/data.Rdata")
lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds")

# load differential expression results from Gandal et al. 2018
rank_bulk = fread(file = "./data/gandal/RNAseq_CMC.csv", data.table = F) %>%
  rename(gene_ensembl = V1) %>%
  dplyr::select(gene_ensembl, starts_with("SCZ"))
names(rank_bulk) = gsub("SCZ.", "", names(rank_bulk))

rank_bulk = rank_bulk %>%
  distinct(gene_ensembl, .keep_all = TRUE) %>%
  mutate(rank = -log10(P.Value)*sign(logFC)) %>%
  arrange(dplyr::desc(rank)) %>%
  dplyr::select(gene_ensembl, rank)

# for each cell type
plot_lst = NULL
for(ct in unique(lcm_df_limma$CT)){
  # prepare ranked gene lists
  rank_lcm = lcm_df_limma %>%
    filter(CT == ct, DX == "SCZ") %>%
    distinct(gene_ensembl, .keep_all = TRUE) %>%
    mutate(rank = -log10(P.Value)*sign(logFC)) %>%
    arrange(dplyr::desc(rank)) %>%
    dplyr::select(gene_ensembl, rank)
  # same genes only
  keep_genes = intersect(rank_lcm$gene_ensembl, rank_bulk$gene_ensembl)
  rank_lcm = rank_lcm %>% filter(gene_ensembl %in% keep_genes)
  rank_bulk = rank_bulk %>% filter(gene_ensembl %in% keep_genes)
  
  # Rank Rank Hypergeometric Overlap
  # http://htmlpreview.github.io/?https://github.com/RRHO2/RRHO2/blob/master/vignettes/RRHO2.html
  rrho_res = RRHO2_initialize(rank_lcm, rank_bulk, labels = c(paste(ct,"LCM-seq"), "bulk tissue RNA-seq"), log10.ind = TRUE, method = "hyper")
  svglite(paste0("./figures/draft/rrho/bulk_scz_rrho",ct,".svg"), scaling = 1.6)
  if(ct %in% c("PyrL2n3","PyrL5n6")){
    threshold = 30
  } else {
    threshold = 10
  }
  RRHO2_heatmap(rrho_res, maximum = threshold, minimum = 0, useRaster = TRUE)
  dev.off()
  plot_lst[[ct]] = svg_to_rasterGrob(paste(readLines(paste0("./figures/draft/rrho/bulk_scz_rrho",ct,".svg")), collapse = "\n"))
}
svglite(filename = "./figures/draft/bulk_scz_rrho.svg", width = 20, height = 5)
grid.arrange(grobs = plot_lst, ncol = 5, nrow = 1)
dev.off()

# correspondence w. single-nucleus mdd ----

load(file = "./data/lcm-seq/data.Rdata")
lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds")

## prepare data ----
# wget https://cells.ucsc.edu/dlpfc-mdd/matrix.mtx.gz # matrix.mtx.gz
# wget https://cells.ucsc.edu/dlpfc-mdd/features.tsv.gz # features.tsv.gz
# wget https://cells.ucsc.edu/dlpfc-mdd/barcodes.tsv.gz # barcodes.tsv.gz
# wget https://cells.ucsc.edu/dlpfc-mdd/meta.tsv # metadata

# load sparse expression matrix
sc_matrix = ReadMtx(
  mtx = "./data/maitra/matrix.mtx.gz", 
  features = "./data/maitra/features.tsv.gz", feature.column = 1,
  cells = "./data/maitra/barcodes.tsv.gz"
)
# load cell metdata
sc_metadata = fread(file = "./data/maitra/meta.tsv", data.table = FALSE) %>%
  mutate(ID = sub("\\..*", "", Cell)) %>%
  mutate(Bulk = as.character(interaction(ID, Cluster))) %>%
  group_by(Bulk) %>%
  mutate(N_cells = n()) %>%
  as.data.frame()
rownames(sc_metadata) = sc_metadata$Cell

# metadata for pseudobulks 
pb_metadata = sc_metadata %>%
  distinct(Bulk, .keep_all = TRUE) %>%
  mutate(Condition = factor(Condition, levels = c("Control","Case"))) %>%
  as.data.frame()
rownames(pb_metadata) = NULL
# filter 
sc_metadata = sc_metadata %>% 
  filter(Bulk %in% pb_metadata$Bulk,
         Cluster %in% pb_metadata$Cluster)
sc_matrix = sc_matrix[,match(rownames(sc_metadata),colnames(sc_matrix))]
all(rownames(sc_metadata) == colnames(sc_matrix))

# get pseudobulks
mm = model.matrix(~ 0 + Bulk, sc_metadata)
pb_matrix = as.matrix(sc_matrix %*% mm)
colnames(pb_matrix) = gsub("Bulk", "", colnames(pb_matrix))

# calculate scaling factors to convert raw library sizes into effective library sizes
dge0 = DGEList(pb_matrix)
dge0 = calcNormFactors(dge0, method = "TMM")
pb_metadata = left_join(pb_metadata, 
                        dge0$samples %>% rownames_to_column(var = "Bulk") %>% dplyr::select(Bulk, lib.size), by = "Bulk") %>%
  mutate(Cluster_comb = case_when(
    Cluster %in% c("InN1_PV","InN9_PV") ~ "PV",
    Cluster %in% c("InN2_SST","InN5_SST") ~ "SST",
    Cluster %in% c("InN4_VIP","InN3_VIP") ~ "VIP",
    Cluster %in% c("ExN1_L23","ExN2_L23","ExN9_L23") ~ "L23 IT",
    Cluster %in% c("ExN11_L56","ExN12_L56","ExL13_L56") ~ "L6b",
    TRUE ~ Cluster
  )) %>%
  rename(N_reads = lib.size)
# quick save
#save(sc_metadata, pb_metadata, pb_matrix, dge0, file = "./data/maitra/pseudobulks.RData")

## celltype stratified limma-voom ----
load(file = "./data/maitra/pseudobulks.RData")
# run for each cell type separately
pb_limma = 
  do.call(rbind, 
          by(pb_metadata, list(pb_metadata$Broad), function(y) {
            # only include subjects with a minimum of 10 cells in the broad cell type (as per authors) 
            y = y %>% filter(N_cells >= 10)
            do.call(rbind,
                    by(y, list(y$Cluster), function(x) {
                      # only include subjects with a minimum of 5 cells in the cluster
                      x = x %>% filter(N_cells >= 5)
                      # keep genes expressed in at least x% of samples
                      dge = dge0[,x$Bulk]
                      dge = dge[rowSums(dge$counts >= 1) >= ncol(dge)*0.3, ]
                      # drop categorical covairates with one level 
                      drop = names(x)[sapply(x, function(z) length(unique(z)) == 1)]
                      model = '~ Condition + Sex + Chemistry + log(N_reads) + log(N_cells)'
                      if (!is.null(drop)) {
                        drop = drop[drop %in% c('Sex', 'Chemistry')]
                        if (length(drop) > 0) {
                          model = as.formula(gsub(paste('+', paste(drop, collapse = " + ")), "", model, fixed = TRUE))
                        }
                      }
                      # note: we are missing metadata on Age and PMI
                      design = model.matrix(as.formula(model), x)
                      colnames(design) = gsub("\\(|\\)","", colnames(design))
                      
                      vm = voom(dge, design = design, plot = TRUE)
                      title(main = unique(x$Cluster), cex.main = 2, line = -2)
                      fit = lmFit(vm, design)
                      fit = eBayes(fit)
                      
                      print(unique(x$Cluster))
                      print(summary(decideTests(fit, p.value = 0.05, adjust.method = "BH"))[,2])
                      
                      tt = topTable(fit, coef = "ConditionCase", n = Inf, sort = "none", adjust.method = "BH")
                      tt$gene_symbol = rownames(tt)
                      tt$Cluster = as.character(unique(x$Cluster))
                      return(tt)
                    })
            )
          })
  )
rownames(pb_limma) = NULL
#write.csv(pb_limma, file = "./data/maitra/de_limma.csv", row.names = FALSE)

## find best cell type matches ----
pb_limma = read.csv(file = "./data/maitra/de_limma.csv")
combs = as.matrix(expand.grid(as.character(unique(lcm_df_limma$CT)), as.character(unique(pb_limma$Cluster))))
combs = split(combs, seq(nrow(combs)))

res = do.call(
  rbind, lapply(combs, function(x){
    v1 = lcm_df_limma %>% 
      filter(DX == "MDD", CT == x[1]) %>% 
      #filter(abs(t) > 1) %>%
      pull(t, gene_symbol)
    v2 = pb_limma %>% 
      filter(Cluster == x[2]) %>% 
      #filter(abs(t) > 1) %>%
      pull(t, gene_symbol)
    v1 = v1[intersect(names(v1), names(v2))]
    v2 = v2[intersect(names(v1), names(v2))]
    data.frame(ct1 = x[1], ct2 = x[2], 
               n = length(v1),
               cor = cor(v1,v2, method = "spearman"), 
               pval = cor.test(v1,v2)$p.value)
  })
)

## RRHO plots ----

# InN1_PV, InN9_PV
# InN2_SST, Inn5_SST
# InN3_VIP, InN4_VIP
pb_limma = read.csv(file = "./data/maitra/de_limma.csv")
pb_limma_filt = pb_limma %>%
  mutate(Cluster = case_match(Cluster,
                              'ExN11_L56'~'PyrL2n3', 
                              'ExN16_L56'~'PyrL5n6', 
                              'InN9_PV'~'PVALB', 
                              'InN5_SST'~"SST", 
                              'InN3_VIP'~'VIP')) %>%
  filter(Cluster %in% c("PVALB","SST","VIP","PyrL2n3","PyrL5n6")) %>%
  mutate(Cluster = factor(Cluster, levels = c("PVALB","SST","VIP","PyrL2n3","PyrL5n6"))) 

plot_lst = NULL
for(ct in unique(lcm_df_limma$CT)){
  # prepare ranked gene lists
  tmp_lcm = lcm_df_limma %>%
    filter(CT == ct, DX == "MDD") %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    mutate(rank = -log10(P.Value)*sign(logFC),
           logFDR = -log10(adj.P.Val)*sign(logFC)) %>%
    arrange(dplyr::desc(rank)) 
  rank_lcm = tmp_lcm %>% dplyr::select(gene_symbol, rank)
  
  tmp_sc = pb_limma_filt %>%
    filter(Cluster == ct) %>%
    mutate(rank = -log10(P.Value)*sign(logFC),
           logFDR = -log10(adj.P.Val)*sign(logFC)) %>%
    arrange(dplyr::desc(rank)) 
  rank_sc = tmp_sc %>% dplyr::select(gene_symbol, rank)
  
  # same genes only
  keep_genes = intersect(rank_lcm$gene_symbol, rank_sc$gene_symbol)
  rank_lcm = rank_lcm %>% filter(gene_symbol %in% keep_genes)
  rank_sc = rank_sc %>% filter(gene_symbol %in% keep_genes)
  
  rrho_res = RRHO2_initialize(rank_lcm, rank_sc, labels = c(paste(ct,"LCM-seq"), paste(ct,"single-nucleus RNA-seq")), log10.ind = TRUE, method = "hyper")
  svglite(paste0("./figures/draft/rrho/sc_mdd_rrho",ct,".svg"), scaling = 1.6)
  if(ct %in% c("PyrL2n3","PyrL5n6")) {
    threshold = 20
  } else {
    threshold = 8
  }
  RRHO2_heatmap(rrho_res, maximum = threshold, minimum = 0, useRaster = TRUE)
  dev.off()
  plot_lst[[ct]] = svg_to_rasterGrob(paste(readLines(paste0("./figures/draft/rrho/sc_mdd_rrho",ct,".svg")), collapse = "\n"))
  
  # plot(
  #   inner_join(tmp_lcm %>% dplyr::select(gene_symbol, logFDR),
  #              tmp_sc %>% dplyr::select(gene_symbol, logFDR), 
  #              by = "gene_symbol", suffix = c("_lcm","_sc")) %>%
  #     mutate(col = case_when(logFDR_lcm > -log10(0.2) & logFDR_sc > -log10(0.2) ~ "#D62728FF",
  #                            logFDR_lcm < log10(0.2) & logFDR_sc < log10(0.2) ~ "#1F77B4FF",
  #                            .default = "grey"),
  #            lab = case_when(abs(logFDR_lcm) > -log10(0.4) | abs(logFDR_sc) > -log10(0.4) ~ gene_symbol)) %>%
  #     ggplot(., aes(x = logFDR_lcm, y = logFDR_sc, color = col, label = lab)) +
  #     geom_point() +
  #     ggrepel::geom_text_repel() + 
  #     scale_color_identity() +
  #     scale_x_continuous(trans = "log10") +
  #     scale_y_continuous(trans = "log10") +
  #     labs(title = paste(ct)) +
  #     theme_classic()
  # )
  # ggsave(filename = paste0("./figures/draft/rrho/scatter_mdd_",ct,".png"), width = 5, height = 5)
}
#svglite(filename = "./figures/draft/sc_mdd_rrho.svg", width = 20, height = 5)
grid.arrange(grobs = plot_lst, ncol = 5, nrow = 1)
#dev.off()

# # create Seurat objects for easier filtering of cells (qc)
# sc_metadata = fread(file = "./data/maitra/meta.tsv", data.table = FALSE) %>%
#   mutate(Split_by = paste(Sex, Chemistry, sep = "."))
# seu_obj = CreateSeuratObject(counts = sc_matrix, meta.data = sc_matrix) 
# seu_obj[["percent.mt"]] = PercentageFeatureSet(seu_obj, pattern = "^MT-")
# # separate thresholds for each sex and chemistry (from the authors)
# seu_lst = SplitObject(seu_obj, split.by = "Split_by")
# dim(seu_lst$Female.v3)
# seu_lst$Female.v3 = subset(seu_lst$Female.v3, subset = nCount_RNA < 120000 & nFeature_RNA > 350 & percent.mt < 10)
# dim(seu_lst$Female.v3)
# seu_lst$Female.v2 = subset(seu_lst$Female.v2, subset = nCount_RNA < 25000 & nFeature_RNA > 250 & percent.mt < 10)
# seu_lst$Male.v2 = subset(seu_lst$Male.v2, subset = nCount_RNA < 35000 & nFeature_RNA > 350 & percent.mt < 10)
# seu_obj = merge(seu_lst$Female.v3, c(seu_lst$Female.v2, seu_lst$Male.v2))
# sc_metadata = seu_obj@meta.data
# sc_matrix = seu_obj@assays$RNA@counts


# correspondence analysis w. bulk mdd ----
load(file = "./data/lcm-seq/data.Rdata")
lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds")

# read male and female DE lists from Labonte et al. 2017
bulk_m = read_excel("./data/labonte/Male_DEG_all.xlsx") %>%
  rename(Brain_region = "Brain Region", p_value = "p-value", gene_ensembl = "ENSG") %>%
  filter(Brain_region == "BA25")
bulk_f = read_excel("./data/labonte/Female_DEG_all.xlsx") %>%
  rename(Brain_region = "Brain region", p_value = "p-value", gene_ensembl = "ENSG") %>%
  filter(Brain_region == "BA25")

#rank_bulk = bulk_m %>%
rank_bulk = bulk_f %>%
  distinct(gene_ensembl, .keep_all = TRUE) %>%
  mutate(rank = -log10(p_value)*sign(logFC)) %>%
  arrange(dplyr::desc(rank)) %>%
  dplyr::select(gene_ensembl, rank) %>%
  as.data.frame()

# for each cell type
plot_lst = NULL
for(ct in unique(lcm_df_limma$CT)){
  # prepare ranked gene lists
  rank_lcm = lcm_df_limma %>%
    filter(CT == ct, DX == "MDD") %>%
    distinct(gene_ensembl, .keep_all = TRUE) %>%
    mutate(rank = -log10(P.Value)*sign(logFC)) %>%
    arrange(dplyr::desc(rank)) %>%
    dplyr::select(gene_ensembl, rank)
  # same genes only
  keep_genes = intersect(rank_lcm$gene_ensembl, rank_bulk$gene_ensembl)
  rank_lcm = rank_lcm %>% filter(gene_ensembl %in% keep_genes)
  rank_bulk = rank_bulk %>% filter(gene_ensembl %in% keep_genes)
  
  # Rank Rank Hypergeometric Overlap
  # http://htmlpreview.github.io/?https://github.com/RRHO2/RRHO2/blob/master/vignettes/RRHO2.html
  rrho_res = RRHO2_initialize(rank_lcm, rank_bulk, labels = c(paste(ct,"LCM-seq"), "bulk tissue RNA-seq"), log10.ind = TRUE, method = "hyper")
  svglite(paste0("./figures/draft/rrho/bulk_mdd_rrho",ct,".svg"), scaling = 1.6)
  if(ct %in% c("PyrL2n3","PyrL5n6")){
    threshold = 15
  } else {
    threshold = 5
  }
  RRHO2_heatmap(rrho_res, maximum = threshold, minimum = 0, useRaster = TRUE)
  dev.off()
  plot_lst[[ct]] = svg_to_rasterGrob(paste(readLines(paste0("./figures/draft/rrho/bulk_mdd_rrho",ct,".svg")), collapse = "\n"))
}
svglite(filename = "./figures/draft/bulk_mdd_rrho_f.svg", width = 20, height = 5)
grid.arrange(grobs = plot_lst, ncol = 5, nrow = 1)
dev.off()


























