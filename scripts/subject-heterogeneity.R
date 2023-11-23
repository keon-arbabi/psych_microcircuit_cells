
# subject-heterogeneity.R 
# load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(xgboost)
  library(caret)
  library(edgeR)
  library(limma)
  library(pROC)
  library(ppcor)
  library(GWENA)
  library(pheatmap)
  library(ggpubr)
  library(patchwork)
  library(RNASeqBits)
  library(scales)
  library(ggsci)
  library(svglite)
})
cols_1 = c("PyrL2n3" = "#D62728FF", "PyrL5n6" = "#FF7F0EFF", "PVALB" = "#2CA02CFF", "SST" = "#9467BDFF", "VIP" = "#1F77B4FF")

# calculate TRS ----
load(file = "./data/lcm-seq/data.Rdata")
# add star qc 
metadata = metadata %>% left_join(., star_qc, by = 'ID')
rownames(metadata) = metadata$ID
# load DE results
lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds")

disorders = as.character(unique(metadata$DX))[-1]
cell_types = as.character(unique(metadata$CT))

trs = do.call(rbind,
              lapply(cell_types, function(ct){
                do.call(rbind,
                        sapply(disorders, simplify = F, function(dx){
                          
                          meta = metadata %>% filter(CT == ct, DX %in% c("CTRL",dx))
                          labs = as.numeric(!meta$DX == "CTRL")
                          
                          dge = DGEList(counts_matrix[,meta$ID])
                          dge = cpm(calcNormFactors(dge))
                          
                          keep = rowSums(dge > 1) >= 15
                          #dge = dge[keep,]
                          dge = t(filter_low_var(t(dge), pct = 0.2, type = "mean"))
                          #dge = dge[rownames(dge) %in% de_genes,]
                          
                          mod = model.matrix(~ Sex + scale(Intergenic_rate), meta)
                          
                          df = data.frame(ID = meta$ID, CT = ct, DX = dx,
                                          pred = sapply(seq(length(meta$ID)), function(ind){
                                            
                                            lfit = lmFit(dge[,-ind], mod[-ind,])
                                            #dat = t(dge) - (mod %*% t(lfit$coefficients))
                                            dat = t(dge)
                                            
                                            # wts = ifelse(labs[-ind], sum(1 - labs[-ind]), sum(labs[-ind]))
                                            # wts = wts/mean(wts)
                                            set.seed(0)
                                            fit = xgboost(data = dat[-ind,], label = labs[-ind], nrounds = 100, objective = "binary:logistic", verbose = F)
                                            
                                            test = matrix(dat[ind,], nrow = 1)
                                            predict(fit, test)
                                          }))
                          df$auc = roc(labs, df$pred, direction = "<")$auc
                          df
                        }))
              })
)
trs_df = trs %>%
  pivot_wider(id_cols = ID, names_from = DX, values_from = pred, names_glue = "{DX}_trs")
#write.csv(trs_df, file = "./output/trs_df.csv", row.names = F)

# combine all metadata ----
load(file = "./data/lcm-seq/data.Rdata")
metadata = metadata %>% left_join(., star_qc, by = 'ID')

# load and tidy CommonMind metadata
metadata_cmc = readRDS(file = "./data/cmc/METADATA.rds") %>%
  filter(Institution == "Pitt") %>%
  dplyr::rename(HU = SampleID, DX = Dx) %>%
  mutate(HU = str_replace(HU, "_BP_PFC", "")) %>%
  mutate(HU = str_replace(HU, "_PFC", "")) %>%
  mutate(HU = sapply(strsplit(HU, "_"), "[[", 3)) %>%
  filter(HU %in% metadata$HU) %>%
  distinct(HU, .keep_all = TRUE) %>%
  dplyr::select(HU, PMI, RIN, IntronicRate, IntragenicRate, IntergenicRate, AlignmentRate, GenesDetected) %>%
  dplyr::rename(PMI_bulk = PMI, RIN_bulk = RIN, Intronic_rate_bulk = IntronicRate, Intragenic_rate_bulk = IntragenicRate, 
                Intergenic_rate_bulk = IntergenicRate, Mapped_pct_bulk = AlignmentRate, Library_size_bulk = GenesDetected)

# load psychiatric PRS scores for CommonMind 
prs_scores = read.csv(file = "./data/cmc/PRSCS_cmcrodrigo.csv") %>%
  dplyr::rename(HU = ID) %>%
  filter(str_detect(HU, 'Pitt|PITT')) %>%
  mutate(HU = str_replace(HU, "_BP", "")) %>%
  mutate(HU = sapply(strsplit(HU, "_"), "[[", 3)) %>%
  dplyr::select(HU, SZ, BD, MDD, IQ) %>%
  dplyr::rename(SCZ_prs = SZ, BD_prs = BD, MDD_prs = MDD, IQ_prs = IQ)

# load and tidy MGP deconvolution of CommonMind bulk RNA-seq 
cell_deconvolution = read.csv(file = "./output/cmc_deconvolution.csv", row.names = 1) %>%
  dplyr::rename(HU = SampleID, subclass_label = CT) %>%
  mutate(HU = str_replace(HU, "_BP_PFC", ""),
         CT = case_when(subclass_label == "IT" ~ "PyrL2n3",
                        subclass_label == "L6b" ~ "PyrL5n6",
                        subclass_label == "PVALB" ~ "PVALB",
                        subclass_label == "VIP" ~ "VIP",
                        subclass_label == "SST" ~ "SST")) %>%
  mutate(HU = str_replace(HU, "_PFC", ""),
         CT = factor(CT, levels = c("PyrL2n3","PyrL5n6","PVALB","VIP","SST")),
         DX = factor(DX, levels = c("CTRL","MDD","BD","SCZ"))) %>%
  mutate(HU = sapply(strsplit(HU, "_"), "[[", 3)) %>%
  drop_na() %>%
  distinct(HU, CT, .keep_all = TRUE) %>%
  dplyr::select(HU, CT, subclass_label, MGP)

# load and tidy FISH cell densities 
cell_densities = read.csv(file = "./output/FISH_densities.csv", row.names = 1) %>%
  dplyr::rename(Resid_density = res) %>%
  mutate(HU = as.character(as.numeric(HU))) %>%
  dplyr::select(HU, CT, Density, Log_density, Resid_density)

# combine all 
metadata_full = metadata %>%
  left_join(., metadata_cmc, by = "HU") %>%
  left_join(., prs_scores, by = "HU") %>%
  left_join(., cell_deconvolution, by = c("HU","CT")) %>%
  left_join(., cell_densities, by = c("HU","CT")) %>%
  left_join(., trs_df, by = "ID") %>%
  mutate(CT = factor(CT, levels =  c("PVALB","SST","VIP","PyrL2n3","PyrL5n6")),
         DX = factor(DX, levels = c("CTRL","MDD","BD","SCZ"))) 

names(metadata_full)


metadata_full %>%
  filter(DX %in% c("BD","CTRL")) %>%
  ggplot(., aes(x = DX, y = BD_trs)) +
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap(~CT) +
  theme_classic2()


p1 = metadata_full %>%
  filter(DX %in% c("SCZ","CTRL")) %>%
  ggplot(., aes(x = SCZ_trs, y = SCZ_prs, color = CT)) +
  geom_smooth(aes(color = CT), method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y.npc = 0.7) +
  scale_color_manual(values = cols_1) +
  theme_classic2() +
  theme(legend.position = "none")
p2 = metadata_full %>%
  mutate(Class = case_when(CT %in% c("PyrL2n3","PyrL5n6") ~ "EXC", TRUE ~ "INH")) %>%
  filter(DX %in% c("SCZ","CTRL")) %>%
  ggplot(., aes(x = SCZ_trs, y = Resid_density, color = CT)) +
  geom_smooth(aes(color = CT), method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y.npc = 0.7) +
  scale_color_manual(values = cols_1) +
  facet_wrap(~Class, scales = "free_y", nrow = 2) + 
  theme_classic2() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())
p3 = metadata_full %>%
  filter(DX %in% c("SCZ","CTRL")) %>%
  ggplot(., aes(x = SCZ_trs, y = Intergenic_rate, color = CT)) +
  geom_smooth(aes(color = CT), method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y.npc = 0.7) +
  scale_color_manual(values = cols_1) +
  theme_classic2()

metadata_full %>%
  filter(DX %in% c("SCZ","CTRL")) %>%
  ggplot(., aes(x = SCZ_trs, y = RIN, color = CT)) +
  geom_smooth(aes(color = CT), method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y.npc = 0.7) +
  scale_color_manual(values = cols_1) +
  theme_classic2()

svglite("./figures/draft/sh_scatter_scz.svg", height = 3.5, width = 9)
p1 + p2 + p3 + plot_layout(nrow = 1)
dev.off()


trs_mat = metadata_full %>%
  filter(DX %in% c("BD","CTRL")) %>%
  dplyr::select(HU, CT, BD_trs) %>%
  pivot_wider(names_from = HU, values_from = BD_trs) %>%
  column_to_rownames(var = "CT")

prs = metadata_full %>% filter(DX %in% c("BD","CTRL"), CT == "PyrL2n3") %>% pull(BD_prs) 
prs = rescale(prs, to = c(min(na.omit(unlist(trs_mat))), max(na.omit(unlist(trs_mat)))))
trs_mat = rbind(trs_mat, 
                Donor_trs = colMeans(trs_mat, na.rm = TRUE),
                Donor_prs = prs)  
trs_mat = trs_mat[,order(trs_mat[6,])]

# color scale
my_breaks = seq(0,1,0.01)
colors = c(rev(colorRampPalette(c("#FFFFFF", "#1F77B4FF"))(length(my_breaks)/2)), colorRampPalette(c("#FFFFFF", "#D62728FF"))(length(my_breaks)/2)) 

anno_col = metadata %>% filter(DX %in% c("BD","CTRL"), CT == "PyrL2n3") %>% arrange(factor(HU, levels = names(trs_mat))) %>% pull(DX) %>% as.character() %>% as.data.frame()
rownames(anno_col) = names(trs_mat)
names(anno_col) = "Diagnosis"

anno_colors = c("lightgrey","black")
names(anno_colors) = unique(anno_col$Diagnosis)
anno_colors = list(group = anno_colors)

svglite("./figures/draft/sh_heatmap_bd.svg", height = 4, width = 9)
pheatmap(as.matrix(trs_mat),
         display_numbers = FALSE,
         color = colors,
         breaks = my_breaks,
         border_color = "black",
         cluster_rows = F, cluster_cols = F,
         show_colnames = F,
         gaps_row = 5,
         annotation_col = anno_col,
         cellheight = 10, cellwidth = 10)
dev.off()


p1 = metadata_full %>%
  filter(DX %in% c("SCZ","CTRL")) %>%
  ggplot(., aes(x = SCZ_trs, y = SCZ_prs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "#FF7F0EFF", lwd = 1.1) +
  stat_cor(method = "pearson") +
  facet_wrap(~CT, scales = "free", nrow = 1) +
  theme_classic2()
p2 = metadata_full %>%
  filter(DX %in% c("SCZ","CTRL")) %>%
  ggplot(., aes(x = SCZ_trs, y = Log_density)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "#FF7F0EFF", lwd = 1.1) +
  stat_cor(method = "pearson") +
  facet_wrap(~CT, scales = "free", nrow = 1) +
  theme_classic2()
p3 = metadata_full %>%
  filter(DX %in% c("SCZ","CTRL")) %>%
  ggplot(., aes(x = SCZ_trs, y = MGP)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "#FF7F0EFF", lwd = 1.1) +
  stat_cor(method = "pearson") +
  facet_wrap(~CT, scales = "free", nrow = 1) +
  theme_classic2()
p4 = metadata_full %>%
  filter(DX %in% c("SCZ","CTRL")) %>%
  ggplot(., aes(x = SCZ_trs, y = Intergenic_rate)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "#FF7F0EFF", lwd = 1.1) +
  stat_cor(method = "pearson") +
  facet_wrap(~CT, scales = "free", nrow = 1) +
  theme_classic2()


p1 / p2 / p3 / p4



















