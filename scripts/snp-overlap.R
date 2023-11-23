# snp-overlap.R 
# load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  #(hypeR)
})

# load data ----
load(file = "./data/lcm-seq/data.Rdata")
lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds") 

# read in psyops score for latest MDD GWAS
# https://doi.org/10.1038/s41591-023-02352-1
psyops_df_mdd = read_tsv("./data/wainberg/Als_MDD_2023_out.tsv") %>%
  as.data.frame() %>%
  dplyr::rename(gene_symbol = gene) %>%
  mutate(DX = "MDD")
pops_df_mdd = read.csv(file = "./data/wainberg/PoPS_res_2023_MDD.csv")[1:2] %>%
  dplyr::rename(gene_ensembl = ENSGID, PoPS_score = PoPS_Score) %>%
  mutate(DX = "MDD")
# read in other psyops scores (from Michael)
# https://doi.org/10.1038/s41588-021-00857-4
# https://doi.org/10.1038/s41586-022-04434-5
psyops_df = read.csv(file = "./data/wainberg/psyops_table.csv") %>%
  as.data.frame() %>%
  dplyr::rename(gene_symbol = gene, DX = trait) %>%
  mutate(DX = case_match(DX, "BIP" ~ "BD", .default = DX)) %>%
  filter(!DX == "MDD") %>%
  bind_rows(., psyops_df_mdd) %>%
  left_join(., lcm_df_limma, by = c("DX","gene_symbol"), multiple = "all", relationship = "many-to-many") %>%
  left_join(., pops_df_mdd, by = c("DX","gene_ensembl")) %>%
  mutate(PoPS_score = coalesce(PoPS_score.x, PoPS_score.y)) %>%
  group_by(DX, lead_variant) %>%
  mutate(PsyOPS_highest = PsyOPS_score == max(PsyOPS_score, na.rm = TRUE),
         PoPS_highest = PoPS_score == max(PoPS_score, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(DX, CT, gene_symbol, lead_variant) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  mutate(logFDR = sign(logFC)*-log10(adj.P.Val)) 
#  left_join(., bulk_rnaseq, by = c("DX","gene_ensembl")) 

genes_select = psyops_df %>% filter(DE == TRUE) %>% pull(gene_symbol, DX)
genes_select = split(genes_select, names(genes_select))

psyops_df_select = do.call(rbind, lapply(seq_along(genes_select), function(i){
  psyops_df %>% filter(DX == names(genes_select)[[i]], gene_symbol %in% genes_select[[i]])
})) %>%
  dplyr::select(DX, CT, lead_variant, gene_symbol, distance, PoPS_highest, PsyOPS_highest, logFDR) %>%
  pivot_wider(names_from = CT, values_from = logFDR) %>%
  filter(PoPS_highest == TRUE | PsyOPS_highest == TRUE) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>% 
  arrange(factor(DX, levels = c("SCZ","BD","MDD"))) %>%
  as.data.frame()
names(psyops_df_select) = c("DX","GWAS SNP","Gene","Distance from SNP","PoPS","PsyOPS","PVALB","SST","VIP","PYR L2/3","PYR L5/6")

write.csv(psyops_df_select, file = "./output/Figure-3-v3.csv", row.names = FALSE)

permute_overlap = function(bkgd, l1, l2) {
  rand_list = sample(bkgd, length(l1), replace = FALSE)
  length(intersect(rand_list, l2))
}
n_perm = 10000
disorders = as.character(unique(metadata$DX))[-1]
cell_types = as.character(unique(metadata$CT))

tmp = do.call(rbind,
              lapply(cell_types, function(ct){
                do.call(rbind,
                        sapply(disorders, simplify = F, function(dx){
                          
                          bkgd = lcm_df_limma %>% pull(gene_symbol) %>% unique()
                          l1 = lcm_df_limma %>% filter(DE == TRUE, CT == ct, DX == dx) %>% pull(gene_symbol) %>% unique()
                          l2 = psyops_df %>% filter(DX == dx) %>% pull(gene_symbol) %>% unique()
                          
                          obs = length(intersect(l1, l2))
                          res = replicate(n_perm, permute_overlap(bkgd, l1, l2))
                          
                          data.frame(CT = ct, DX = dx, n_degs = length(l1), n_gwas = length(l2), overlap = obs, 
                                     p_val = sum(res >= obs) / n_perm)}))
                })
)





