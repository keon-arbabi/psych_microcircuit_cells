
# prepare-data-lcm-rna-seq.R
# load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(biomaRt)
  library(edgeR)
})

# prepare metadata ----
# load metadata and format 
metadata = read.csv(file = "data/lcm-seq/metadata.csv", row.names = 1) %>%
  # define factors and numeric variables explicitly
  mutate_at(vars("Age","PMI","pH","RNA.Ratio","RIN"), as.numeric) %>%
  # add scaled numeric variables
  mutate_at(vars("HU","Tetrad","MOD","Sex","Race","Tob.ATOD","Hand.Preference","Pool"), factor) %>%
  # set controls as reference level
  mutate(DX = factor(DX, levels = c("CTRL","MDD","BD","SCZ")),
  # standard order for cell types
         CT = factor(CT, levels =  c("PVALB","SST","VIP","PyrL2n3","PyrL5n6")))
  
# prepare counts matrix ----
# load processed introns and exons data 
introns = fread(file = "./data/lcm-seq/raw_counts/introns.csv", data.table = F)[,-1] %>% column_to_rownames(var = "Gene_ID")
exons = fread(file = "./data/lcm-seq/raw_counts/exons.csv", data.table = F)[,-1] %>% column_to_rownames(var = "Gene_ID")
# combine
counts_matrix = introns + exons
# sum technical replicates (will take a while)
counts_matrix = counts_matrix %>% t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") %>%
  separate(col = id, c("sample_id"), "_", remove = T) %>%
  group_by(sample_id) %>%
  summarize_each(funs(sum)) %>%
  column_to_rownames(var = "sample_id") %>%
  t() %>% as.data.frame()
names(counts_matrix) = gsub("\\.","-", names(counts_matrix))
# fix sample naming issue 
names(counts_matrix)[names(counts_matrix) %in% c("PV-Hu1195","PV-Hu1222")] = c("PVALB-Hu1195", "PVALB-Hu1222")
# match and check 
counts_matrix = counts_matrix[, match(rownames(metadata), colnames(counts_matrix))]
all(rownames(metadata) == colnames(counts_matrix)) 
#write.csv(counts_matrix, file = "./data/lcm-seq/counts_matrix.csv")

# log cpm matrix
lcpm_matrix = cpm(counts_matrix, prior.count = 2, log = TRUE) 

# get gene symbols and keep separately
# mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# saveRDS(mart, file = "./output/ensembl_mart.rds")
mart = readRDS(file = "./output/ensembl_mart.rds")
ensembl = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id","description","gene_biotype","percentage_gene_gc_content"), mart = mart) 
names(ensembl) = c("gene_ensembl","gene_symbol","gene_entrez","gene_description","gene_biotype","gene_percentage_gc_content")

# remove unwanted genes 
rm_genes = ensembl %>% dplyr::filter(str_detect(gene_symbol, "^LOC|^MT-|^RP[0-9]|^BC[0-9]|-PS")) %>% pull(gene_ensembl)
counts_matrix = counts_matrix[!row.names(counts_matrix) %in% rm_genes,]

# prepare star quality control metrics ----
star_qc_raw = read.csv(file = "./data/lcm-seq/STAR_qc.csv", row.names = 1) %>%
  separate(col = sample, "ID", "_", remove = T) %>%
  mutate(ID = gsub("\\.","-", ID)) 
# collapse across technical replicates, by 
# summing reads
star_qc = star_qc_raw %>%
  dplyr::select(-c(Exon_only_genes, Intron_only_genes, All_detected_genes, Unique_align_pct, Unmapped_pct, Pool)) %>%
  group_by(ID) %>%
  summarize_all(funs(sum)) %>%
  mutate(Intronic_rate = Intronic_reads/Uniquely_aligned_reads,
         Intragenic_rate = Genic_reads/Uniquely_aligned_reads,
         Intergenic_rate = Intergenic_reads/Uniquely_aligned_reads,
         Mapped_pct = ((Genic_reads + Intergenic_reads + Multi_mapped_reads)/Sequenced_reads)*100,
         Uniquely_mapped_pct = (Uniquely_aligned_reads/Sequenced_reads)*100,
         Multi_mapped_pct = (Multi_mapped_reads/Sequenced_reads)*100,
         Unmapped_pct = (Unmapped_reads/Sequenced_reads)*100
  ) %>%
  left_join(., data.frame(ID = colnames(counts_matrix), Detected_genes = colSums(counts_matrix > 0)), by = "ID") %>%
  left_join(., data.frame(ID = rownames(DGEList(counts_matrix)$samples), Library_size = DGEList(counts_matrix)$samples$lib.size))
rownames(metadata) = metadata$ID

# save all
#save(metadata, star_qc, counts_matrix, lcpm_matrix, ensembl, file = "./data/lcm-seq/data.Rdata")

# plots ----
# get ID order by unmapped reads
my_order = star_qc %>%
  arrange(Unmapped_reads)
# reformatting for plots
star_qc_long = star_qc %>%
  dplyr::select(ID, Exonic_reads, Intronic_reads, Intergenic_reads, Multi_mapped_reads, Unmapped_reads) %>%
  pivot_longer(cols = -ID, names_to = "Read_type", values_to = "Read_amount") %>%
  mutate(Read_type = gsub("_reads", "", Read_type)) %>%
  mutate(Read_type = factor(Read_type)) %>%
  left_join(., metadata %>% dplyr::select(c("ID","DX","CT")), by = "ID")

# alignment percent breakdown
star_qc_long %>%
  mutate(ID = factor(ID, levels = my_order$ID)) %>%
  ggplot() +
  geom_bar(aes(y = Read_amount, x = ID, fill = Read_type), stat = "identity", position = 'fill', width = 1.2) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_viridis_d() +
  labs(y = "Read type percentage", x = "Samples", fill = "Type of read") +
  facet_nested_wrap(~CT + DX, nrow = 1, scales = "free_x") +
  theme_classic2() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing = unit(0, 'lines'), panel.border = element_rect(colour = "black", fill = NA, size = 1))
ggsave(filename = "./figures/draft/star_breakdown.svg", width = 10, height = 4)
