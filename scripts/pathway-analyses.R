
################################################ Functions #########################################################
# Code adapted from: https://github.com/djunamay/apoe4myelin

read.geneset = function(path_to_gset){
  bp = GSA.read.gmt(path_to_gset)
  out = bp$genesets
  out = lapply(1:length(out), function(x) out[[x]][out[[x]]!=''])
  names(out) = bp$geneset.names
  return(out)
}

filter_lowly_exp_genes = function(expressed, all_paths){
  gsets_per_celltype = list()
  
  for(i in names(expressed)){
    index = unlist(lapply(names(all_paths), function(x) (sum(all_paths[[x]]%in%expressed[[i]])/length(all_paths[[x]]))>0.33))
    p = all_paths[(index)]
    x = lapply(names(p), function(y) intersect(expressed[[i]], p[[y]]))
    names(x) = names(p)
    x = x[!duplicated(x)]
    gsets_per_celltype[[i]] = x
  }
  return(gsets_per_celltype)
}

get_gset_names_by_category = function(cat, gsets){
  gset = lapply(gsets, function(x){ 
    y = sapply(cat, grepl, x)
    if(sum(y)>0){
      df = data.frame(parentTerm = names(y[y == TRUE]))
      df$term = x
      df$n = nrow(df)
      return(df)
    } else {
      NULL
    }
  })
  gset = gset[!sapply(gset, is.null)]
  gset = do.call(rbind, gset)
  return(gset)
}

get_pathway_fits = function(order, av_expression, pathways, top_40, summary, gsva_pathway_summarys){
  # run GSVA on GO pathways for each cell type
  all_data = list()
  if(missing(gsva_pathway_summarys)){
    print('running GSVA...')
    out_bp_terms = lapply(order, function(x){ 
      t(gsva(expr = av_expression[[x]], gset.idx.list = pathways[[x]],
             mx.diff = TRUE, verbose = TRUE, 
             kcdf = c("Gaussian"), min.sz = 5, max.sz = 150, parallel.sz = 12))
      })
    names(out_bp_terms) = order
    # split gsva matrix into case-control groups for each disorder 
    # facilitates linear modeling step
    res = lapply(out_bp_terms, function(x){
      l1 = dlply(summary[summary$Group == "CASE" & summary$ID %in% rownames(x),], .(DX), function(y) x[y$ID,])
      l2 = dlply(summary[summary$Group == "CTRL" & summary$ID %in% rownames(x),], .(DX), function(y) rep(list(x[y$ID,]),3))
      l2 = lapply(rapply(l2, enquote, how = "unlist"), eval)
      tmp = lapply(seq_along(l1), function(x) rbind(l1[[x]], l2[[x]]))
      names(tmp) = names(l1)
      return(tmp)
    })
    res = unlist(res, recursive = FALSE)
    
    all_data[['gsva_out']] = res
  } else {
    out_bp_terms = gsva_pathway_summarys
    all_data[['gsva_out']] = out_bp_terms
  }
  
  # get linear model fits
  print('getting linear model fits...')
  fits = get_fits(out_bp_terms, summary)
  all_data[['fits_all']] = fits
  
  # get matrix of scores for heatmap
  print('get matrix of scores')
  scores = get_scores(fits)
  all_data[['scores_all']] = scores
  
  print('filter by score 1.3')
  names = unique(unname(unlist(lapply(names(scores$all), function(x) rownames(scores$all[[x]])))))
  mat = get_matrix(scores$all, names)
  mat = mat[unname(rowSums(abs(mat)>1.3)>0),]

  if(top_40 == TRUE){
    index = c(unique(unname(unlist(lapply(colnames(mat), function(x) order(mat[[x]], decreasing = F)[1:20])))),
              rev(unique(unname(unlist(lapply(colnames(mat), function(x) order(mat[[x]], decreasing = T)[1:20]))))))
    mat = mat[index,]
  }
  all_data[['scores_filtered']] = mat
  return(all_data)
}

get_fits = function(gsva_out, meta){
  fits = list()
  for(i in names(gsva_out)){
    predict = meta[as.character(rownames(gsva_out[[i]])),]
    mod = model.matrix(~Group + Sex + scale(Intergenic_rate), predict)
    fits[[i]] = fit.gsva(mod, i, gsva_out, 'GroupCASE')
  }
  return(fits)
}
fit.gsva = function(mod1, i, gsva.per.celltype, coef){
  fit = lmFit(t(gsva.per.celltype[[i]]), design = mod1)
  fit = eBayes(fit)
  allgenesets = topTable(fit, coef = coef, number = Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
  allgenesets$celltype = i
  allgenesets$names = rownames(allgenesets)
  return(allgenesets)
}
get_scores = function(fits){
  all = list()
  for(i in names(fits)){
    df = fits[[i]]
    df$celltype = i
    df = df[order(abs(df$logFC),decreasing = T),]
    all[[i]] = df[,c('celltype','logFC','P.Value','adj.P.Val','names')]
  }
  return(list('all' = all))
}
get_matrix = function(scores, top_paths){
  df = do.call('rbind',scores)
  df$score = sign(df$logFC) * -log10(df$P.Value)
  df = as.data.frame(pivot_wider(df[,c('celltype', 'names', 'score')], values_from = 'score', names_from = 'celltype'))
  df[is.na(df)] = 0
  rownames(df) = df$names
  df$names = NULL
  return(df[top_paths,])
}

################################################# Script ###########################################################
# load libraries 
suppressPackageStartupMessages({
  library(plyr)
  library(tidyverse)
  library(data.table)
  library(cba)
  library(GSVA)
  library(GSA)
  library(edgeR)
  library(limma)
  library(gprofiler2)
  library(naniar)
  library(parallel)
  library(ComplexHeatmap)
  library(circlize)
  library(rrvgo)
  library(dynamicTreeCut)
  library(seriation)
  library(openxlsx)
  library(ggforce)
  library(ggpubr)
  library(svglite)
  library(scales)
  library(ggsci)
  library(igraph)
  library(ggraph)
  library(tidygraph)
})

# combined exon and intron counts matrix, sample metadata, log2 cpm matrix, gene annotations 
load(file = "./data/lcm-seq/data.Rdata")
# add alignment qc
metadata = metadata %>% left_join(., star_qc, by = 'ID')
rownames(metadata) = metadata$ID
# add case-control group
metadata$Group = as.character(metadata$DX)
metadata$Group[metadata$DX %in% c("MDD","BD","SCZ")] = "CASE"
metadata$Group = factor(metadata$Group, levels = c("CTRL","CASE"))

# load database pathway terms
db = c(read.geneset('./data/pathway-dbs/GO_Biological_Process_2021.txt')
       #read.geneset('./data/pathway-dbs/GO_Cellular_Component_2021.txt'),
       #read.geneset('./data/pathway-dbs/GO_Molecular_Function_2021.txt')
       )
# db = read.geneset('./data/pathway_databases/HumanCyc_2016.txt')
# db = read.geneset('./data/pathway_databases/KEGG_2021_Human.txt')
# db = read.geneset('./data/pathway_databases/Reactome_2022.txt')

# convert to ensembl identifiers 
db = lapply(db, mapvalues, from = ensembl$gene_symbol, to = ensembl$gene_ensembl, warn_missing = FALSE)

# calculate scaling factors to convert raw library sizes into effective library sizes
dge0 = DGEList(counts_matrix)
dge0 = calcNormFactors(dge0, method = "TMM")

# get filtered genes by cell type
expressed = dlply(metadata, .(CT), function(x){
  keep = rowSums(cpm(dge0[,x$ID]) > 1) >= ncol(dge0[,x$CT])*0.30
  row.names(dge0$counts[keep,])
})
# get filtered expression by cell type 
av_expression = dlply(metadata, .(CT), function(x){
  keep = rowSums(cpm(dge0[,x$ID]) > 1) >= ncol(dge0[,x$CT])*0.30
  voom(dge0[keep, x$ID], plot = T)$E
})
# get gene filtered GO pathways 
pathways = filter_lowly_exp_genes(expressed, db)
pathways = pathways[order(match(names(pathways), as.character(unique(metadata$CT))))]

# get gene filtered disease-associated pathways
disease_keys = c("transcription","translation","mitochondria","ATP","respira","metabol","histone","MAPK","acetylation","endoplasmic",
                 "ER","Golgi","vesicle","lysosom","vacuol","organelle","lumen","cytoskel","actin","matrix","axon","neuro","dendrit","synap","channel",
                 "receptor","growth","hormone","angio","inflamma","immune","cytokine","chemokine","adaptive","innate")

disease_pathways_df = get_gset_names_by_category(disease_keys, names(db))
disease_pathways = filter_lowly_exp_genes(expressed, db[disease_pathways_df$term])

gsva_out = readRDS(file = "./output/gsva_out_go_disease_all.rds")
order = names(av_expression)
pathway_scores = get_pathway_fits(order, av_expression, disease_pathways, top_40 = TRUE, metadata, gsva_out)

#gsva_out = pathway_scores$gsva_out
#saveRDS(gsva_out, file = "./output/gsva_out_go_disease_bp.rds")

################################################## Figures #########################################################

## Heatmap ----

# select pathways significantly associated with DX in 1 or more CT*DX pairs
mat = pathway_scores[['scores_filtered']]
mat = mat[unname(rowSums(abs(mat) > -log10(0.01)) >= 1),]
mat = mat[!grepl("\\.1$", row.names(mat)),]
mat = mat[,c("PVALB.SCZ","PVALB.BD","PVALB.MDD","SST.SCZ","SST.BD","SST.MDD","VIP.SCZ","VIP.BD","VIP.MDD",
             "PyrL2n3.SCZ","PyrL2n3.BD","PyrL2n3.MDD","PyrL5n6.SCZ","PyrL5n6.BD","PyrL5n6.MDD")]

# # add biological themes
# anno_terms_full = disease_pathways_df %>%
#   filter(term %in% rownames(mat)) %>%
#   mutate(group = case_match(parentTerm,
#                             c("translation","transcription") ~ "Transcription/Translation",
#                             c("mitochondria","ATP","respira") ~ "Energetics",
#                             c("metabol") ~ "Metabolic/Macromolecular processes",
#                             c("histone","MAPK","acetylation") ~ "Protein modification",
#                             c("endoplasmic","ER","Golgi","vesicle","lysosom","vacuol","organelle","lumen") ~ "Cellular trafficking",
#                             c("cytoskel","actin","matrix") ~ "Cytoskeleton",
#                             c("axon","neuro","dendrit","synap","channel","receptor","transport") ~ "Neuron Structure/Function",
#                             c("growth","hormone") ~ "Growth factors",
#                             c("angio") ~ "Angiogenesis",
#                             c("inflamma","immune","cytokine","chemokine","adaptive","innate") ~ "Immune process"))

# read manually-adjusted biological themes 
anno_terms_full = read.csv(file = "./output/go_groups.csv") %>% filter(!term == "")
anno_terms_full = anno_terms_full %>% mutate(group = factor(group, levels = unique(anno_terms_full$group))) 
anno_terms = anno_terms_full %>% dplyr::select(term, group)

# cluster terms within each theme
group_lst = split(anno_terms$term, anno_terms$group)
mat_reorder = lapply(group_lst, function(x){
  
  tmp = as.matrix(mat[rownames(mat) %in% x,])
  tmp[seriate(tmp, margin = 1)[[1]],]
})
# gaps for biological themes 
gaps = cumsum(sapply(mat_reorder, nrow))

mat_reorder_full = do.call(rbind, mat_reorder)
mat_reorder = mat_reorder_full
anno_terms = anno_terms[match(rownames(mat_reorder), anno_terms$term),]
anno_terms_full = anno_terms_full[match(rownames(mat_reorder), anno_terms_full$term),]
rownames(anno_terms) = NULL
anno_terms = column_to_rownames(anno_terms, var = "term")
rownames(anno_terms) = gsub("\\s*\\([^\\)]+\\)", "", rownames(anno_terms)) # remove go ids 
rownames(mat_reorder) = gsub("\\s*\\([^\\)]+\\)", "", rownames(mat_reorder))

# annotation colors
n = length(unique(anno_terms$group)) 
anno_colors = ggsci::pal_d3("category20")(n)
names(anno_colors) = unique(anno_terms$group)
anno_colors_list = list(group = anno_colors)

# color scale
my_breaks = seq(-2,2,0.05)
scale_colors = c(rev(colorRampPalette(c("#000000", "#1F77B4FF"))(length(my_breaks)/2)), colorRampPalette(c("#000000", "#D62728FF"))(length(my_breaks)/2)) 
wpos = which(scale_colors %in% "#000000")
scale_colors = append(scale_colors, rep("#000000", 16), wpos[1]) # more white

# draw heatmap for GO pathways
#svglite("./figures/draft/gobp_heatmap.svg", height = 11, width = 9)
pheatmap(as.matrix(mat_reorder),
         color = scale_colors,
         breaks = my_breaks,
         display_numbers = F,
         border_color = NA,
         cluster_rows = F, cluster_cols = F,
         show_rownames = T,
         row_names_side = "left",
         fontsize_row = 2.8,
         #labels_row = path_names$selected_name
         gaps_col = c(3,6,9,12),
         gaps_row = gaps,
         cellheight = 2.2, cellwidth = 10,
         annotation_row = anno_terms,
         annotation_colors = anno_colors_list
)
#dev.off()

out_xl = mat_reorder_full %>%
  as.data.frame() %>%
  rownames_to_column(var = "term") %>%
  left_join(., anno_terms_full, by = "term")
#write.csv(out_xl, file = "./output/supplementary_table_2.csv")

## Network ----

lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds")
disorders = as.character(unique(lcm_df_limma$DX))
ct_pairs = combn(as.character(unique(lcm_df_limma$CT)), 2, simplify = FALSE)
names(group_lst) = sub("/", " & ", names(group_lst))

# for each biological theme
cor_df = data.frame(do.call(rbind, lapply(names(group_lst), function(x){
  # for each disorder
  res_d = data.frame(do.call(rbind, lapply(disorders, function(dx){
    # for each cell type pair 
    df = lcm_df_limma %>% filter(gene_ensembl %in% unique(unlist(db[group_lst[[x]]])))
    res_c = lapply(ct_pairs, function(ct){
      v1 = df %>% filter(DX == dx, CT == ct[1]) %>% pull(t, gene_ensembl)
      v2 = df %>% filter(DX == dx, CT == ct[2]) %>% pull(t, gene_ensembl)
      v1 = v1[intersect(names(v1), names(v2))]
      v2 = v2[intersect(names(v1), names(v2))]
      data.frame(cor = cor(v1,v2, method = "pearson"), pval = cor.test(v1,v2)$p.value)
    })
    res_c = data.frame(do.call(rbind, res_c)) %>%
      mutate(adj_pval = p.adjust(pval, method = "fdr"), 
             group = as.character(lapply(ct_pairs, function(x) paste(x[1],x[2]))),
             DX = as.character(dx),
             theme = as.character(x)) %>%
      separate(group, into = c("c1", "c2"), sep = " ")
    return(res_c)
        
  #   df = lcm_df_limma %>%
  #     filter(DX == dx, gene_ensembl %in% unique(unlist(db[group_lst[[x]]]))) %>%
  #     dplyr::select(CT, gene_ensembl, t) %>%
  #     pivot_wider(id_cols = gene_ensembl, names_from = CT, values_from = t, values_fn = mean) %>%
  #     column_to_rownames(var = "gene_ensembl")
  #   
  #   png(filename = paste0("./figures/draft/go_networks/",x,"_",dx,".png"), width = 8, height = 8, units = "in", res = 300)
  #   PerformanceAnalytics::chart.Correlation(df, histogram = F, pch = 19)
  #   mtext(paste(x, dx), side = 3, line = 3, cex = 1.5)
  #   dev.off()
  })))
  return(res_d)
})))

plt_lst = list()
for (x in names(group_lst)) {
  res_d = cor_df %>% filter(theme == x)
  
  graph = res_d[c("c1","c2","pval","DX")]
  graph$pval[graph$pval > 0.05] = NA
  graph$pval[graph$cor < 0] = NA
  graph = graph.data.frame(graph, directed = F)
  g = as_tbl_graph(graph)
  
  layout = as.matrix(data.frame(x = c(0, -1, 1, -0.5, 0.5), y = c(-sqrt(4), -sqrt(3)/2, -sqrt(3)/2, sqrt(3)/2, sqrt(3)/2)))
  
  plt_lst[[x]] = ggraph(g, layout = layout) +
      geom_edge_parallel(aes(width = -log10(pval), color = DX), show.legend = TRUE) +
      scale_edge_width(range = c(0.9, 2.2)) + 
      geom_node_point(size = 25, shape = 16, color = "grey") +
      geom_node_text(aes(label = name), family = "sans", vjust = 0.5, hjust = 0.5) +
      scale_edge_color_manual(values = c("SCZ" = "#393B79FF", "BD" = "#7B4173FF", "MDD" = "#AD4941FF")) +
      theme_void() +
      theme(plot.margin = margin(1, 1, 1, 1, "cm"),
            plot.title = element_text(hjust = 0.5, vjust = 10, size = 16),
            legend.position = "none",
            legend.margin = margin(1, 1, 1, 1, "cm"),
            legend.box.margin = margin(0, 0, 0, 0)) +
      coord_cartesian(clip = "off") +
      ggtitle(paste(x)) 
}
svglite("./figures/draft/go_networks/combined.svg", height = 8, width = 16)
do.call(gridExtra::grid.arrange, c(plt_lst, ncol = 4))
dev.off()

## Dotplots ----

go_summary = cbind(anno_terms, mat_reorder) %>%
  rownames_to_column(var = "term") %>%
  pivot_longer(cols = -c(term, group), values_to = "score") %>%
  mutate(dir = sign(score)) %>%
  separate(name, c("CT","DX")) %>%
  mutate(sig = abs(score) > -log10(0.05),
         term = factor(term, levels = unique(term)),
         DX = factor(DX, levels = c("SCZ","BD","MDD"), labels = c("SCZ","BD","MDD")),
         CT = factor(CT, levels =  c("PVALB","SST","VIP","PyrL2n3","PyrL5n6"), labels = c("PVALB","SST","VIP","L2/3 PYR","L5/6 PYR")))

select_terms = c("positive regulation of endoplasmic reticulum unfolded protein response",
                 "positive regulation of ATPase activity",
                 "regulation of NLRP3 inflammasome complex assembly",
                 "positive regulation of lipid metabolic process",
                 "glial cell-derived neurotrophic factor receptor signaling pathway",
                 "regulation of glucocorticoid receptor signaling pathway",
                 "glutamine metabolic process",
                 "tetrahydrobiopterin metabolic process",
                 "dendritic spine maintenance",
                 "negative regulation of translation")

select_terms = c("positive regulation of ATPase activity",
                 "negative regulation of cellular amide metabolic process",
                 "negative regulation of translation",
                 "glial cell-derived neurotrophic factor receptor signaling pathway",
                 "positive regulation of response to endoplasmic reticulum stress",
                 "negative regulation of macrophage cytokine production",
                 "cell-matrix adhesion mediator activity",
                 "regulation of retinoic acid receptor signaling pathway", 
                 "long-term synaptic depression",
                 "neuropeptide hormone activity",
                 "lysosomal protein catabolic process")

plot_df = go_summary %>% filter(term %in% select_terms) 
ggplot() + 
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(log10(0.05), -log10(0.05)), linetype = "dashed") +
  geom_point(data = plot_df %>% filter(sig == FALSE), aes(x = score, y = CT, color = group, shape = DX), alpha = 0.6, size = 3) +
  geom_point(data = plot_df %>% filter(sig == TRUE), aes(x = score, y = CT, color = group, shape = DX), alpha = 1, size = 3) +
  scale_color_manual(values = anno_colors_list$group) + 
  scale_y_discrete(limits = rev) + 
  labs(x = "signed log(P-value)", y = "") +
  guides(color = guide_legend(nrow = 4, title.position = "top", override.aes = list(shape = 15, size = 6)),
         shape = guide_legend(nrow = 1, title.position = "top")) +
  facet_wrap(~ reorder(term, group), ncol = 1) + 
  theme_classic2() +
  theme(strip.text = element_text(face = "bold"),
        title = element_text(face = "bold", size = 12), 
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(face = "bold"), 
        axis.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(),
        legend.position = "bottom")

ggsave(filename = "./figures/draft/tmp.svg", width = 5, height = 15)

go_top = cbind(anno_terms, mat_reorder) %>%
  rownames_to_column(var = "term") %>%
  pivot_longer(cols = -c(term, group), values_to = "score") %>%
  #slice_max(order_by = abs(score), n = 250) %>%
  pull(term) %>% unique()

plot_df = go_summary %>% filter(term %in% go_top) 
p = ggplot() + 
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(log10(0.05), -log10(0.05)), linetype = "dashed") +
  geom_point(data = plot_df %>% filter(sig == FALSE), aes(x = score, y = CT, color = group, shape = DX), alpha = 0.6, size = 3) +
  geom_point(data = plot_df %>% filter(sig == TRUE), aes(x = score, y = CT, color = group, shape = DX), alpha = 1, size = 3) +
  scale_color_manual(values = anno_colors_list$group) + 
  scale_y_discrete(limits = rev) + 
  labs(x = "signed log(P-value)", y = "") +
  guides(color = guide_legend(override.aes = list(size = 6)), size = "none") +
  theme_classic2() +
  theme(strip.text = element_text(face = "bold"),
        title = element_text(face = "bold", size = 12), 
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(face = "bold"), 
        axis.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(),
        legend.position = "right") +
  facet_wrap_paginate(~ term, ncol = 2, nrow = 7, labeller = label_wrap_gen(width = 50)) 
  
for(i in 1:n_pages(p)){
  p + facet_wrap_paginate(~ term, ncol = 2, nrow = 7, labeller = label_wrap_gen(width = 50), page = i) 
  ggsave(filename = paste0("./figures/draft/go_terms/page_",i,".png"), width = 12, height = 12)
}



################################################## Tables #########################################################

# summarize pathway pathway_summarys
pathway_summary = cbind(anno_terms_full[,1:2], mat_reorder) %>%
  pivot_longer(cols = -c(term, group), values_to = "score", names_to = "x") %>%
  separate(x, c("CT","DX")) %>%
  mutate(dir = sign(score), DE = abs(score) > -log10(0.01))
  
# which DEGs drive pathway changes 
lcm_df_limma = read_rds(file = "./output/lcm_df_limma.rds") 
genes = lcm_df_limma %>% pull(gene_ensembl) %>% unique()
  
tmp = lapply(db, function(x) intersect(x, genes))
tmp = tmp[sapply(tmp, length) > 0]
pathway_degs = data.frame("term" = names(tmp), 
                          "gene_ensembl" = sapply(tmp, function(x) paste(x, collapse = ", ")), 
                          stringsAsFactors = FALSE) %>%
  separate_rows("gene_ensembl", sep = ", ") %>%
  right_join(., pathway_summary %>% filter(DE == TRUE), "term") %>%
  dplyr::select(-DE) %>%
  left_join(., lcm_df_limma %>% filter(DE == TRUE), by = c("CT","DX","gene_ensembl")) %>%
  filter(!gene_symbol == "")
  

# dot_plot = function(go_summary, select_term){
#   
#   colors = c("SCZ" = "black", "BD" = "darkgrey", "MDD" = "lightgrey")
#   go_summary %>%
#     filter(term == select_term) %>%
#     ggplot(., aes(x = score, y = CT, color = DX, size = abs(score))) + 
#     geom_vline(xintercept = 0) +
#     geom_vline(xintercept = c(log10(0.05), -log10(0.05)), linetype = "dashed") +
#     geom_point(alpha = 0.95, shape = 16) +
#     scale_size(range = c(4, 9)) +
#     scale_color_manual(values = colors) + 
#     #scale_x_continuous(limits = c(-2.7,2.7)) + 
#     scale_y_discrete(limits = rev) + 
#     labs(title = gsub("\\([^()]*\\)", "", select_term), x = "signed log(P-value)", y = "") +
#     guides(color = guide_legend(override.aes = list(size = 6)), size = "none") +
#     theme_classic2() +
#     theme(strip.text = element_text(face = "bold"),
#           title = element_text(face = "bold", size = 12), 
#           axis.text = element_text(color = "black"),
#           axis.text.y = element_text(face = "bold"), 
#           axis.title = element_text(face = "bold"),
#           legend.position = "none")
# }
#   
# p1 = dot_plot(go_summary, "positive regulation of endoplasmic reticulum unfolded protein response")
# p2 = dot_plot(go_summary, "positive regulation of ATPase activity")
# p3 = dot_plot(go_summary, "regulation of NLRP3 inflammasome complex assembly")
# p4 = dot_plot(go_summary, "positive regulation of lipid metabolic process")
# p5 = dot_plot(go_summary, "glial cell-derived neurotrophic factor receptor signaling pathway")
# p6 = dot_plot(go_summary, "regulation of glucocorticoid receptor signaling pathway")
# p7 = dot_plot(go_summary, "glutamine metabolic process")
# p8 = dot_plot(go_summary, "tetrahydrobiopterin metabolic process")
# p9 = dot_plot(go_summary, "dendritic spine maintenance")
# p10 = dot_plot(go_summary, "negative regulation of translation")
# 
# svglite("./figures/draft/go_dotplots.svg", height = 11, width = 9)
# ggarrange(plotlist = list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10), ncol = 2, nrow = 5)
# dev.off()

  
# # pie chart
# prop = table(anno_terms$parentTerm)
# prop = prop[match(unique(anno_terms$parentTerm)[-1],names(prop))]
# svglite("./figures/draft/gobp_pie.svg", height = 3, width = 3)
# pie(prop, labels = "", col = anno_colors, border = "grey")
# dev.off()

# # clean-up pathway names
# path_names = pathway_scores$fits_all %>%
#   bind_rows() %>% as.data.frame() %>%
#   filter(names %in% rownames(mat)) %>%
#   mutate(Score = sign(logFC)*-log10(P.Value)) %>%
#   group_by(names) %>%
#   top_n(1, abs(Score)) %>%
#   dplyr::select(names, celltype, logFC, P.Value, adj.P.Val, Score)
# path_names = path_names[match(rownames(mat), path_names$names),]
# #write.csv(path_names, file = "./output/path_names.csv", row.names = F)
# 
# path_names = read.csv(file = "./output/path_names.csv")
# path_names[which(is.na(path_names$names)), 1] = rownames(mat)[which(is.na(path_names$names))]
# rownames(path_names) = path_names$names


# # parent terms from semantic similarity 
# re = "\\(([^()]+)\\)"
# go_term_ids = gsub(re, "\\1", str_extract_all(rownames(mat), re))
# sim_matrix = calculateSimMatrix(go_term_ids, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
# reduced_terms = reduceSimMatrix(sim_matrix, scores = NULL, threshold = 0.7, orgdb = "org.Hs.eg.db")
# unique(reduced_terms$parentTerm)
# 
# heatmapPlot(sim_matrix, reduced_terms, annotateParent = TRUE, annotationLabel = "parentTerm", fontsize = 6)
# # treemapPlot(reduced_terms, size = "score") 
# # scatterPlot(sim_matrix, reduced_terms)
# 
# anno_terms = reduced_terms %>%
#   mutate(term = paste0(reduced_terms$term," (",reduced_terms$go,")")) %>%
#   filter(term %in% rownames(mat)) %>%
#   group_by(parentTerm) %>%
#   mutate(n = n()) 
# anno_terms = anno_terms[match(rownames(mat), anno_terms$term),]
# anno_terms[which(is.na(anno_terms$term)), "term"] = rownames(mat)[which(is.na(anno_terms$term))]
# anno_terms = anno_terms %>%
#   arrange(desc(n), parentTerm)
# anno_terms = anno_terms[,c("term","parentTerm")]
# rownames(anno_terms) = NULL
# anno_terms = column_to_rownames(anno_terms, var = "term")
# 
# mat_reorder = mat[match(rownames(anno_terms), rownames(mat)),]
# 
# # annotation colors
# n = length(unique(anno_terms$parentTerm)) - 1
# anno_colors = ggsci::pal_d3("category20")(n)
# names(anno_colors) = unique(anno_terms$parentTerm)[-1]
# anno_colors_list = list(parentTerm = anno_colors)



# selected_parent_terms = c("calcium-mediated signaling using intracellular calcium source","dendritic spine maintenance","dopamine uptake","glutamine metabolic process",
#                           "integrin activation","lysosome organization","monoacylglycerol catabolic process","negative regulation of amyloid precursor protein biosynthetic process",
#                           "negative regulation of CD4-positive, alpha-beta T cell activation","negative regulation of macrophage cytokine production",
#                           "negative regulation of tumor necrosis factor-mediated signaling pathway","noradrenergic neuron differentiation","positive regulation of histone acetylation",
#                           "protein deacetylation","protein localization to chromatin","regulation of early endosome to late endosome transport")
# anno_terms = reduced_terms %>% 
#   replace_with_na_at(.vars = "parentTerm", condition = ~!.x %in% selected_parent_terms) %>%
#   mutate(parentTerm = dplyr::recode(parentTerm, 
#                                     "calcium-mediated signaling using intracellular calcium source" = "calcium-mediated signaling",
#                                     "monoacylglycerol catabolic process" = "monoacylglycerol metabolism",
#                                     "negative regulation of amyloid precursor protein biosynthetic process" = "amyloid biosynthesis",
#                                     "negative regulation of CD4-positive, alpha-beta T cell activation" = "immune process",
#                                     "negative regulation of macrophage cytokine production" = "immune process",
#                                     "negative regulation of tumor necrosis factor-mediated signaling pathway" = "immune process",
#                                     "noradrenergic neuron differentiation" = "neuron differentiation",
#                                     "positive regulation of histone acetylation" = "chromatin remodeling",
#                                     "protein deacetylation" = "chromatin remodeling",
#                                     "regulation of early endosome to late endosome transport" = "endosome transport")) %>%
#   mutate(term = paste0(reduced_terms$term," (",reduced_terms$go,")")) %>%
#   filter(term %in% rownames(mat)) %>%
#   dplyr::select(term, parentTerm) 
# anno_terms = anno_terms[match(rownames(mat), anno_terms$term),]
# anno_terms[which(is.na(anno_terms$term)), 1] = rownames(mat)[which(is.na(anno_terms$term))]
# rownames(anno_terms) = NULL
# anno_terms = column_to_rownames(anno_terms, var = "term")


