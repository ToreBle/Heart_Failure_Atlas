library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(multinichenetr)

organism = "human"
if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()}

scRNA <- readRDS("./multinichenet/Heart_atlas_reduced.rds")

sce <- as.SingleCellExperiment(scRNA)
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

SummarizedExperiment::colData(sce)$Donor = SummarizedExperiment::colData(sce)$Donor %>% make.names()
SummarizedExperiment::colData(sce)$Condition = SummarizedExperiment::colData(sce)$Condition %>% make.names()
SummarizedExperiment::colData(sce)$Annotation_2 = SummarizedExperiment::colData(sce)$Annotation_2 %>% make.names()
sample_id = "Donor"
group_id = "Condition"
celltype_id = "Annotation_2"
covariates = NA
batches = NA #Pseudobulking of Samples

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()

min_cells = 10

##current state

abundance_expression_info = get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network, batches = batches)

pdf("MultiNicheNet_Annotation_2_Evaluation.pdf")
abundance_expression_info$abund_plot_sample
dev.off()

#contrasts_oi = c("'M-(S+A)/2','S-(M+A)/2','A-(S+M)/2'")
contrasts_oi = c("'HCM-Control','DCM-Control', 'Ischemic-Control', 'Post.MI-Control'")

#contrast_tbl = tibble(contrast = 
#                        c('ACM-Control','HCM-Control','Fibrotic-Control', 'Post.MI-Control'), 
#                      group = c("ACM","HCM", "Fibrotic", "Post.MI"))
contrast_tbl = tibble(contrast = 
                        c('DCM-Control','HCM-Control','Ischemic-Control', 'Post.MI-Control'), 
                      group = c("DCM","HCM", "Ischemic", "Post.MI"))



DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)

DE_info$celltype_de$de_output_tidy %>% arrange(p_adj) %>% head()

celltype_de = DE_info$celltype_de$de_output_tidy


sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)


sender_receiver_de %>% head(20)

print("Step 3")
logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05


p_val_adj = FALSE
top_n_target = 250
verbose = TRUE
cores_system = 4
n.cores = min(cores_system, sender_receiver_de$receiver %>% unique() %>% length()) # use one core per receiver cell type

print("check1")


ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose, 
  n.cores = n.cores
)))

ligand_activities_targets_DEgenes$de_genes_df %>% head(20)

print("check2")

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

print("Step 4")

prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)

prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
))

prioritization_tables$group_prioritization_tbl %>% head(20)

print("Step 5")

lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)

path = "./multinichenet/"

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output_Annotation_2.rds"))
  
}

contrast_tbl = tibble(contrast = 
                        c('DCM-Control','HCM-Control','Ischemic-Control', 'Post.MI-Control'), 
                      group = c("DCM","HCM", "Ischemic", "Post.MI"))

pdf(paste0(path, "MultiNicheNet_Annotation_2.pdf"), width = 30, height = 20)
prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, rank_per_group = TRUE)

prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = colorRampPalette(brewer.pal(11, "Spectral"))(length(senders_receivers)) %>% magrittr::set_names(senders_receivers)
colors_receiver = colorRampPalette(brewer.pal(11, "Spectral"))(length(senders_receivers)) %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
print(circos_list)
dev.off()
