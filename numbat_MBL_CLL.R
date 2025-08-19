###################################################################################
#                                                                                 #
# -------------------------------CVA_RNA_MBL_CLL----------------------------------#
#                                                                                 #
###################################################################################

# Goal: Identify allele-specific CNVs from scRNA data with single cell resolution #

###########################
#Load libraries
###########################
library(Seurat)
library(numbat)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

###########################
#Set paths
###########################
source("config_paths.R")
data_dir <- file.path(project_dir, "data")
result_dir <- file.path(project_dir_ext, "numbat")
GEX_analysis_dir <- file.path(project_dir, "GEX_analysis")
GEX_RDS_dir <- file.path(GEX_analysis_dir, "RDS_objects")
GEX_plotting_dir <- file.path(GEX_analysis_dir, "Plots")

###########################
#Load data
###########################
source(file.path(project_dir, "scripts/publication_theme.R"))

metadata <- read_tsv(file = file.path(data.dir, "20240723_Tidy_MBL_CLL_metadata_gex.tsv"))

label2pool <- metadata %>%
  dplyr::select(label, pool) %>%
  deframe()

label2plottingid <- metadata %>%
  dplyr::select(label, plotting_id) %>%
  deframe()

gex_scaled_all <- readRDS(file.path(GEX_RDS_dir, "gex_scaled_all.RDS"))

allele_lst <- metadata$label %>%
  purrr::map(~ read_tsv(paste0(result_dir, "preprocessing/allowed_ALTREFswap/", label2pool[.], "/", label2pool[.], "_allele_counts.tsv.gz"))) %>%
  purrr::set_names(metadata$label)

allele_lst <- imap(allele_lst, ~ dplyr::mutate(.x, label = .y)) %>%
  map(~ tidyr::unite(., "cell", c(label, cell), sep = "_"))

###########################
#Run numbat
###########################
for (s in metadata$label) {
  # create internal reference
  count_mat <- GetAssayData(subset(gex_scaled_all, subset = orig.ident == s & My_annotation %in% c("T", "other", "Mono", "NK", "DC")), assay = "RNA", slot = "counts")
  cell_annot <- tibble("cell" = colnames(subset(gex_scaled_all, subset = orig.ident == s & My_annotation %in% c("T", "other", "Mono", "NK", "DC")))) %>%
    dplyr::mutate("group" = "ref")
  ref_internal <- aggregate_counts(count_mat, cell_annot)

  # subset cells to malignent cells
  counts <- GetAssayData(subset(x = gex_scaled_all, subset = orig.ident == s & My_annotation %in% c("MBL B Cell", "CLL B Cell")), assay = "RNA", slot = "counts")
  alleles <- allele_lst[[s]] %>%
    dplyr::filter(cell %in% intersect(colnames(counts), allele_lst[[s]]$cell))
  counts <- counts[, unique(alleles$cell)]

  # run numbat for CNV identification
  out <- run_numbat(
    count_mat = counts,
    # ref_hca,
    ref_internal,
    alleles,
    genome = "hg38",
    t = 1e-5,
    ncores = 6,
    plot = TRUE,
    out_dir = file.path(result_dir,label2pool[s])
  )
}

# ###########################
# #Visualisation and plotting
# ###########################
nb_lst <- map(metadata$label[c(7,13, 18,23, 24,29, 32, 34,35)], ~ Numbat$new(out_dir = file.path(result_dir, label2pool[.]))) %>%
  purrr::set_names(metadata$label[c(7,13, 18,23, 24,29, 32, 34,35)])

genotype_colour <- c("1" = "gray", brewer.pal(n = 8, name = "Set1") %>%
                       set_names(c(1:length(brewer.pal(n = 8, name = "Set1")) + 1)))

phyloheatmaps <- nb_lst %>%
  map(~ .$plot_phylo_heatmap(
    clone_bar = TRUE,
    p_min = 0.9,
    pal_clone = genotype_colour
  ))

imap(phyloheatmaps, ~ ggsave(.x,
              file = file.path(GEX_plotting_dir, "CNV_RNA_phyloheatmaps", paste0("CNV_RNA_phyoloheatmaps_", .y, ".png")),
              device = "png",
              bg = "white",
              width = 6,
              height = 3
))

scatterplots_bulk <- nb_lst %>%
  map(~ .$bulk_clones %>%
        filter(n_cells > 50) %>%
        plot_bulks(
          min_LLR = 10, # filtering CNVs by evidence
          legend = TRUE
        ))

imap(~ ggsave(scatterplots_bulk,
              file = file.path(GEX_plotting_dir, "CNV_RNA_bulk_scatterplots", paste0("CNV_RNA_bulk_scatterplot_", .y, ".png")),
              device = "png",
              bg = "white",
              width = 10,
              height = 7))

all_cnv_umaps <- list()
all_genotype_umaps <- list()

for (s in metadata$label[c(3)]) {
  plist <- list()
  muts <- nb_lst[[s]]$joint_post$seg %>% unique()
  cnv_type <- nb_lst[[s]]$joint_post %>%
    dplyr::select(seg, cnv_state) %>%
    unique() %>%
    deframe()

  gex_scaled_patient <- subset(gex_scaled_all, subset = orig.ident == s)

  clone_opt <- nb_lst[[s]]$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::mutate(clone_opt = as.character(clone_opt)) %>%
    deframe()

  gex_scaled_patient <- AddMetaData(
    object = gex_scaled_patient,
    metadata = clone_opt,
    col.name = "clone_opt"
  )

  genotype_umap <- Seurat::FeaturePlot(gex_scaled_patient,
                                       features = "clone_opt", pt.size = 1
  ) +
    labs(x = "UMAP1", y = "UMAP2", title = "") +
    scale_colour_manual(name = "Genotype", values = genotype_colour) +
    theme_pub() +
    labs(title = label2plottingid[s]) +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "right",
      aspect.ratio = 1
    )
  all_genotype_umaps[[s]] <- genotype_umap

  for (mut in muts) {
    mut_probabilites <- nb_lst[[s]]$joint_post %>%
      dplyr::filter(seg == mut) %>%
      dplyr::select(cell, p_cnv) %>%
      deframe()

    gex_scaled_patient_cnv <- AddMetaData(
      object = gex_scaled_patient,
      metadata = mut_probabilites,
      col.name = "p_cnv"
    )

    cnv_umap <- Seurat::FeaturePlot(gex_scaled_patient_cnv,
                                    features = "p_cnv", , pt.size = 1
    ) +
      labs(x = "UMAP1", y = "UMAP2", title = "") +
      scale_colour_gradientn(
        colours = c("midnightblue", "white", "darkred"),
        name = "Probability",
        na.value = "lightgrey",
        values = c(0, 0.5, 1),
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25)
      ) +
      theme_pub() +
      labs(title = paste0(cnv_type[mut], " ", mut)) +
      theme(
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        aspect.ratio = 1
      )
    plist[[mut]] <- cnv_umap
  }

  plist_all <- ggarrange(plotlist = plist, ncol = ceiling(length(plist) / 3), nrow = ceiling(length(plist) / 3), common.legend = TRUE, legend = "right") %>%
    annotate_figure(top = text_grob(label2plottingid[s], color = "black", face = "plain", size = 12))

  all_cnv_umaps[[s]] <- plist_all

  write_tsv(nb_lst[[s]]$joint_post, file.path(GEX_analysis_dir, "CNV_tables", paste0("numbat_joint_post_", s, ".tsv")))
  write_tsv(nb_lst[[s]]$clone_post, file.path(GEX_analysis_dir, "CNV_tables", paste0("numbat_clone_post_", s, ".tsv")))
}


all_cnv_umaps %>%
  imap(~ ggsave(.,
                file = file.path(GEX_plotting_dir, "CNV_RNA_UMAPs", paste0("CNV_RNA_umap_", .y, ".png")),
                device = "png",
                bg = "white"
  ))

ggarrange(plotlist = all_genotype_umaps, ncol = 5, nrow = 5, legend = "right") %>%
  ggsave(
    file = file.path(project_dir, "plots", paste0("Genotype_RNA_umap.png")),
    device = "png",
    bg = "white", width = 18, height = 18
  )