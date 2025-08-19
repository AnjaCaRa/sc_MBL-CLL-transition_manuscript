
source("config_paths.R")

# Colour palettes
my_unpaired_colours_named <- c(
  "CLL1" = "darkgoldenrod1",
  "CLL2" = "darkorange1",
  "CLL3" = "indianred1",
  "CLL4" = "plum1",
  "CLL5" = "orchid",
  "CLL6" = "cadetblue1",
  "CLL7" = "turquoise",
  "CLL8" = "darkseagreen1",
  "CLL9" = "burlywood1",
  "CLL10" = "lightsalmon",
  "MBL1" = "darkslategrey",
  "MBL2" = "rosybrown3",
  "MBL3" = "thistle4",
  "MBL4" = "olivedrab",
  "MBL6" = "steelblue4",
  "MBL7" = "slateblue2",
  "MBL7dup" = "lavender",
  "MBL8" = "lightcoral",
  "MBL9" = "yellow3",
  "MBL10" = "lightgoldenrod1",
  "MBL11" = "antiquewhite2",
  "MBL12" = "azure4",
  "MBL13" = "aquamarine4",
  "MBL14" = "darkolivegreen1",
  "MBL15" = "lightsteelblue"
)

plottingid2colours <- c(
  "P15" = "darkgoldenrod1",
  "P17" = "darkorange1",
  "P23" = "indianred1",
  "P20" = "plum1",
  "P18" = "orchid",
  "P14" = "cadetblue1",
  "P16" = "turquoise",
  "P19" = "darkseagreen1",
  "P22" = "burlywood1",
  "P21" = "lightsalmon",
  "P8" = "darkslategrey",
  "P12" = "rosybrown3",
  "P10" = "thistle4",
  "P24" = "olivedrab",
  "P7" = "steelblue4",
  "P9" = "slateblue2",
  "P9dup" = "lavender",
  "P11" = "lightcoral",
  "P13" = "yellow3",
  "P6" = "lightgoldenrod1",
  "P2" = "antiquewhite2",
  "P1" = "azure4",
  "P4" = "aquamarine4",
  "P3" = "darkolivegreen1",
  "P5" = "lightsteelblue"
)

colours <- list(
  "Cell type" = c("Physiologic B Cells" = "darkred", "LC-MBL Cells" = "#fcbda4", "HC-MBL Cells" = "#FC8D62", "CLL Cells" = "#66C2A5"),
  "Patient" = my_unpaired_colours_named,
  "IGHV" = c("mutated" = "steelblue", 
             "unmutated" = "orange", 
             "NA" = "lightgrey"),
  "processing_date" = c(
    "2022-01-10" = "#1B9E77", "2022-10-20" = "#D95F02",
    "2022-11-04" = "#7570B3", "2023-03-24" = "#66A61E",
    "2024-03-05" = "#A6761D"
  ),
  "technology" = c(
    "v2 scATAC-seq" = "#666666", "multiome kit" = "#E7298A",
    "v2 5'scRNA-seq" = "coral3", "v2 HT 5'scRNA-seq" = "lightblue"
  ),
  "sequencing_date" = c("2022-02-28" = "#E6AB02", 
                        "2022-11-18" = "#7FC97F", 
                        "2023-05-02" = "#FFFF99", 
                        "2024-03-18" = "#386CB0",
                        "2022-02-26" = "#FFEF40"),
  
  "sequencing_date_?" = c("2022-01-14" = "midnightblue", 
                          "2022-10-20" = "steelblue4", 
                          "2022-10-26" = "magenta4", 
                          "2022-10-27" = "purple4", 
                          "2024-03-18" = "#386CB0", 
                          "2023-05-02" = "violet", 
                          "2022-11-18" = "plum4",
                          "2022-02-28" = "deeppink4"),
  "sex" = c("male" = "lightblue", "female" = "lightpink")
)



# downsample cell barcodes to equal number per celltype
downsample_cells <- function(df, # one column with label_cellbarcodes and one with celltype name and one with patient label
                             max_cells # cells per celltype
) {
  set.seed(123)
  num_patients <- n_distinct(df$label)
  cells_per_patient <- max_cells %/% num_patients
  remaining_cells <- max_cells %% num_patients
  
  # Sample cells per patient with even distribution
  downsampled <- df %>%
    group_by(label) %>%
    sample_n(size = min(cells_per_patient, n()), replace = FALSE) %>%
    ungroup()
  
  # If there are remaining cells, sample them from the already downsampled data
  if (remaining_cells > 0) {
    remaining <- df %>%
      anti_join(downsampled, by = "label_cellbarcode") %>%
      sample_n(size = min(remaining_cells, n()), replace = FALSE)
    
    downsampled <- bind_rows(downsampled, remaining)
  }
  
  return(downsampled)
}

# create complex heatmap with columns separated by celltype according to patient clustering
build_celltype_heatmap <- function(bcs_per_celltype, # named list with celltype name, rows are cells, columns are label_cellbarcode and B_like_annotation
                                   data_matrix, # rownames are genes, colnames are label_cellbarcodes, rownames are already subsetted for GOIs
                                   metadata_all,
                                   modality, # ATAC or GEX?
                                   include_all_metadata, # TRUE, FALSE; if false only celltype, patient and IGHV mutational status plotted
                                   show_heatmap_legend, # TRUE, FALSE,
                                   show_annotation_legend, # TRUE FALSE
                                   quantile_capping, # values between 0 and 1
                                   absolute_capping, # absolute values to cap on should be a numeric vector with min and max c(min, max)
                                   high_score_colour, # "#colourvalue"
                                   mid_score_colour, # "#colourvalue"
                                   low_score_colour, # "#colourvalue"
                                   heatmap_legend_name # z-score, expression or similar
) {
  metadata_all <- read_tsv(file = file.path(project_dir, "data/20250403_Tidy_MBL_CLL_metadata.tsv"))
  
  
  
  # Calculate mean score per patient and cluster to determine patient and cell order in plot
  clustering_input <- bcs_per_celltype %>%
    map(~ data_matrix[, .$label_cellbarcode] %>%
          as_tibble(rownames = NA) %>%
          rownames_to_column("TF") %>%
          pivot_longer(cols = -TF, names_to = "label_cellbarcode", values_to = "score") %>%
          tidyr::separate(label_cellbarcode, c("label", "cellbarcode"), sep = "#") %>%
          group_by(TF, label) %>%
          dplyr::summarise("mean_score" = mean(`score`)) %>%
          ungroup() %>%
          pivot_wider(id_cols = label, names_from = TF, values_from = mean_score) %>%
          column_to_rownames(var = "label")) %>%
    set_names(names(bcs_per_celltype))
  
  set.seed(20)
  distance_mat <- clustering_input %>%
    map(~ dist(., method = "euclidean"))
  clusters <- distance_mat %>%
    map(~ hclust(., method = "average"))
  
  patients_ordered <- clustering_input %>%
    map2(clusters, ~ rownames(.x)[.y %>%
                                    as.dendrogram() %>%
                                    order.dendrogram()])
  
  bcs_per_celltype_ordered <- bcs_per_celltype %>%
    map2(patients_ordered, ~ dplyr::mutate(.x, label = factor(label, levels = .y)) %>%
           arrange(label))
  
  # Top annotation rows
  barcodes2annotation <- bcs_per_celltype %>%
    map(~ dplyr::select(., label_cellbarcode, B_like_annotation) %>%
          dplyr::mutate(B_like_annotation = factor(B_like_annotation, levels = c("Physiologic B Cells", "LC-MBL Cells", "HC-MBL Cells", "CLL Cells"))) %>%
          deframe())
  
  barcodes2label <- bcs_per_celltype %>%
    map2(patients_ordered, ~ left_join(.x, tibble(
      "label_cellbarcode" = colnames(data_matrix[, .x$label_cellbarcode]),
      "label" = colnames(data_matrix[, .x$label_cellbarcode]) %>%
        as_tibble() %>%
        dplyr::rename("label_cellbarcode" = value) %>%
        tidyr::separate(label_cellbarcode, c("label", "cellbarcode"), sep = "#", remove = FALSE) %>%
        pull(label)
    )) %>%
      dplyr::select(label_cellbarcode, label) %>%
      dplyr::mutate(label = factor(label, levels = .y)) %>%
      deframe())
  
  barcodes2patient <- bcs_per_celltype %>%
    map2(patients_ordered, ~ left_join(.x, tibble(
      "label_cellbarcode" = colnames(data_matrix[, .x$label_cellbarcode]),
      "label" = colnames(data_matrix[, .x$label_cellbarcode]) %>%
        as_tibble() %>%
        dplyr::rename("label_cellbarcode" = value) %>%
        tidyr::separate(col = label_cellbarcode, into = c("label", "cellbarcode"), sep = "#", remove = FALSE) %>%
        pull(label)
    )) %>%
      left_join(metadata_all %>%
                  dplyr::filter(data_type == modality) %>%
                  dplyr::select("label", "patient_id")) %>%
      dplyr::select(label_cellbarcode, patient_id) %>%
      deframe())
  
  barcodes2IGHV <- bcs_per_celltype %>%
    map2(patients_ordered, ~ left_join(.x, metadata_all %>%
                                         dplyr::filter(data_type == modality) %>%
                                         dplyr::select("label", "IGHV")) %>%
           dplyr::mutate(
             label = factor(label, levels = .y),
             IGHV = if_else(is.na(IGHV), "NA", IGHV)
           ) %>%
           dplyr::select(label_cellbarcode, IGHV) %>%
           deframe())
  
  barcode2processing_date <- bcs_per_celltype %>%
    map2(patients_ordered, ~ left_join(.x, metadata_all %>%
                                         dplyr::filter(data_type == modality) %>%
                                         dplyr::select("label", "processing_date")) %>%
           dplyr::mutate(
             label = factor(label, levels = .y),
             processing_date = as.character(processing_date)
           ) %>%
           dplyr::select(label_cellbarcode, processing_date) %>%
           deframe())
  
  barcode2technology <- bcs_per_celltype %>%
    map2(patients_ordered, ~ left_join(.x, metadata_all %>%
                                         dplyr::filter(data_type == modality) %>%
                                         dplyr::select("label", "technology")) %>%
           dplyr::mutate(label = factor(label, levels = .y)) %>%
           dplyr::select(label_cellbarcode, technology) %>%
           deframe())
  
  barcode2sequencing_date <- bcs_per_celltype %>%
    map2(patients_ordered, ~ left_join(.x, metadata_all %>%
                                         dplyr::filter(data_type == modality) %>%
                                         dplyr::select("label", "sequencing_date")) %>%
           dplyr::mutate(
             label = factor(label, levels = .y),
             sequencing_date = as.character(sequencing_date)
           ) %>%
           dplyr::select(label_cellbarcode, sequencing_date) %>%
           deframe())
  
  barcode2sex <- bcs_per_celltype %>%
    map2(patients_ordered, ~ left_join(.x, metadata_all %>%
                                         dplyr::filter(data_type == modality) %>%
                                         dplyr::select("label", "sex")) %>%
           dplyr::mutate(label = factor(label, levels = .y)) %>%
           dplyr::select(label_cellbarcode, sex) %>%
           deframe())
  
  "barcodes2sequencing_date_?" <- bcs_per_celltype %>%
    map2(patients_ordered, ~ left_join(.x, metadata_all %>%
                                         dplyr::filter(data_type == modality) %>%
                                         dplyr::select("label", "sequencing_date_?")) %>%
           dplyr::mutate(
             label = factor(label, levels = .y),
             sequencing_date = as.character(`sequencing_date_?`)
           ) %>%
           dplyr::select(label_cellbarcode, `sequencing_date_?`) %>%
           deframe())
  
  if (include_all_metadata == TRUE) {
    anno_df <- list(
      barcodes2annotation,
      barcodes2patient,
      barcodes2IGHV,
      barcode2processing_date,
      barcode2technology,
      barcode2sequencing_date,
      `barcodes2sequencing_date_?`,
      barcode2sex,
      bcs_per_celltype
    ) %>%
      pmap(~ tibble(
        "Cell type" = ..1[colnames(data_matrix[, ..9$label_cellbarcode])],
        "Patient" = ..2[colnames(data_matrix[, ..9$label_cellbarcode])],
        "IGHV" = ..3[colnames(data_matrix[, ..9$label_cellbarcode])],
        "processing_date" = ..4[colnames(data_matrix[, ..9$label_cellbarcode])],
        "technology" = ..5[colnames(data_matrix[, ..9$label_cellbarcode])],
        "sequencing_date" = ..6[colnames(data_matrix[, ..9$label_cellbarcode])],
        "sequencing_date_?" = ..7[colnames(data_matrix[, ..9$label_cellbarcode])],
        "sex" = ..8[colnames(data_matrix[, ..9$label_cellbarcode])]
      ))
  } else {
    anno_df <- list(
      barcodes2annotation,
      barcodes2patient,
      barcodes2IGHV,
      barcode2processing_date,
      barcode2technology,
      barcode2sequencing_date,
      `barcodes2sequencing_date_?`,
      barcode2sex,
      bcs_per_celltype
    ) %>%
      pmap(~ tibble(
        "Cell type" = ..1[colnames(data_matrix[, ..9$label_cellbarcode])],
        "Patient" = ..2[colnames(data_matrix[, ..9$label_cellbarcode])],
        "IGHV" = ..3[colnames(data_matrix[, ..9$label_cellbarcode])]
      ))
  }
  
  column_ha <- anno_df %>%
    map(~ HeatmapAnnotation(
      df = .,
      col = colours,
      show_legend = show_annotation_legend,
      annotation_name_side = "left",
      annotation_name_gp = gpar(fontsize = 8),
      border = TRUE
    ),
    annotation_height = unit(4, "mm"),
    gap = unit(0.5, "mm"),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10,
                              fontfamily ="Arial"))
  
  colour_breaks <- c()
  if (is.numeric(quantile_capping)) {
    colour_breaks <- colorRamp2(
      breaks = c(-quantile(as.matrix(data_matrix), quantile_capping),
                 0,
                 quantile(as.matrix(data_matrix), quantile_capping)),
      c(low_score_colour, mid_score_colour, high_score_colour)
    )
  } else{
    colour_breaks <-colorRamp2(breaks = c(absolute_capping[1], 0, absolute_capping[2]),
                               c(low_score_colour, mid_score_colour, high_score_colour))
  }
  
  hm <- list(
    bcs_per_celltype,
    column_ha,
    bcs_per_celltype_ordered %>% 
      map(~ pull(., label_cellbarcode))
  ) %>%
    pmap(~ Heatmap(as.matrix(data_matrix[, ..1$label_cellbarcode]),
                   col = colour_breaks,
                   top_annotation = ..2,
                   column_order = ..3,
                   use_raster = TRUE,
                   cluster_rows = TRUE,
                   cluster_columns = FALSE,
                   raster_by_magick = TRUE,
                   show_column_names = FALSE,
                   show_heatmap_legend = show_heatmap_legend,
                   na_col = "black",
                   row_title_side = "left",
                   row_names_side = "left",
                   row_names_gp = gpar(fontsize = 8),
                   row_title_gp = gpar(fontsize = 8),
                   border_gp = gpar(col = "black", lty = 1),
                   heatmap_legend_param = list(
                     title = heatmap_legend_name,
                     direction = "horizontal",
                     legend_direction = "vertical",
                     legend_gp = gpar(
                       fontsize = 8,
                       fontface = "plain"
                     )
                   )
    ))
  
  return(hm)
}


# create complex heatmap without columns being separated and clustering by heatmap values
build_clustered_heatmap <- function(bcs_cells_df, # dataframe, rows are cells, columns are label_cellbarcode and B_like_annotation
                                    data_matrix, # rownames are genes, colnames are label_cellbarcodes, rownames are already subsetted for GOIs
                                    modality, # ATAC or GEX?
                                    include_all_metadata, # TRUE, FALSE; if false only celltype, patient and IGHV mutational status plotted
                                    show_heatmap_legend, # TRUE, FALSE,
                                    show_annotation_legend, # TRUE FALSE
                                    quantile_capping, # values between 0 and 1
                                    absolute_capping, # absolute values to cap on should be a numeric vector with min and max c(min, max)
                                    high_score_colour, # "#colourvalue"
                                    mid_score_colour, # "#colourvalue"
                                    low_score_colour, # "#colourvalue"
                                    heatmap_legend_name # z-score, expression or similar
) {
  
  metadata_all <- read_tsv(file = file.path(project_dir, "data/20250403_Tidy_MBL_CLL_metadata.tsv"))
  
  
  barcodes2annotation <- bcs_cells_df %>%
    dplyr::select(label_cellbarcode, B_like_annotation) %>%
    dplyr::mutate(B_like_annotation = factor(B_like_annotation, levels = c("Physiologic B Cells", "LC-MBL Cells", "HC-MBL Cells", "CLL Cells"))) %>%
    deframe()
  
  barcodes2patient <- bcs_cells_df %>%
    left_join(tibble(
      "label_cellbarcode" = colnames(data_matrix[, .$label_cellbarcode]),
      "label" = colnames(data_matrix[, .$label_cellbarcode]) %>%
        as_tibble() %>%
        dplyr::rename("label_cellbarcode" = value) %>%
        tidyr::separate(label_cellbarcode, c("label", "cellbarcode"), sep = "#", remove = FALSE) %>%
        pull(label)
    )) %>%
    left_join(metadata_all %>%
                dplyr::filter(data_type == modality) %>%
                dplyr::select("label", "patient_id")) %>%
    dplyr::select(label_cellbarcode, patient_id) %>%
    deframe()
  
  barcodes2IGHV <- bcs_cells_df %>%
    left_join(metadata_all %>%
                dplyr::filter(data_type == modality) %>%
                dplyr::select("label", "IGHV")) %>%
    dplyr::mutate(IGHV = if_else(is.na(IGHV), "NA", IGHV)) %>%
    dplyr::select(label_cellbarcode, IGHV) %>%
    deframe()
  
  barcode2processing_date <- bcs_cells_df %>%
    left_join(metadata_all %>%
                dplyr::filter(data_type == modality) %>%
                dplyr::select("label", "processing_date")) %>%
    dplyr::mutate(
      processing_date = as.character(processing_date)
    ) %>%
    dplyr::select(label_cellbarcode, processing_date) %>%
    deframe()
  
  barcode2technology <- bcs_cells_df %>%
    left_join(metadata_all %>%
                dplyr::filter(data_type == modality) %>%
                dplyr::select("label", "technology")) %>%
    dplyr::select(label_cellbarcode, technology) %>%
    deframe()
  
  barcode2sequencing_date <- bcs_cells_df %>%
    left_join(metadata_all %>%
                dplyr::filter(data_type == modality) %>%
                dplyr::select("label", "sequencing_date")) %>%
    dplyr::mutate(
      sequencing_date = as.character(sequencing_date)
    ) %>%
    dplyr::select(label_cellbarcode, sequencing_date) %>%
    deframe()
  
  barcode2sex <- bcs_cells_df %>%
    left_join(metadata_all %>%
                dplyr::filter(data_type == modality) %>%
                dplyr::select("label", "sex")) %>%
    dplyr::select(label_cellbarcode, sex) %>%
    deframe()
  
  if (include_all_metadata == TRUE) {
    
    anno_df <- list(
      list(barcodes2annotation),
      list(barcodes2patient),
      list(barcodes2IGHV),
      list(barcode2processing_date),
      list(barcode2technology),
      list(barcode2sequencing_date),
      list(barcode2sequencing_date),
      list(barcode2sex),
      list(bcs_cells_df)
    ) %>%
      pmap(~ tibble(
        "Cell type" = ..1[colnames(data_matrix[, ..9$label_cellbarcode])],
        "Patient" = ..2[colnames(data_matrix[, ..9$label_cellbarcode])],
        "IGHV" = ..3[colnames(data_matrix[, ..9$label_cellbarcode])],
        "processing_date" = ..4[colnames(data_matrix[, ..9$label_cellbarcode])],
        "technology" = ..5[colnames(data_matrix[, ..9$label_cellbarcode])],
        "sequencing_date" = ..6[colnames(data_matrix[, ..9$label_cellbarcode])],
        "sequencing_date_?" = ..7[colnames(data_matrix[, ..9$label_cellbarcode])],
        "sex" = ..8[colnames(data_matrix[, ..9$label_cellbarcode])]
      ))} else{
        
        anno_df <- list(
          list(barcodes2annotation),
          list(barcodes2patient),
          list(barcodes2IGHV),
          list(barcode2processing_date),
          list(barcode2technology),
          list(barcode2sequencing_date),
          list(barcode2sequencing_date),
          list(barcode2sex),
          list(bcs_cells_df)
        ) %>%
          pmap(~ tibble(
            "Cell type" = ..1[colnames(data_matrix[, ..9$label_cellbarcode])],
            "Patient" = ..2[colnames(data_matrix[, ..9$label_cellbarcode])],
            "IGHV" = ..3[colnames(data_matrix[, ..9$label_cellbarcode])]
          ))}
  
  
  column_ha <- anno_df %>%
    map(~ HeatmapAnnotation(
      df = .x,
      col = colours,
      show_legend = show_annotation_legend,
      annotation_name_side = "left",
      annotation_name_gp = gpar(fontsize = 8),
      border=TRUE
    ))
  
  colour_breaks <- c()
  if (is.numeric(quantile_capping)) {
    colour_breaks <- colorRamp2(
      breaks = c(-quantile(as.matrix(data_matrix), quantile_capping),
                 0,
                 quantile(as.matrix(data_matrix), quantile_capping)),
      c(low_score_colour, mid_score_colour, high_score_colour)
    )
  } else{
    colour_breaks <-colorRamp2(breaks = c(absolute_capping[1], 0, absolute_capping[2]),
                               c(low_score_colour, mid_score_colour, high_score_colour))
  }
  
  hm <- list(
    list(bcs_cells_df),
    column_ha
  ) %>%
    pmap(~ Heatmap(as.matrix(data_matrix[, ..1$label_cellbarcode]),
                   col = colour_breaks,
                   top_annotation = ..2,
                   use_raster = TRUE,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   raster_by_magick = TRUE,
                   show_column_names = FALSE,
                   show_heatmap_legend = show_heatmap_legend,
                   na_col = "black",
                   row_title_side = "left",
                   row_names_side = "left",
                   row_names_gp = gpar(fontsize = 8),
                   row_title_gp = gpar(fontsize = 8),
                   border_gp = gpar(col = "black", lty = 1),
                   heatmap_legend_param = list(
                     title = heatmap_legend_name,
                     direction = "horizontal",
                     legend_direction = "vertical",
                     legend_gp = gpar(
                       fontsize = 8,
                       fontface = "plain"
                     )
                   )
    ))
  
  return(hm)
}