# Load libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(Seurat)
library(tidyr)
library(SeuratWrappers)
library(progress)
library(reticulate)
library(sceasy)

# Set paths
source("config_paths.R")
results_dir <- file.path(project_dir_remote, "results")
RDS_dir <- file.path(results_dir, "RDS_objects")

# Load data
print("Loading RDS file")
gex_scaled_all <- readRDS(file.path(RDS_dir, "gex_scaled_all.RDS"))
print("Successfully loaded Seurat object")

# If required try first with a subset
#gex_scaled_all <- subset(gex_scaled_all, idents = c("CLL2.1", "CLL2.2","CLL3.1", "CLL3.2", "CLL9.1", "CLL9.2", "MBL10", "MBL15")) # choose subset 
#gex_scaled_all <- FindVariableFeatures(gex_scaled_all, selection.method = "vst", nfeatures = 2000) adapt later to more of only interessted in latent space or all if batch correction on expresssion matrix is desired
#print("Variable Features successfully identified")
#top2000 <- head(VariableFeatures(gex_scaled_all), 200) # adapt later to more of only interessted in latent space or all if batch correction on expresssion matrix is desired
#gex_scaled_all <- subset(gex_scaled_all, features=top2000) # remove once running on all samples and cells
#saveRDS(gex_scaled_all, file=file.path(RDS_dir, "gex_scaled_subsetted.RDS"))
#gex_scaled_all <- readRDS(file=file.path(RDS_dir, "gex_scaled_subsetted.RDS")
Idents(gex_scaled_all) <- "orig.ident"
print("Successfully loaded Seurat object")


# load libraries
use_condaenv("scvi_env", required = TRUE)
print("Successfully loaded scVI conda environment")
scvi <- import("scvi", convert = FALSE)
print("Successfully loaded scVI")
jax <- import("jax", convert = FALSE)
print("Successfully loaded JAX")

# Convert Seurat to AnnData
aData <- convertFormat(gex_scaled_all, from = "seurat", to = "anndata", main_layer = "counts", drop_single_values = FALSE)
aData$X<-aData$X$tocsr() # enables faster processing in row sparse format
print("Successfully converted Seurat to AnnData object")

# Setup model
scvi$model$SCVI$setup_anndata(aData, batch_key="sequencing_date") # could be also added: categorical_covariate_keys="orig.idents" if correction for interpatient heterogeneity desired
model <- scvi$model$SCVI(aData)
print("Successfully created SCVI model object")

# Train model
model$train(accelerator = 'gpu', 
            devices = 'auto', 
            max_epochs=as.integer(20), 
            early_stopping=TRUE) # Should be adapted depending on size of data set
print("Successfully trained model")

# Get latent space
latent <- as.matrix(model$get_latent_representation())
rownames(latent)<- colnames(gex_scaled_all)
saveRDS(latent, file=file.path(RDS_dir, "scvi_latent.RDS"))
print("Successfully pulled latent representation and saved")

# Get batch-projected expression data
expression <- model$get_normalized_expression(aData, transform_batch="2022-11-18") # choose ideally batch which represents as many celltypes as possible for expression to be projected to
expression <- py_to_r(expression) # this is a dataframe with names
saveRDS(expression, file=file.path(RDS_dir, "scvi_expression.RDS"))
print("Successfully pulled batch-corrected expression matrix and saved")
rm(expression)

# Save final Seurat object
#saveRDS(gex_scaled_all, file = file.path(RDS_dir, "gex_scaled_all_scvi.RDS"))
#print("Successfully saved final Seurat object with new UMAP embedding")