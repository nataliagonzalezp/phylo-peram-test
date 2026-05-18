################################################################################
# Script: Phylogenetic Peramorphosis Test with Bootstrap Recovery
#
# Description:
# This script implements a phylogenetic extension of the peramorphosis test (Piras et al. 2011)
# to infer heterochronic events (PERAMORPHOSIS, PAEDOMORPHOSIS) across
# ancestor-descendant pairs in a phylogenetic tree using geometric morphometric
# data. 
# The workflow includes:
#   (1) Direct extraction of coordinates and age/developmental progress indicator from .TPS
#   (2) Estimation of a juvenile archetype shape
#   (3) Estimation of species-specific offset shape via Procrustes regression
#   (4) Reconstruction of ancestral shapes using Maximum Likelihood on phylogenetic PCA scores
#   (5) Inference of heterochronic events based on distance comparisons
#   (6) Bootstrap resampling to assess the robustness of inferred patterns
#
# Input:
#   - .TPS file containing landmark coordinates and an age/developmental progress indicator
#   - Phylogenetic tree file in Newick format (.tre)
#
# Output:
#   - PDF plots of heterochronic events onto the phylogeny
#   - CSV table with inferred events and bootstrap recovery percentages
#   - CSV table with raw bootstrap results
#
# Key Parameters:
#   - resampling_proportion: proportion of specimens retained per species during
#   resampling iteration.
#       default value: 0.66
#
#   - diff_ratio_threshold: threshold for inferring heterochronic events
#       default value: 0.10
#
#   - n_iter: number of valid bootstrap replicates
#       default value: 100
#
# Recovery percentage:
#    number of bootstrap replicates matching the observed trend /
#    number of valid bootstrap replicates * 100
#
# Dependencies:
#   abind, ape, borealis, devtools, dplyr, geiger, geomorph, ggplot2,
#   ggpubr, ggthemes, gridExtra, Morpho, phytools, readr, tidyr, tidyverse, viridis
#
# Developed with:
#   R version 4.5.3 (2026-03-11) — "Reassured Reassurer"
#
# Maintainers:
#   Natalia González-Piñeres - Becaria doctoral (UEL-FML-CONICET)
#
# Affiliation:
#   Unidad Ejecutora Lillo - Fundación Miguel Lillo (UEL-FML)
#   Consejo Nacional de Investigaciones Científicas y Técnicas (CONICET)
#
# Notes:
#   - The script assumes consistent species naming between TPS and tree file.
#   - Bootstrap replicates include safeguards for low sample sizes and invalid fits.
#   - Species represented by exactly two specimens are kept fixed during resampling.
#   A warning message is delivered in that case.
#   - Bootstrap samples must retain at least two specimens per species and size
#     variation within every species. (Enough for pPCA to run, in absence of variation the function crashes.)
#   - Incomplete bootstrap replicates are not counted as valid.
#   - Output filenames are dynamically generated based on parameter values.
#
# Citation:
#   Please cite this script along with:
#   González-Piñeres & Catalano (2026, in prep.)
#   Piras et al. (2011), and relevant geomorph references.
#
# Version:
#   v1.0 
#
################################################################################


# ============================================================
# CLEAN WORKING ENVIRONMENT
# ============================================================

rm(list = ls())


# ============================================================
# LOAD LIBRARIES
# ============================================================

library(abind)
library(ape)
library(borealis)
library(devtools)
library(dplyr)
library(geiger)
library(geomorph)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(gridExtra)
library(Morpho)
library(phytools)
library(readr)
library(tidyr)
library(tidyverse)
library(viridis)


# ============================================================
# USER-DEFINED PARAMETERS
# ============================================================

# Complete the names of the corresponding input files.
TPS_file <- "YourCoordinatesFile.tps"
tree_file <- "YourTreeFile.tre"

tps_specID <- "ID"
tps_devtime_keyword <- "AGE"

# coords_aligned: set to TRUE if the TPS file contains Procrustes-aligned
# (superimposed) coordinates. Set to FALSE if the TPS contains raw (unaligned)
# landmark coordinates — gpagen() will be run automatically.
coords_aligned <- TRUE

# has_age_label: set to TRUE if the TPS file contains an AGE= field for each
# specimen. Set to FALSE if no AGE field is present - log(centroid size) from
# gpagen() will be used as the developmental indicator instead.
# NOTE: has_age_label = FALSE requires coords_aligned = FALSE.
# If coords_aligned = TRUE and has_age_label = FALSE, the script will stop
# and ask you to provide a TPS file with AGE= labels.
has_age_label <- TRUE

resampling_proportion <- 0.66
diff_ratio_threshold <- 0.10
n_iter <- 100

set.seed(5422)

sampling_percent <- resampling_proportion * 100
diff_threshold_percent <- diff_ratio_threshold * 100

output_prefix <- paste0(
 "Bootstrap",
 sampling_percent,
 "_",
 diff_threshold_percent,
 "_DThreshold"
)

observed_pdf_name <- paste0(
 "HeterochronicChanges_",
 diff_threshold_percent,
 "_DThreshold.pdf"
)

bootstrap_pdf_name <- paste0(
 "HeterochronicChanges_",
 output_prefix,
 ".pdf"
)

results_csv_name <- paste0(
 sampling_percent,
 "Keep",
 diff_threshold_percent,
 "Threshold_TotalSample_Heterochrony_Results_WithRecovery.csv"
)

raw_bootstrap_csv_name <- paste0(
 sampling_percent,
 "Keep",
 diff_threshold_percent,
 "Threshold_RawBootstrapTrends.csv"
)

observed_plot_title <- paste0(
 "Heterochronic Changes with ",
 diff_threshold_percent,
 "% distance threshold"
)

bootstrap_plot_title <- paste0(
 "Heterochronic changes with Bootstrap Recovery ",
 sampling_percent,
 "% keep, ",
 diff_threshold_percent,
 "% distance threshold"
)


# ============================================================
# FUNCTION: EXTRACT DATA DIRECTLY FROM TPS FILE
# ============================================================
#
# This function does NOT write .RData files.
# This function does NOT write input metadata CSV files.
#
# It only reads the .TPS file and creates an in-memory geomorph.data.frame
# containing:
#   - Data    = landmark coordinates
#   - Devtime    = age/developmental indicator extracted from AGE=
#   - Species = species labels extracted from TPS IDs
#   - ID      = original TPS IDs
#
# If your TPS IDs are already species names, keep:
#   Species <- ID
#
# If your TPS IDs contain specimen-level suffixes, modify the Species line.
# Example:
#   Species <- sub("_IND_[0-9]+$", "", ID)

extract_tps_input_data <- function(
  tps_file,
  specID = "ID",
  devtime_keyword = "AGE",
  read_devtime = TRUE
) {

 coords <- geomorph::readland.tps(
  file = tps_file,
  specID = specID
 )

 file_lines <- readLines(tps_file)

 if (read_devtime) {

  devtime_pattern <- paste0("^", devtime_keyword, "=")
  devtime_lines <- grep(devtime_pattern, file_lines, value = TRUE)

  if (length(devtime_lines) != dim(coords)[3]) {
   stop(
    "The number of ",
    devtime_keyword,
    " values does not match the number of specimens in the TPS file.\n",
    "Number of specimens in coordinates: ",
    dim(coords)[3],
    "\nNumber of ",
    devtime_keyword,
    " lines: ",
    length(devtime_lines)
   )
  }

  Devtime <- as.numeric(sub(devtime_pattern, "", devtime_lines))

  if (any(is.na(Devtime))) {
   stop(
    "Some ",
    devtime_keyword,
    " values could not be converted to numeric."
   )
  }

 } else {
  Devtime <- NULL
 }

 specimen_names <- dimnames(coords)[[3]]

 if (is.null(specimen_names)) {
  stop("No specimen names were found in the TPS file.")
 }

 ID <- specimen_names

 # Original logic:
 # Species are extracted directly from the TPS ID.
 # Modify this line only if your TPS IDs include individual suffixes.
 Species <- ID

 if (!is.null(Devtime)) {
  gdf <- geomorph::geomorph.data.frame(
   Data = coords,
   Devtime = Devtime,
   Species = Species,
   ID = ID
  )
 } else {
  gdf <- geomorph::geomorph.data.frame(
   Data = coords,
   Species = Species,
   ID = ID
  )
 }

 message("TPS file successfully read: ", tps_file)
 message("Number of specimens: ", dim(coords)[3])
 message("Number of landmarks: ", dim(coords)[1])
 message("Number of dimensions: ", dim(coords)[2])
 message("Number of species/labels: ", length(unique(Species)))

 return(gdf)
}


# ============================================================
# READ TPS FILE DIRECTLY
# ============================================================

gdf <- extract_tps_input_data(
 tps_file = TPS_file,
 specID = tps_specID,
 devtime_keyword = tps_devtime_keyword,
 read_devtime = has_age_label
)


# ============================================================
# ALIGN COORDINATES AND ASSIGN DEVELOPMENTAL INDICATOR
# ============================================================
#
# Behavior is controlled by coords_aligned and has_age_label:
#
#   coords_aligned = TRUE,  has_age_label = TRUE
#     Coordinates are pre-aligned and AGE labels are present.
#     No alignment is performed. AGE values are used as Devtime.
#
#   coords_aligned = FALSE, has_age_label = FALSE
#     Raw (unaligned) coordinates, no AGE field.
#     gpagen() aligns the coordinates and centroid sizes are extracted.
#     log(centroid size) is assigned as the developmental indicator (Devtime).
#
#   coords_aligned = FALSE, has_age_label = TRUE
#     Raw coordinates but AGE labels are present.
#     gpagen() aligns the coordinates. The provided AGE values are used as Devtime.
#
#   coords_aligned = TRUE,  has_age_label = FALSE
#     Aligned coordinates but no developmental indicator. Cannot proceed.
#     The script stops and asks the user to supply a TPS file with AGE= labels.

if (!coords_aligned && !has_age_label) {

 message(
  "coords_aligned = FALSE, has_age_label = FALSE.\n",
  "Running gpagen() on raw landmark coordinates..."
 )

 gpagen_result <- geomorph::gpagen(
  gdf$Data,
  print.progress = FALSE
 )

 gdf$Data <- gpagen_result$coords

 message("Coordinates successfully aligned with gpagen().")
 message("Computing log(centroid size) as developmental indicator...")

 gdf$Devtime <- log(gpagen_result$Csize)

 message("log(centroid size) assigned as Devtime for all specimens.")

} else if (coords_aligned && !has_age_label) {

 stop(
  "coords_aligned = TRUE but has_age_label = FALSE.\n",
  "The script cannot proceed without a developmental indicator (AGE= field).\n",
  "Please provide a TPS file that includes AGE= values for each specimen,\n",
  "or set coords_aligned = FALSE to let the script align the coordinates\n",
  "and derive log(centroid size) as the developmental indicator."
 )

} else if (!coords_aligned && has_age_label) {

 message(
  "coords_aligned = FALSE: running gpagen() on raw landmark coordinates.\n",
  "The AGE labels provided in the TPS file will be used as the developmental indicator."
 )

 gpagen_result <- geomorph::gpagen(
  gdf$Data,
  print.progress = FALSE
 )

 gdf$Data <- gpagen_result$coords

 message("Coordinates successfully aligned with gpagen().")

} else {

 message(
  "coords_aligned = TRUE, has_age_label = TRUE.",
  " Proceeding with pre-aligned coordinates and provided AGE labels."
 )

}


# ============================================================
# LOAD TREE
# ============================================================

tree <- ape::read.tree(tree_file)

phylo_tree <- tree
phylo_tree$node.label <- NULL

if (is.null(phylo_tree$edge.length)) {
 warning(
  "The tree does not contain branch lengths.\n",
  "All branch lengths have been set to 1.\n",
  "The analysis will proceed with an ultrametric tree (equal unit branch lengths)."
 )
 phylo_tree$edge.length <- rep(1, nrow(phylo_tree$edge))
}

tree_edges <- phylo_tree$edge
tip_labels <- phylo_tree$tip.label


# ============================================================
# LOAD DATA FROM TPS-DERIVED GDF
# ============================================================

Coords <- gdf$Data
CoordinatesData <- gdf$Data

if (is.null(CoordinatesData)) {
 stop("The geomorph object must contain coordinate data named `Data`.")
}

Nlandmarks <- dim(CoordinatesData)[1]
Ndimensions <- dim(CoordinatesData)[2]

SpeciesData <- as.factor(gdf$Species)

Species <- levels(as.factor(gdf$Species))
SpecimenNames <- character(length(gdf$Species))

for (species in Species) {
 
 indices <- which(gdf$Species == species)
 
 new_names <- paste0(
  species,
  "_IND_",
  sprintf("%02d", seq_along(indices))
 )
 
 SpecimenNames[indices] <- new_names
}

gdf$ID <- SpecimenNames
SpecimensData <- as.factor(gdf$ID)


# ============================================================
# STANDARDIZE CENTROID SIZE / DEVELOPMENTAL INDICATOR
# ============================================================

DevtimeTwo <- gdf$Devtime

gdf$OrigDevtime <- as.numeric(DevtimeTwo)

if (any(is.na(gdf$OrigDevtime))) {
 stop("The Devtime variable contains NA values or non-numeric values.")
}

scaled_devtime <- numeric(length(gdf$OrigDevtime))

for (sp in unique(gdf$Species)) {
 
 idx <- which(gdf$Species == sp)
 sp_devtimes <- gdf$OrigDevtime[idx]
 
 if (length(idx) < 2) {
  stop(
   "Species `",
   sp,
   "` has fewer than 2 specimens. ",
   "Species-specific regressions require at least 2 specimens per species."
  )
 }
 
 if (length(unique(sp_devtimes)) < 2) {
  stop(
   "Species `",
   sp,
   "` has no variation in Devtime/AGE. ",
   "The regression cannot be fitted for this species."
  )
 }
 
 scaled_devtime[idx] <- (sp_devtimes - min(sp_devtimes)) /
  (max(sp_devtimes) - min(sp_devtimes))
}

gdf$Devtime <- scaled_devtime
gdf$Devtime[gdf$Devtime == 0] <- 1e-7

DevtimeTwo <- gdf$Devtime


# ============================================================
# SORT DATA BY SPECIES
# ============================================================

sort_df <- data.frame(
 species = SpeciesData,
 specimen = SpecimensData,
 index = seq_along(SpeciesData)
)

sort_df_sorted <- sort_df[order(sort_df$species), ]
new_order <- sort_df_sorted$index

CoordinatesData_sorted <- CoordinatesData[, , new_order, drop = FALSE]
DevtimeTwo_sorted <- DevtimeTwo[new_order]
SpeciesData_sorted <- SpeciesData[new_order]
SpecimensData_sorted <- SpecimensData[new_order]

gdf_sorted <- geomorph::geomorph.data.frame(
 Data = CoordinatesData_sorted,
 Devtime = DevtimeTwo_sorted,
 Species = SpeciesData_sorted,
 ID = SpecimensData_sorted
)

gdf2 <- gdf_sorted

CoordinatesData <- gdf2$Data
DevtimeTwo <- gdf2$Devtime
SpeciesData <- as.factor(gdf2$Species)
SpecimensData <- as.factor(gdf2$ID)


# ============================================================
# CHECK TREE AND DATA MATCH
# ============================================================

species_in_data <- levels(SpeciesData)
species_in_tree <- phylo_tree$tip.label

missing_from_data <- setdiff(species_in_tree, species_in_data)
missing_from_tree <- setdiff(species_in_data, species_in_tree)

if (length(missing_from_data) > 0) {
 stop(
  "The following tree tips are missing from the morphometric data:\n",
  paste(missing_from_data, collapse = ", ")
 )
}

if (length(missing_from_tree) > 0) {
 warning(
  "The following species are present in the data but absent from the tree:\n",
  paste(missing_from_tree, collapse = ", "),
  "\nThey will not be used in the phylogenetic reconstruction."
 )
}


# ============================================================
# FUNCTIONS
# ============================================================

classify_heterochrony <- function(
  dist_A,
  dist_B,
  diff_ratio_threshold
) {
 
 diff_ratio <- abs(dist_A - dist_B) / mean(c(dist_A, dist_B))
 
 trend <- if (diff_ratio < diff_ratio_threshold) {
  "NC"
 } else if (dist_A > dist_B) {
  "PAEDOMORPHOSIS"
 } else {
  "PERAMORPHOSIS"
 }
 
 return(
  list(
   trend = trend,
   diff_ratio = diff_ratio
  )
 )
}


estimate_max_shapes <- function(
  coords,
  devtimes,
  species,
  max_devtimes
) {
 
 shape_list <- lapply(levels(species), function(sp) {
  
  idx <- which(species == sp)
  
  if (length(idx) < 2) {
   stop(
    "Species `",
    sp,
    "` has fewer than 2 specimens in this sample."
   )
  }
  
  coords_sp <- coords[, , idx, drop = FALSE]
  devtime_sp <- devtimes[idx]
  
  if (length(unique(devtime_sp)) < 2) {
   stop(
    "Species `",
    sp,
    "` has no size variation in this sample."
   )
  }
  
  gdf_sp <- geomorph::geomorph.data.frame(
   coords = coords_sp,
   DevtimeTwo = devtime_sp
  )
  
  fit_sp <- geomorph::procD.lm(
   coords ~ DevtimeTwo,
   data = gdf_sp
  )
  
  max_devtime <- max_devtimes$MaxDevtime[
   as.character(max_devtimes$Species) == sp
  ]
  
  if (length(max_devtime) == 0 || is.na(max_devtime)) {
   stop(
    "No maximum size found for species `",
    sp,
    "`."
   )
  }
  
  predicted_vector <- predict(
   fit_sp,
   newdata = data.frame(DevtimeTwo = max_devtime)
  )$mean
  
  matrix(
   predicted_vector,
   ncol = Ndimensions,
   byrow = TRUE
  )
 })
 
 names(shape_list) <- levels(species)
 
 arr <- array(
  NA,
  dim = c(
   Nlandmarks,
   Ndimensions,
   length(shape_list)
  )
 )
 
 for (j in seq_along(shape_list)) {
  arr[, , j] <- shape_list[[j]]
 }
 
 dimnames(arr)[[3]] <- names(shape_list)
 
 return(arr)
}


reconstruct_ancestral_shapes <- function(
  predicted_shapes,
  phylo_tree
) {
 
 pca_obj <- geomorph::gm.prcomp(
  predicted_shapes,
  phy = phylo_tree
 )
 
 if (is.null(pca_obj$ancestors) || nrow(pca_obj$ancestors) == 0) {
  stop("No ancestral states were recovered by gm.prcomp().")
 }
 
 anc_mat <- as.matrix(pca_obj$ancestors)
 
 ancestral_coords <- array(
  t(anc_mat),
  dim = c(
   Ndimensions,
   Nlandmarks,
   nrow(anc_mat)
  )
 )
 
 ancestral_coords <- aperm(
  ancestral_coords,
  c(2, 1, 3)
 )
 
 dimnames(ancestral_coords)[[3]] <- rownames(anc_mat)
 
 return(ancestral_coords)
}


infer_edge_trends <- function(
  shape_0,
  predicted_shapesMAX,
  AncestralCoordinatesMAX,
  phylo_tree,
  diff_ratio_threshold
) {
 
 tree_edges_local <- phylo_tree$edge
 tip_labels_local <- phylo_tree$tip.label
 
 results_list <- list()
 
 for (k in 1:nrow(tree_edges_local)) {
  
  ancestor_label <- as.character(tree_edges_local[k, 1])
  descendant <- tree_edges_local[k, 2]
  
  descendant_label <- if (descendant <= length(tip_labels_local)) {
   tip_labels_local[descendant]
  } else {
   as.character(descendant)
  }
  
  if (!(ancestor_label %in% dimnames(AncestralCoordinatesMAX)[[3]])) {
   stop(
    "Ancestor node `",
    ancestor_label,
    "` was not found in reconstructed ancestral shapes."
   )
  }
  
  ancestor_shape <- AncestralCoordinatesMAX[, , ancestor_label]
  
  max_shape <- if (descendant_label %in% dimnames(predicted_shapesMAX)[[3]]) {
   
   predicted_shapesMAX[, , descendant_label]
   
  } else {
   
   if (!(descendant_label %in% dimnames(AncestralCoordinatesMAX)[[3]])) {
    stop(
     "Descendant node `",
     descendant_label,
     "` was not found in reconstructed ancestral shapes."
    )
   }
   
   AncestralCoordinatesMAX[, , descendant_label]
  }
  
  if (!all(dim(ancestor_shape) == dim(shape_0)) ||
      !all(dim(max_shape) == dim(shape_0))) {
   warning(
    paste(
     "Shape dimension mismatch at edge",
     k,
     "skipping..."
    )
   )
   next
  }
  
  dist_A <- sqrt(sum((ancestor_shape - shape_0)^2))
  dist_B <- sqrt(sum((shape_0 - max_shape)^2))
  
  classification <- classify_heterochrony(
   dist_A = dist_A,
   dist_B = dist_B,
   diff_ratio_threshold = diff_ratio_threshold
  )
  
  results_list[[length(results_list) + 1]] <- data.frame(
   Ancestor = ancestor_label,
   Descendant = descendant_label,
   DistanceA = dist_A,
   DistanceB = dist_B,
   DiffRatio = classification$diff_ratio,
   HeterochronyTrend = classification$trend,
   stringsAsFactors = FALSE
  )
 }
 
 dplyr::bind_rows(results_list)
}


draw_species_sample <- function(
  inds,
  devtimes,
  resampling_proportion,
  min_keep = 2
) {
 
 if (length(inds) == 2) {
  return(seq_along(inds))
 }
 
 n_keep <- max(
  min_keep,
  round(length(inds) * resampling_proportion)
 )
 
 n_keep <- min(n_keep, length(inds))
 
 tochs <- sample(
  seq_along(inds),
  size = n_keep,
  replace = FALSE
 )
 
 return(tochs)
}


prepare_plot_variables <- function(
  results,
  phylo_tree
) {
 
 tip_map <- setNames(
  seq_along(phylo_tree$tip.label),
  phylo_tree$tip.label
 )
 
 results$DescNum <- sapply(results$Descendant, function(x) {
  if (x %in% names(tip_map)) {
   tip_map[x]
  } else {
   as.integer(x)
  }
 })
 
 results$AncNum <- as.integer(results$Ancestor)
 
 results$TrendAbbrev <- ifelse(
  results$HeterochronyTrend == "PERAMORPHOSIS",
  "PERA",
  ifelse(
   results$HeterochronyTrend == "PAEDOMORPHOSIS",
   "PAEDO",
   "NC"
  )
 )
 
 results$TrendColor <- ifelse(
  results$TrendAbbrev == "PERA",
  "blue",
  ifelse(
   results$TrendAbbrev == "PAEDO",
   "red",
   "gray70"
  )
 )
 
 return(results)
}


plot_heterochrony_tree <- function(
  phylo_tree,
  results,
  pdf_name,
  plot_title,
  show_recovery = FALSE
) {
 
 # IMPORTANT:
 # Equal branch lengths are used only here, for plotting.
 # This object is local and does not overwrite the calculation tree.
 phylo_tree_plot <- phylo_tree
 phylo_tree_plot$edge.length <- rep(
  1,
  nrow(phylo_tree_plot$edge)
 )
 
 results_plot <- prepare_plot_variables(
  results,
  phylo_tree_plot
 )
 
 desc_nodes <- phylo_tree_plot$edge[, 2]
 
 edge_label_main <- setNames(
  results_plot$TrendAbbrev,
  results_plot$DescNum
 )[as.character(desc_nodes)]
 
 edge_colors <- setNames(
  results_plot$TrendColor,
  results_plot$DescNum
 )[as.character(desc_nodes)]
 
 if (show_recovery) {
  
  results_plot$LabelCombo <- paste0(
   round(results_plot$DiffRatio, 3),
   " -> ",
   results_plot$RecoveryPercentage,
   "%"
  )
  
  edge_label_sub <- setNames(
   results_plot$LabelCombo,
   results_plot$DescNum
  )[as.character(desc_nodes)]
  
 } else {
  
  edge_label_sub <- setNames(
   round(results_plot$DiffRatio, 3),
   results_plot$DescNum
  )[as.character(desc_nodes)]
 }
 
 edge_label_main[is.na(edge_label_main)] <- ""
 edge_label_sub[is.na(edge_label_sub)] <- ""
 edge_colors[is.na(edge_colors)] <- "gray85"
 
 pdf(
  pdf_name,
  width = 16,
  height = 10
 )
 
 plot(
  phylo_tree_plot,
  show.tip.label = TRUE,
  cex = 1.1,
  main = plot_title
 )
 
 nodelabels(
  cex = 0.9,
  frame = "circle",
  col = "darkblue"
 )
 
 edgelabels(
  text = edge_label_main,
  edge = 1:nrow(phylo_tree_plot$edge),
  col = edge_colors,
  adj = c(0.5, -0.4),
  frame = "none",
  cex = 1.1
 )
 
 edgelabels(
  text = edge_label_sub,
  edge = 1:nrow(phylo_tree_plot$edge),
  col = "black",
  adj = c(0.5, 1.4),
  frame = "none",
  cex = 0.75
 )
 
 legend(
  "topright",
  legend = c(
   "PERAMORPHOSIS",
   "PAEDOMORPHOSIS",
   "NC"
  ),
  col = c(
   "blue",
   "red",
   "gray70"
  ),
  lty = 1,
  lwd = 2,
  cex = 0.9,
  bty = "n"
 )
 
 dev.off()
}


# ============================================================
# ESTIMATE SHARED ORIGIN SHAPE
# ============================================================

gdf_all <- geomorph::geomorph.data.frame(
 Data = CoordinatesData,
 Devtime = DevtimeTwo,
 SpeciesData = SpeciesData
)

fit_all <- geomorph::procD.lm(
 Data ~ Devtime * SpeciesData,
 data = gdf_all
)

predicted_all_origin <- predict(
 fit_all,
 newdata = data.frame(Devtime = 0)
)$mean

shape_0 <- matrix(
 predicted_all_origin,
 ncol = Ndimensions,
 byrow = TRUE
)


# ============================================================
# ESTIMATE MAXIMUM SHAPES
# ============================================================

species_devtime_df <- data.frame(
 Species = SpeciesData,
 DevtimeTwo = DevtimeTwo
)

max_base_devtimes <- species_devtime_df %>%
 dplyr::group_by(Species) %>%
 dplyr::summarise(
  MaxDevtime = max(DevtimeTwo),
  .groups = "drop"
 )

predicted_shapesMAX <- estimate_max_shapes(
 coords = CoordinatesData,
 devtimes = DevtimeTwo,
 species = SpeciesData,
 max_devtimes = max_base_devtimes
)


# ============================================================
# RECONSTRUCT OBSERVED ANCESTRAL SHAPES
# ============================================================

AncestralCoordinatesMAX <- reconstruct_ancestral_shapes(
 predicted_shapes = predicted_shapesMAX,
 phylo_tree = phylo_tree
)


# ============================================================
# OBSERVED HETEROCHRONIC TRENDS
# ============================================================

DistanceResults <- infer_edge_trends(
 shape_0 = shape_0,
 predicted_shapesMAX = predicted_shapesMAX,
 AncestralCoordinatesMAX = AncestralCoordinatesMAX,
 phylo_tree = phylo_tree,
 diff_ratio_threshold = diff_ratio_threshold
)

DistanceResults$DiffRatio <- round(
 DistanceResults$DiffRatio,
 3
)

DistanceResults <- prepare_plot_variables(
 DistanceResults,
 phylo_tree
)

DistanceResults$ObservedTrend <- DistanceResults$HeterochronyTrend


# ============================================================
# OBSERVED TREE PLOT
# ============================================================

plot_heterochrony_tree(
 phylo_tree = phylo_tree,
 results = DistanceResults,
 pdf_name = observed_pdf_name,
 plot_title = observed_plot_title,
 show_recovery = FALSE
)


# ============================================================
# BOOTSTRAP TRENDS
# ============================================================

set.seed(5422)

bootstrap_trends <- list()

species_counts <- table(SpeciesData)

if (any(species_counts < 2)) {
 stop(
  "Some species have fewer than 2 specimens and cannot be used:\n",
  paste(names(species_counts[species_counts < 2]), collapse = ", ")
 )
}

fixed_species <- names(species_counts[species_counts == 2])

attempts <- 0
max_attempts <- max(
 1000,
 n_iter * 20
)

replicate_id <- 0

while (replicate_id < n_iter && attempts < max_attempts) {
 
 attempts <- attempts + 1
 
 cat(
  "Bootstrap attempt:",
  attempts,
  "\n"
 )
 
 coords_boot_list <- list()
 species_boot <- NULL
 devtime_boot <- NULL
 
 sample_failed <- FALSE
 
 for (sp in levels(SpeciesData)) {
  
  inds <- which(SpeciesData == sp)
  
  if (length(inds) < 2) {
   cat(
    " -> Skipping: species",
    sp,
    "has fewer than 2 specimens.\n"
   )
   sample_failed <- TRUE
   break
  }
  
  if (sp %in% fixed_species) {
   
   tochs <- seq_along(inds)
   
  } else {
   
   if (length(inds) < 3) {
    cat(
     " -> Skipping: species",
     sp,
     "has fewer than 3 specimens for resampling.\n"
    )
    sample_failed <- TRUE
    break
   }
   
   tochs <- draw_species_sample(
    inds = inds,
    devtimes = DevtimeTwo,
    resampling_proportion = resampling_proportion,
    min_keep = 2
   )
  }
  
  if (length(tochs) < 2) {
   cat(
    " -> Skipping: fewer than 2 sampled specimens for",
    sp,
    "\n"
   )
   sample_failed <- TRUE
   break
  }
  
  samp_Coords <- CoordinatesData[, , inds[tochs], drop = FALSE]
  samp_LogCS <- DevtimeTwo[inds[tochs]]
  samp_SpeciesData <- as.character(SpeciesData[inds[tochs]])
  
  if (length(unique(samp_LogCS)) < 2) {
   cat(
    " -> Skipping: no size variation retained for",
    sp,
    "\n"
   )
   sample_failed <- TRUE
   break
  }
  
  coords_boot_list[[length(coords_boot_list) + 1]] <- samp_Coords
  devtime_boot <- c(devtime_boot, samp_LogCS)
  species_boot <- c(species_boot, samp_SpeciesData)
 }
 
 if (sample_failed) {
  next
 }
 
 coords_boot <- do.call(
  abind::abind,
  c(
   coords_boot_list,
   list(along = 3)
  )
 )
 
 if (is.null(coords_boot) ||
     length(devtime_boot) != dim(coords_boot)[3]) {
  cat(" -> Skipping: sampling error.\n")
  next
 }
 
 # IMPORTANT:
 # Keep original factor behavior.
 # Do not force levels = tree tips.
 species_boot <- factor(species_boot)
 
 if (length(unique(species_boot)) < 2 || sd(devtime_boot) < 1e-4) {
  cat(" -> Skipping: insufficient diversity or size variation.\n")
  next
 }
 
 if (!all(phylo_tree$tip.label %in% levels(species_boot))) {
  cat(" -> Skipping: not all tree tips are present in this bootstrap replicate.\n")
  next
 }
 
 gdf_boot <- geomorph::geomorph.data.frame(
  coords = coords_boot,
  DevtimeTwo = devtime_boot,
  Species = species_boot
 )
 
 fit_boot_shared <- tryCatch(
  geomorph::procD.lm(
   coords ~ DevtimeTwo * Species,
   data = gdf_boot
  ),
  error = function(e) {
   cat(
    " -> Skipping: procD.lm failed:",
    e$message,
    "\n"
   )
   return(NULL)
  }
 )
 
 if (is.null(fit_boot_shared)) {
  next
 }
 
 shape_0_boot <- tryCatch({
  
  matrix(
   predict(
    fit_boot_shared,
    newdata = data.frame(DevtimeTwo = 0)
   )$mean,
   ncol = Ndimensions,
   byrow = TRUE
  )
  
 }, error = function(e) {
  cat(
   " -> Skipping: origin prediction failed:",
   e$message,
   "\n"
  )
  return(NULL)
 })
 
 if (is.null(shape_0_boot)) {
  next
 }
 
 predicted_shapesMAX_boot <- tryCatch({
  
  estimate_max_shapes(
   coords = coords_boot,
   devtimes = devtime_boot,
   species = species_boot,
   max_devtimes = max_base_devtimes
  )
  
 }, error = function(e) {
  cat(
   " -> Skipping: maximum shape estimation failed:",
   e$message,
   "\n"
  )
  return(NULL)
 })
 
 if (is.null(predicted_shapesMAX_boot) ||
     any(is.na(predicted_shapesMAX_boot))) {
  cat(" -> Skipping: invalid predicted maximum shapes.\n")
  next
 }
 
 # IMPORTANT:
 # Use the original phylo_tree for calculation.
 # Equal branch lengths are NOT used here.
 pca_boot <- tryCatch({
  
  geomorph::gm.prcomp(
   predicted_shapesMAX_boot,
   phy = phylo_tree
  )
  
 }, error = function(e) {
  cat(
   " -> Skipping: PCA failed:",
   e$message,
   "\n"
  )
  return(NULL)
 })
 
 if (is.null(pca_boot) ||
     is.null(pca_boot$ancestors) ||
     nrow(pca_boot$ancestors) == 0) {
  cat(" -> Skipping: no ancestral states recovered.\n")
  next
 }
 
 anc_mat <- as.matrix(pca_boot$ancestors)
 
 AncestralCoordinatesMAX_boot <- array(
  t(anc_mat),
  dim = c(
   Ndimensions,
   Nlandmarks,
   nrow(anc_mat)
  )
 )
 
 AncestralCoordinatesMAX_boot <- aperm(
  AncestralCoordinatesMAX_boot,
  c(2, 1, 3)
 )
 
 dimnames(AncestralCoordinatesMAX_boot)[[3]] <- rownames(anc_mat)
 
 trends_this_boot <- list()
 edge_failed <- FALSE
 
 for (k in 1:nrow(tree_edges)) {
  
  ancestor <- as.character(tree_edges[k, 1])
  descendant <- tree_edges[k, 2]
  
  descendant_label <- if (descendant <= length(tip_labels)) {
   tip_labels[descendant]
  } else {
   as.character(descendant)
  }
  
  if (!(ancestor %in% dimnames(AncestralCoordinatesMAX_boot)[[3]])) {
   edge_failed <- TRUE
   break
  }
  
  ancestor_shape <- AncestralCoordinatesMAX_boot[, , ancestor]
  
  max_shape <- if (descendant_label %in% dimnames(predicted_shapesMAX_boot)[[3]]) {
   
   predicted_shapesMAX_boot[, , descendant_label]
   
  } else {
   
   if (!(descendant_label %in% dimnames(AncestralCoordinatesMAX_boot)[[3]])) {
    edge_failed <- TRUE
    break
   }
   
   AncestralCoordinatesMAX_boot[, , descendant_label]
  }
  
  if (any(is.na(ancestor_shape)) ||
      any(is.na(max_shape)) ||
      !all(dim(ancestor_shape) == dim(shape_0_boot)) ||
      !all(dim(max_shape) == dim(shape_0_boot))) {
   edge_failed <- TRUE
   break
  }
  
  dist_A <- sqrt(sum((ancestor_shape - shape_0_boot)^2))
  dist_B <- sqrt(sum((shape_0_boot - max_shape)^2))
  
  classification <- classify_heterochrony(
   dist_A = dist_A,
   dist_B = dist_B,
   diff_ratio_threshold = diff_ratio_threshold
  )
  
  trends_this_boot[[length(trends_this_boot) + 1]] <- data.frame(
   Ancestor = ancestor,
   Descendant = descendant_label,
   DistanceA = dist_A,
   DistanceB = dist_B,
   DiffRatio = classification$diff_ratio,
   HeterochronyTrend = classification$trend,
   stringsAsFactors = FALSE
  )
 }
 
 if (edge_failed) {
  cat(" -> Skipping: incomplete edge evaluation.\n")
  next
 }
 
 boot_df <- dplyr::bind_rows(trends_this_boot)
 
 if (nrow(boot_df) != nrow(DistanceResults)) {
  cat(" -> Skipping: incomplete bootstrap replicate.\n")
  next
 }
 
 boot_df$Replicate <- replicate_id + 1
 
 bootstrap_trends[[replicate_id + 1]] <- boot_df
 
 replicate_id <- replicate_id + 1
 
 cat(
  " -> Valid replicate saved! Total:",
  replicate_id,
  "\n"
 )
}

cat(
 "Finished with",
 replicate_id,
 "valid replicates after",
 attempts,
 "attempts.\n"
)

if (replicate_id == 0) {
 stop("No valid bootstrap replicates were obtained.")
}

valid_replicates <- replicate_id

if (valid_replicates < n_iter) {
 warning(
  "Only ",
  valid_replicates,
  " valid bootstrap replicates were obtained out of the requested ",
  n_iter,
  ". Recovery percentages will be calculated using valid_replicates as denominator."
 )
}


# ============================================================
# AGGREGATE BOOTSTRAP RESULTS
# ============================================================

all_trends <- dplyr::bind_rows(
 bootstrap_trends,
 .id = "BootstrapListID"
)

all_trends$DiffRatio <- round(
 all_trends$DiffRatio,
 3
)

readr::write_csv(
 all_trends,
 raw_bootstrap_csv_name
)

DistanceResults$BootstrapMatches <- mapply(function(anc, desc, trend) {
 
 sum(
  all_trends$Ancestor == anc &
   all_trends$Descendant == desc &
   all_trends$HeterochronyTrend == trend
 )
 
},
DistanceResults$Ancestor,
DistanceResults$Descendant,
DistanceResults$ObservedTrend
)

DistanceResults$RecoveryPercentage <- round(
 100 * DistanceResults$BootstrapMatches / valid_replicates,
 2
)

DistanceResults <- DistanceResults %>%
 dplyr::select(-BootstrapMatches)

readr::write_csv(
 DistanceResults,
 results_csv_name
)

print(DistanceResults)


# ============================================================
# EXPORT BOOTSTRAP TREND DISTRIBUTION BY BRANCH
# ============================================================

bootstrap_trend_distribution <- all_trends %>%
 dplyr::group_by(
  Ancestor,
  Descendant,
  HeterochronyTrend
 ) %>%
 dplyr::summarise(
  Count = dplyr::n(),
  .groups = "drop"
 ) %>%
 dplyr::group_by(
  Ancestor,
  Descendant
 ) %>%
 dplyr::mutate(
  Percentage = round(
   100 * Count / sum(Count),
   2
  )
 ) %>%
 dplyr::ungroup()


# ============================================================
# FINAL TREE PLOT WITH BOOTSTRAP RECOVERY
# ============================================================

plot_heterochrony_tree(
 phylo_tree = phylo_tree,
 results = DistanceResults,
 pdf_name = bootstrap_pdf_name,
 plot_title = bootstrap_plot_title,
 show_recovery = TRUE
)


# ============================================================
# END
# ============================================================

message("Analysis completed successfully.")
message("Observed plot: ", observed_pdf_name)
message("Bootstrap plot: ", bootstrap_pdf_name)
message("Results table: ", results_csv_name)
message("Raw bootstrap trends: ", raw_bootstrap_csv_name)
