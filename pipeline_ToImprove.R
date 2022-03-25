
# set your working directory to the folder where the files were downloaded
# using setwd()

# execute
# source('installation.R')
# source('functions.R')
source(file = "../CytoQP_toImprove/functions_ToImprove.R")

# create data folder where all analysis will be stored
if(!dir.exists("data")) dir.create("data")

# set it as a main directory
dir <- file.path(getwd(), "data")

# ------------------------------------------------------------------------------
# Bead normalization -----------------------------------------------------------
#-------------------------------------------------------------------------------

# set input directory (pathway to the files that are going to be normalized)
raw_data_dir <- file.path(dir, "RawFiles")

# set a directory where bead-normalized fcs files and plots will be saved
bead_norm_dir <- file.path(dir, "BeadNorm")

# define full pathway to the files that you want to normalize
files <- list.files(raw_data_dir,
                    pattern = ".FCS$",
                    full.names = T)

# create baseline file to which all the files will be normalized
set.seed(2)
ref_sample <- baseline_file(fcs_files = files,
                            beads = "dvs",
                            out_dir = bead_norm_dir)

# Normalize files
bead_normalize(files, cores = 1,
               out_dir = bead_norm_dir,
               non_mass_channel = NULL,
               norm_to_ref = ref_sample,
               to_plot = TRUE,
               remove_beads = TRUE,
               k = 80,
               markers_to_keep = c("CD", "HLA", "IgD", "TCR", "Ir",
                                   "Viability","IL", "IFNa",
                                   "TNF", "TGF", "MIP", "MCP", "Granz"))


# ------------------------------------------------------------------------------
# Visualized files after bead normalization  -----------------------------------
#-------------------------------------------------------------------------------

# Define files for visualization
# Before normalization
raw_data_dir <- file.path(dir, "RawFiles")
files_b <- list.files(raw_data_dir,
                      pattern = ".FCS$",
                      ignore.case = T,
                      full.names = TRUE)

# After normalization
bead_norm_dir <- file.path(dir, "BeadNorm")
files_a <- list.files(bead_norm_dir,
                      pattern = "_beadNorm.fcs$",
                      ignore.case = T,
                      full.names = TRUE)

# Define batch id and sample id for each file
batch_pattern <- stringr::str_match(basename(files_b), "(?i).*(day[0-9]*).*.FCS")[,2]

plot_marker_quantiles(files_after_norm = files_a,
                      files_before_norm = files_b,
                      batch_pattern = batch_pattern,
                      transform_flowframe = TRUE,
                      remove_beads = TRUE,
                      bead_channel = "140",
                      uncommon_prefix = "_beadNorm.fcs|.FCS",
                      markers_to_plot = c("CD", "HLA", "IgD", "IL", "TNF",
                                          "TGF", "GR", "IFNa"),
                      manual_colors = c("darkorchid4", "darkorange", "darkgreen"),
                      out_dir = bead_norm_dir)

# ------------------------------------------------------------------------------
# Signal Cleaning --------------------------------------------------------------
#-------------------------------------------------------------------------------

# set and create the directory where cleaned fcs files will be saved
clean_dir <- file.path(dir, "Cleaned")

# Define which files will be cleaned
files <- list.files(bead_norm_dir,
                    ignore.case = TRUE,
                    pattern = "_beadNorm.fcs$",
                    full.names = TRUE)

# Clean files
clean_files(files, cores = 1,
            out_dir = clean_dir,
            to_plot = "Flagged Only",
            data_type = "MC",
            Segment = 1000,
            transform_flowframe = TRUE,
            non_used_bead_ch = "140")

# ------------------------------------------------------------------------------
# File outliers detection -----------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory
clean_dir <- file.path(dir, "Cleaned")

# Define files for visualization
files <- list.files(clean_dir,
                    pattern = "_cleaned.fcs$",
                    full.names = TRUE)

# Define batch_id for each file
file_batch_id <- stringr::str_match(basename(files),
                                    "(day[0-9]*).*.fcs")[,2]

file_quality_check(fcs_files = files,
                   file_batch_id = file_batch_id,
                   phenotyping_markers = c("Ir","CD", "HLA", "IgD", "Pt"),
                   arcsine_transform = TRUE,
                   nClus = 10,
                   sd = 3)

# ------------------------------------------------------------------------------
# Files debarcoding ------------------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory
clean_dir <- file.path(dir, "Cleaned")

# Define files for debarcoding
files <- list.files(clean_dir,
                    pattern = "_cleaned.fcs$",
                    full.names = TRUE)

# Read in file scores if calculated
file_scores <- readRDS(list.files(path = dir,
                                  recursive = TRUE,
                                  full.names = TRUE,
                                  pattern = "Quality_AOF_score.RDS"))

# Define file batch ID for each file
file_batch_id <- stringr::str_match(basename(files),
                                    "(day[0-9]*).*.fcs")[,2]

# Read in metadata
md <- utils::read.csv(file.path(dir, "RawFiles", "meta_data.csv"))

# read in barcode key
sample_key <- CATALYST::sample_key

# Extract information about barcodes used in each batch
barcodes_list <- list()
for (batch in unique(file_batch_id)){
  idx <- md[md[,"BATCH"] == batch, "BARCODE"]
  barcodes_list[[batch]] <- rownames(sample_key)[idx]
}

# Debarcode files
debarcode_files(fcs_files = files,
                out_dir = NULL,
                file_score = file_scores,
                min_threshold = TRUE,
                barcodes_used = barcodes_list,
                file_batch_id = file_batch_id,
                less_than_th = TRUE,
                barcode_key = sample_key)

# ------------------------------------------------------------------------------
# Files aggregation and file name deconvolution --------------------------------
# ------------------------------------------------------------------------------
# Set input directory
debarcode_dir <- file.path(dir, "Debarcoded")

# Define files for debarcoding
files <- list.files(debarcode_dir,
                    pattern = "_debarcoded.fcs$",
                    full.names = TRUE, recursive = T)

# Define out_dir for aggregated files
aggregate_dir <- file.path(dir, "Aggregated")

# Bring metadata
md <- utils::read.csv(file.path(dir, "RawFiles", "meta_data.csv"))

# Assign barcodes names
md$barcode_name <- paste0(rownames(CATALYST::sample_key)[md$BARCODE])

# Assign new sample names specifying patient id and its batch name
md$fcs_new_name <- paste0(md$ID, "_", md$STIM, "_", md$BATCH, ".fcs")

# Aggregate and deconvolute file names
aggregate_files(fcs_files = files,
                md,
                barcode_column = "barcode_name",
                batch_column = "BATCH",
                cores = 1,
                out_dir = aggregate_dir,
                write_agg_file = TRUE)

# ------------------------------------------------------------------------------
# Files gating -----------------------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
aggregate_dir <- file.path(dir, "Aggregated")

# List files for gating 
files <- list.files(path = aggregate_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)

# Gate the files and plot the gating strategy for each file 
n_plots <- 3  
png(file.path(gate_dir, paste0("gating.png")), 
    width = n_plots * 300, height = length(files) * 300)
layout(matrix(1:(length(files) * n_plots), ncol = n_plots, byrow = TRUE))

for (file in files){
  
  ff <- flowCore::read.FCS(filename = file, 
                           transformation = FALSE)
  
  ff <- gate_intact_cells(flow_frame = ff, 
                          file_name = basename(file), save_gated_flow_frame = FALSE)
  
  ff <- gate_singlet_cells(flow_frame = ff,
                           channels = "Event_length",
                           file_name = basename(file), save_gated_flow_frame = FALSE)
  
  ff <- gate_live_cells(flow_frame = ff, 
                        viability_channel = "Pt195Di", save_gated_flow_frame = TRUE, 
                        file_name = basename(file), suffix = "_gated")
}

dev.off()

# ------------------------------------------------------------------------------
# Normalization using reference sample -----------------------------------------
#-------------------------------------------------------------------------------

# Set input directory
gate_dir <- file.path(dir, "Gated")

# Define reference samples
files_ref <- list.files(gate_dir,
                        pattern = "*_gated.fcs$",
                        full.names = TRUE,
                        recursive = T)

df <- data.frame("file_paths" = files_ref,
                 "batch_labels" = stringr::str_match(files_ref, "day[0-9]*")[,1],
                 "ref_ids" = grepl("REF", files_ref))


model <- train_REF_model(df = df, 
                         markers_to_normalize = c("CD", "HLA", "IgD", 
                                                  "IL", "TN", "MCP", "MIP",
                                                  "Gran", "IFNa"), 
                         arcsine_transform = TRUE,
                         nQ = 2,
                         limit = c(0,8), 
                         quantileValues = c(0.05, 0.95), 
                         goal = "mean",
                         norm_with_clustering = FALSE, 
                         save_model = TRUE, 
                         clustering_markers = c("CD", "HLA", "IgD"))

# Normalize files
normalize_REF(model = model, df = df, arcsine_transform = TRUE, 
              norm_with_clustering = FALSE)

# ------------------------------------------------------------------------------
# Plot batch effect ------------------------------------------------------------
#-------------------------------------------------------------------------------

# Define files before normalization 
gate_dir <- file.path(dir, "Gated")
files_before_norm <- list.files(gate_dir,
                                pattern = ".fcs",
                                full.names = T)

# Define files after normalization 
norm_dir <- file.path(dir, "CytoNormed")
files_after_norm <- list.files(norm_dir,
                               pattern = ".fcs",
                               full.names = T)

# files needs to be in the same order, check and order if needed
test_match_order(x = basename(gsub("Norm_","",files_after_norm)), 
                 basename(files_before_norm))

batch_labels <- stringr::str_match(basename(files_before_norm), "day[0-9]*")[,1]

# Plot batch effect
set.seed(789)
plot_batch(files_before_norm = files_before_norm,
           files_after_norm = files_after_norm,
           batch_labels = batch_labels,
           cores = 1,
           out_dir = norm_dir,
           clustering_markers = c("CD", "IgD", "HLA"),
           manual_colors = c("darkorchid4", "darkorange", "chartreuse4"))

batch_pattern <- "day[0-9]*"
plot_marker_quantiles(files_after_norm = files_after_norm,
                      files_before_norm = files_before_norm,
                      batch_labels = batch_labels,
                      arcsine_transform = TRUE,
                      markers_to_plot = c("CD", "HLA", "IgD", "IL", "TNF",
                                          "TGF", "GR", "IFNa", "MCP", "MIP"),
                      manual_colors = c("darkorchid4", "darkorange", "darkgreen"),
                      out_dir = norm_dir)

# Extract cell frequency and MSI
mx <- extract_pctgs_msi_per_flowsom(files_after_norm = files_after_norm,
                                    files_before_norm = files_before_norm,
                                    nCells = 50000,
                                    phenotyping_markers =  c("CD", "HLA", "IgD"),
                                    functional_markers = c("MIP", "MCP", "IL",
                                                           "IFNa", "TNF", "TGF",
                                                           "Gr"),
                                    xdim = 10,
                                    ydim = 10,
                                    n_metaclusters = 35,
                                    out_dir = norm_dir,
                                    arcsine_transform = TRUE, 
                                    save_matrix = TRUE, 
                                    seed = 343)


# create the list to store the plots
plots <- list()
for (name in names(mx[[1]])){
  df_plot <- prepare_data_for_plotting(frequency_msi_list = mx,
                                       matrix_type = name,
                                       n_neighbours = 11, seed = 1)
  
  
  batch <- stringr::str_match(rownames(df_plot), "day[0-9]*")[,1]
  samples_id <- ifelse(grepl("p1", rownames(df_plot)),"p1",
                      ifelse(grepl("p2", rownames(df_plot)), "p2", "ref"))
  stimulation <- stringr::str_match(rownames(df_plot), "UNS|RSQ|IMQ|LPS|CPG")[,1]
 
  plots[[name]] <- plot_batch_using_freq_msi(df_plot = df_plot, fill = batch, 
                                             shape = samples_id, color = batch,
                                             split_by_normalization = TRUE, title = name)
  
}

gg_a <- grid_arrange_common_legend(plot_lists = plots, nrow = 2, ncol = 2, 
                                   position = "right")

ggplot2::ggsave(filename ="batch_effect_frequency_MSI.png",
                device = "png",
                path = norm_dir,
                plot = gg_a,
                units = "cm",
                width = 22,
                height = 14, dpi = 300)



