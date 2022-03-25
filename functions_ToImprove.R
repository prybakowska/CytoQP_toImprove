#' Finds all the mass/fluorochrome channels for the flow frame
#'
#' @description Finds all the mass channels
#'
#' @param flow_frame Untransformed flow frame
#' @param channels Pattern for non-mass channels, default is
#' "Time|Event_length|Center|Offset|Width|Residual|SSC|FSC|File_scattered"
#' @param ... Additional arguments to pass to grep
#'
#' @return Logical vector with TRUE values for mass channels
#'
#' @export

find_mass_ch <- function(flow_frame,
                         channels = "Time|Event_length|Center|Offset|Width|Residual|SSC|FSC|File_scattered",
                         ...){
  non_mass_ch <- grep(c(channels),
       flowCore::colnames(flow_frame),
       invert = TRUE, ...)
  return(non_mass_ch)
}


' flow_rate_bin_adapted
#'
#' @param x flow frame
#' @param second_fraction the fraction of the seconds used in the data.
#' @param timeCh Time channel.
#' @param timestep
#'
#' @return timeFlowData, the cell assignment to the bins.
#'
#' @references this code is strongly based on flowAI::flow_rate_bin.
.flow_rate_bin_adapted <- function (x, second_fraction = 0.1, timeCh = timeCh, timestep = timestep)
{
  xx <- flowCore::exprs(x)[, timeCh]
  idx <- c(1:nrow(x))
  endsec <- ceiling(timestep * max(xx))
  lenx <- length(xx)
  secbegin <- as.numeric(gsub("(.*)(\\.)(.{0}).*", "\\1\\2\\3", xx[1]))
  tbins <- seq(secbegin, endsec/timestep, by = as.numeric(second_fraction)/timestep)
  if (tail(tbins, n=1) < endsec/timestep){
    tbins <- c(tbins, tail(tbins, n=1) + 10)
  }
  if (secbegin == 0){
    secbegin2 <- 0
  } else {
    secbegin2 <- as.numeric(gsub("(.*)(\\.)(.{1}).*", "\\1\\2\\3", xx[1]/100))
  }

  secbin <- seq(secbegin2, endsec, by = as.numeric(second_fraction))
  minbin <- round(secbin/60, 3)
  nrBins <- length(tbins) - 1
  tbCounts <- c(0, hist(xx, tbins, plot = FALSE)$counts)
  expEv <- lenx/(nrBins)
  binID <- do.call(c, mapply(rep, x = 1:length(tbCounts),
                             times = tbCounts, SIMPLIFY = FALSE))
  if (length(idx) != length(binID))
    stop("length of cell ID not equal length of bin ID")
  timeFlowData <- list(frequencies = cbind(tbins, minbin,
                                           secbin, tbCounts),
                       cellBinID = data.frame(cellID = idx,
                                              binID = binID),
                       info = data.frame(second_fraction = second_fraction,
                                         expFrequency = expEv, bins = nrBins))
  return(timeFlowData)
}


.clean_flow_rate_ind <- function(flow_frame, to_plot = TRUE,
                            out_dir = getwd(), alpha = 0.01, data_type = "MC") {

  if (data_type == "MC"){
    time_division <- 100
    timestep <- 0.01
  }
  else if (data_type == "FC") {
    time_division <- 1
    word <- which(grepl("TIMESTEP", names(flowCore::keyword(flowCore::flowSet(ff)[[1]])),
                        ignore.case = TRUE))
    timestep <- as.numeric(flowCore::keyword(flowCore::flowSet(ff)[[1]])[[word[1]]])
  }
  else{
    stop("type of data MC or FC needs to be specified")
  }

  if (to_plot == "None") {
    to_plot <- FALSE
  }
  else {
    to_plot <- TRUE

  }

  flow_frame@exprs[, "Time"] <- flow_frame@exprs[, "Time"]/time_division


  FlowRateData <- .flow_rate_bin_adapted(flow_frame,
                                         timeCh = "Time",
                                         timestep = timestep)

  FlowRateQC <- .flow_rate_check_adapted(x = flow_frame,
                                        FlowRateData = FlowRateData,
                                        alpha = alpha,
                                        use_decomp = TRUE)

  if (to_plot == TRUE){

    out_dir <- file.path(out_dir, "FlowRateCleaning")
    if(!dir.exists(out_dir)){
      dir.create(out_dir)
    }

    png(file.path(out_dir,
                  gsub(".fcs", "_flowAI.png",
                       basename(flow_frame@description$FILENAME),
                       ignore.case = TRUE)),
        width = 800,
        height = 600)
    if(data_type == "MC"){
      FlowRateQC$res_fr_QC[,1] <- timestep
    }
    p <- .plot_flowrate(FlowRateQC, data_type = data_type)
    print(p)
    dev.off()
  }

  flow_frame_cl <- flow_frame[FlowRateQC$goodCellIDs,]
  flow_frame_cl@exprs[,"Time"] <- flow_frame_cl@exprs[,"Time"]*time_division

  return(flow_frame_cl)
}


#' plot_flowrate
#'
#' @description plots flow rate for .fcs files
#'
#' @param data_type
#' @param FlowRateQC list obtained using flowAI:::flow_rate_check function
#'
#' @return xgraph plot
#'
#' @references this code is adapted from the flowAI:::flow_rate_plot()
#' Monaco, G., Chen, H., Poidinger, M., Chen, J., de Magalhães, J.P.,
#' and Larbi, A. (2016). flowAI: automatic and interactive anomaly discerning
#' tools for flow cytometry data. Bioinformatics 32, 2473–2480.
.plot_flowrate <- function (FlowRateQC, data_type = "MC")
{
  if (data_type == "MC"){
    lab <- "Time (10 * Seconds)"
  } else {
    lab <- "Time (Seconds)"
  }
  second_fraction <- FlowRateQC$res_fr_QC$second_fraction
  num_obs = FlowRateQC$res_fr_QC$num_obs
  frequencies = as.data.frame(FlowRateQC$frequencies)
  anoms = as.data.frame(FlowRateQC$anoms)
  anoms_points = as.data.frame(cbind(sec_anom = frequencies$secbin[anoms$index],
                                     count_anom = anoms$anoms))
  xgraph <- ggplot2::ggplot(frequencies, ggplot2::aes_string(x = "secbin", y = "tbCounts")) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), text = element_text(size = 30)) +
    ggplot2::geom_line(colour = "darkblue")
  xgraph <- xgraph + ggplot2::labs(x = lab, y = paste0("Number of events per 1/",
                                                       1/second_fraction, " of a second"))
  if (!is.null(anoms_points)) {
    xgraph <- xgraph + ggplot2::geom_point(data = anoms_points, aes_string(x = "sec_anom",
                                                                           y = "count_anom"),
                                           color = "green4", size = 5,
                                           shape = 1, stroke = 3)
  }
  return(xgraph)
}


.clean_signal_ind <- function(flow_frame,
                         channels_to_clean = NULL,
                         to_plot = "All",
                         Segment = 1000,
                         out_dir = getwd(),
                         arcsine_transform = TRUE,
                         non_used_bead_ch = NULL,
                         MaxPercCut = 0.5,
                         UseOnlyWorstChannels = TRUE,
                         AllowFlaggedRerun = TRUE,
                         AlwaysClean = TRUE,
                         data_type = "MC",
                         ...){

  channels_to_transform <- find_mass_ch(flow_frame, value = FALSE)

  if (arcsine_transform){

    if(data_type == "MC"){
      ff_t <- flowCore::transform(flow_frame,
                                  flowCore::transformList(flowCore::colnames(flow_frame)[channels_to_transform],
                                                CytoNorm::cytofTransform))
    }
    else if (data_type == "FC"){
      ff_t <- flowCore::transform(flow_frame,
                                  flowCore::transformList(flowCore::colnames(flow_frame)[channels_to_transform],
                                                          flowCore::arcsinhTransform(a = 0, b = 1/150, c = 0)))

    }
    else {
      stop("specify data type MC or FC")
    }

  }
  else {
    ff_t <- flow_frame
  }

  if (!is.null(channels_to_clean)){

    ch_to_clean <- which(flowCore::colnames(flow_frame) %in% channels_to_clean)

    if(!("TIME" %in% toupper(flowCore::colnames(flow_frame)[ch_to_clean]))){
      ind_Time <- grep("TIME", toupper(flowCore::colnames(flow_frame)))
      channels <- unique(sort(c(ch_to_clean, ind_Time)))
    }

  }
  else {

    if (!is.null(non_used_bead_ch)) {
      non_bead_ch <- "140"
    }
    else {
      non_bead_ch <- paste(non_used_bead_ch, collapse="|")
    }

    ind_Time <- grep("TIME", flowCore::colnames(flow_frame), value = T, ignore.case = T)
    ch_to_clean <- c(ind_Time, find_mass_ch(flow_frame, value = TRUE))
    ind_nonbeads <- grep(non_bead_ch, flowCore::colnames(flow_frame), value = TRUE)
    channels <- ch_to_clean[!(ch_to_clean %in% ind_nonbeads)]
    channels <- grep(paste(channels, collapse = "|"), flowCore::colnames(flow_frame))
  }

  out_dir <- file.path(out_dir, "SignalCleaning")
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }

  cleaned_data <- flowCut::flowCut(f = ff_t,
                                   Segment = Segment,
                                   MaxPercCut = MaxPercCut,
                                   Channels = channels,
                                   FileID = gsub("_beadNorm", "_flowCutCleaned",
                                                 basename(flow_frame@description$FILENAME)),
                                   Plot = to_plot,
                                   Directory = out_dir,
                                   UseOnlyWorstChannels = UseOnlyWorstChannels,
                                   AllowFlaggedRerun = AllowFlaggedRerun,
                                   AlwaysClean = AlwaysClean)

  ff_t_clean <- cleaned_data$frame

  if (arcsine_transform){

    if(data_type == "MC"){
      ff_clean <- flowCore::transform(ff_t_clean,
                                  flowCore::transformList(flowCore::colnames(ff_t_clean)[channels_to_transform],
                                                          CytoNorm::cytofTransform.reverse))
    }
    else if (data_type == "FC"){
      ff_clean <- flowCore::transform(ff_t_clean,
                                  flowCore::transformList(flowCore::colnames(ff_t_clean)[channels_to_transform],
                                                          CytoNorm::cytofTransform.reverse(x = 150)))

    }
    else {
      stop("specify data type MC or FC")
    }

  }
  else {
    ff_clean <- ff_t_clean
  }

  return(ff_clean)
}

.save_bead_clean <- function(file,
                             to_plot = "All",
                             clean_flow_rate = TRUE,
                             clean_signal = TRUE,
                             out_dir = getwd(),
                             alpha = 0.01,
                             data_type = "MC",
                             channels_to_clean = NULL,
                             Segment = 1000,
                             arcsine_transform = TRUE,
                             non_used_bead_ch = NULL,
                             MaxPercCut = 0.5,
                             UseOnlyWorstChannels = TRUE,
                             AllowFlaggedRerun = TRUE,
                             AlwaysClean = TRUE,
                             ...){
  # read fcs file
  ff <- flowCore::read.FCS(filename = file,
                           transformation = FALSE)

  if(clean_flow_rate){
    # clean flow rate
    message("cleaning flowrate for ", basename(file))
    ff <- .clean_flow_rate_ind(flow_frame = ff,
                               out_dir = out_dir,
                               to_plot = to_plot,
                               data_type = data_type)

  }

  if(clean_signal){
    # clean signal
    message("cleaning signal for ", basename(file))
    ff <- .clean_signal_ind(flow_frame = ff,
                            to_plot = to_plot,
                            out_dir = out_dir,
                            Segment = Segment,
                            arcsine_transform = arcsine_transform,
                            data_type = data_type,
                            non_used_bead_ch = non_used_bead_ch)
  }



  # Write FCS files
  flowCore::write.FCS(ff,
                      file = file.path(out_dir, gsub("_beadNorm","_cleaned",
                                                       basename(file))))
}

#' Clean flow rate and signal
#'
#' @description Cleans the flow rate using functions from flowAI package and
#' the signal using flowCut package.
#'
#' @param files Character, full path to fcs_file.
#' @param cores Number of cores to be used.
#' @param to_plot Character variable that indicates if plots should be generated.
#' The default is "All", which generates plots for flow rate and all channels.
#' Other options are "Flagged Only", plots the flow rate and channels that were
#' spotted with flowCut as incorrect and "None", does not plots anything.
#' @param clean_flow_rate Logical, if flow rate should be cleaned.
#' @param clean_signal, Logical, if signal should be cleaned.
#' @param out_dir Character, pathway to where the plots should be saved,
#' only if argument to_plot = TRUE, default is set to file.path(getwd(), Cleaned).
#' @param Alpha numeric, as in flowAI::flow_auto_qc. The statistical
#' significance level used to accept anomalies. The default value is 0.01.
#' @param data_type Character, if MC (mass cytometry) of FC (flow cytometry)
#' data are analyzed.
#' @param channels_to_clean Character vector of the channels that needs
#' to be cleaned.
#' @param Segment As in flowCut, an integer value that specifies the
#' number of events in each segment to be analyzed.Default is 1000 events.
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5, default is set to TRUE.
#' @param non_used_bead_ch Character vector, bead channels that does not contain
#' any marker information, thus do not need to be cleaned and used
#' for further analysis.
#' @param MaxPercCut As in flowCut, numeric between 0-1 the maximum percentage of
#' event that will be removed form the data.
#' @param UseOnlyWorstChannels As in flowCut, logical, automated detection of the
#' worst channel that will be used for cleaning.
#' @param AllowFlaggedRerun as in flowCut, logical, specify if flowCut will run
#  second time in case the file was flagged.
#' @param AlwaysClean as in flowCut, logical. The file will be cleaned even if
#' it has a relatively stable signal. The segments that are 7 SD away
#' from the mean of all segments are removed.
#' @param ... Additional arguments to pass to flowcut.
#'
#' @return Cleaned, untransformed flow frame if arcsine_transform argument
#' set to TRUE, otherwise transformed flow frame is returned. Save plots
#' with prefix "_beadNorm_flowAI.png" and "flowCutCleaned.png" to out_dir
#' if parameter to_plot set to "All" or "Flagged Only".
#'
#' @examples
#' # Set and create the directory where cleaned fcs files will be saved
#'clean_dir <- file.path(dir, "Cleaned")
#'
#'# Define which files will be cleaned
#'files <- list.files(bead_norm_dir,
#'                    ignore.case = TRUE,
#'                    pattern = "_beadNorm.fcs$",
#'                    full.names = TRUE)
#'
#'# Clean files
#'clean_files(files, cores = 1,
#'            out_dir = clean_dir,
#'            to_plot = "All",
#'            data_type = "MC",
#'            Segment = 1000,
#'            arcsine_transform = TRUE,
#'            non_used_bead_ch = "140")
#'
#' @import ggplot2
#'
#' @export
clean_files <- function(files,
                        cores = 1,
                        to_plot = "All",
                        clean_flow_rate = TRUE,
                        clean_signal = TRUE,
                        out_dir = NULL,
                        alpha = 0.01,
                        data_type = "MC",
                        channels_to_clean = NULL,
                        Segment = 1000,
                        arcsine_transform = TRUE,
                        non_used_bead_ch = NULL,
                        MaxPercCut = 0.5,
                        UseOnlyWorstChannels = TRUE,
                        AllowFlaggedRerun = TRUE,
                        AlwaysClean = TRUE,
                        ...) {
  # Check parameters
  if(!is(files, "character") & !is(files, "list")) {
    stop("files must be a character vector or a list")
  }

  if (any(!is(cores, "numeric") | cores < 1)){
    stop("cores must be a positive number")
  }

  if(!is(clean_flow_rate, "logical")){
    stop("clean_flow_rate must be logical")
  }

  if(!is(clean_signal, "logical")){
    stop("clean_signal must be logical")
  }

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Cleaned")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

  if (!all(file.exists(files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  # Parallelized analysis
  BiocParallel::bplapply(files, function(x) {
    .save_bead_clean(x,
                     to_plot = to_plot,
                     clean_flow_rate = clean_flow_rate,
                     clean_signal = clean_signal,
                     out_dir = out_dir,
                     alpha = alpha,
                     data_type = data_type,
                     channels_to_clean = channels_to_clean,
                     Segment = Segment,
                     arcsine_transform = arcsine_transform,
                     non_used_bead_ch = non_used_bead_ch,
                     MaxPercCut = MaxPercCut,
                     UseOnlyWorstChannels = UseOnlyWorstChannels,
                     AllowFlaggedRerun = AllowFlaggedRerun,
                     AlwaysClean = AlwaysClean)},
    BPPARAM = BiocParallel::MulticoreParam(workers = cores))

}

#' Creates baseline file for bead normalization
#'
#' @description Creates the reference flow frame for which mean beads
#' values will be computed and used during the normalization.
#'
#' @param fcs_files Character, path to fcs files to be normalized.
#' @param beads Character, as in CATALYST::normCytof, "dvs"
#' (for bead masses 140, 151, 153 ,165, 175)
#' or "beta" (for bead masses 139, 141, 159, 169, 175)
#' or a numeric vector of masses. Default is set to "dvs".
#' @param to_plot Logical, indicates if plots should be generated,
#' default set to FALSE
#' @param out_dir Character, pathway to where the plots should be saved,
#' only if argument to_plot = TRUE, default is set to working directory.
#' @param k The same as in CATALYST::normCytof, integer width of the
#' median window used for bead smoothing (affects visualizations only).
#' @param ncells number of cells to be aggregated per each file, defaults is
#' set to 25000 per file.
#' @param ... Additional arguments to pass to CATALYST::normCytof.
#'
#' @return Returns reference, aggregated flow frame.
#'
#' @examples
#' # set input directory (pathway to the files that are going to be normalized)
#' raw_data_dir <- file.path(dir, "RawFiles")
#'
#' # set a directory where bead-normalized fcs files and plots will be saved
#' bead_norm_dir <- file.path(dir, "BeadNorm")
#'
#' # define full pathway to the files that you want to normalize
#' files <- list.files(raw_data_dir,
#'                     pattern = ".FCS$",
#'                     full.names = TRUE)
#'
#' # create baseline file to which all the files will be normalized
#' set.seed(2)
#' ref_sample <- baseline_file(fcs_files = files,
#'                             beads = "dvs",
#'                             out_dir = bead_norm_dir)
#'
#' @export
baseline_file <- function(fcs_files, beads = "dvs", to_plot = FALSE,
                       out_dir = getw(), k = 80, ncells = 25000, ...){

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "BeadNorm")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  ff <- FlowSOM::AggregateFlowFrames(fileNames = fcs_files,
                            cTotal = length(fcs_files)*ncells)

  dat <- CATALYST::prepData(ff)

  dat_norm <- CATALYST::normCytof(x = dat,
                        beads = beads,
                        remove_beads = TRUE,
                        norm_to = NULL,
                        k = k,
                        plot = to_plot,
                        verbose = FALSE,
                        transform = FALSE,
                        ...)
  ff_ref <- CATALYST::sce2fcs(dat_norm$beads)
  rm(ff)

  if (to_plot == TRUE){

    # plot and save diagnostic plots
    plot_dir <- file.path(out_dir, "Plots_BeadNormalization")
    if(!dir.exists(plot_dir))(dir.create(plot_dir))

    p <- dat_norm$scatter
    ggplot2::ggsave(filename = file.path(plot_dir,"RefBeadGate.png"),
           plot = p, limitsize = FALSE)

    p <- dat_norm$lines
    ggplot2::ggsave(filename = file.path(plot_dir,"RefBeadLines.png"),
           plot = p, limitsize = FALSE)
  }
  return(ff_ref)
}


# Función interna para leer, normalizar y guardar
.bead_normalize_ind <- function(flow_frame,
                           markers_to_keep = NULL,
                           non_mass_channel = NULL,
                           beads = "dvs",
                           norm_to_ref = NULL,
                           remove_beads = TRUE,
                           to_plot = TRUE,
                           out_dir = getwd(),
                           k = 80,
                           ...){

  if (!is.null(markers_to_keep)){

    matches <- paste(markers_to_keep, collapse="|")

    m_to_keep <- grep(matches, FlowSOM::GetMarkers(flow_frame, flowCore::colnames(flow_frame)),
                             ignore.case = TRUE, value = FALSE)

    if(is.null(non_mass_channel)){
      non_mass_ch <- grep("Time|length|Ce140|151|153|165|175|Center|Offset|Width|
                        |Residual|Pd",
                          flowCore::colnames(flow_frame),
                          ignore.case = TRUE, value = FALSE)

    } else {
      matches_ch <- paste(c(non_mass_channel, "Time", "length"), collapse="|")
      non_mass_ch <- grep(matches_ch,
                          flowCore::colnames(flow_frame),
                          ignore.case = TRUE, value = FALSE)

    }


    channels_to_keep <- c(m_to_keep, non_mass_ch)
    channels_to_keep <- flowCore::colnames(flow_frame)[sort(unique(channels_to_keep))]

    flow_frame <- flow_frame[, channels_to_keep]
  }

  if (is.null(norm_to_ref)){
    warning("the reference file is not defined. Each file will be normalized
       to its own bead mean but not across all files")
  }

  dat <- CATALYST::prepData(flow_frame)

  # normalize the data and remove beads
  dat_norm <- CATALYST::normCytof(x = dat,
                                  beads = beads,
                                  remove_beads = remove_beads,
                                  norm_to = norm_to_ref,
                                  k = k,
                                  plot = TRUE,
                                  transform = FALSE,
                                  ...)

  # convert back to .fcs files and save
  f <- CATALYST::sce2fcs(dat_norm$data)

  filename <- basename(flow_frame@description$FILENAME)

  if (to_plot == TRUE){

    # plot and save diagnostic plots
    plot_dir <- file.path(out_dir, "Plots_BeadNormalization")
    if(!dir.exists(plot_dir))(dir.create(plot_dir))

    p <- dat_norm$scatter
    ggplot2::ggsave(filename = file.path(plot_dir,
                                         gsub(".FCS|.fcs","_beadGate.png",
                                              filename)),
           plot = p, limitsize = FALSE)

    p <- dat_norm$lines
    ggplot2::ggsave(filename = file.path(plot_dir,
                                         gsub(".FCS|.fcs","_beadLines.png",
                                              filename)),
           plot = p, limitsize = FALSE)

  }

  f@description$FILENAME <- basename(flow_frame@description$FILENAME)
  f@description$FIL <- basename(flow_frame@description$FILENAME)

  return(f)
}

.save_bead_normalize <- function(file,
                                 markers_to_keep,
                                 beads,
                                 non_mass_channel,
                                 norm_to_ref,
                                 remove_beads,
                                 to_plot,
                                 out_dir,
                                 k,
                                 ...) {

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "BeadNorm")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  # read fcs file
  ff <- flowCore::read.FCS(file, transformation = FALSE,
                           truncate_max_range = FALSE)

  # bead normalize the files
  ff_norm <- .bead_normalize_ind(flow_frame = ff,
                                 out_dir = out_dir,
                                 non_mass_channel = non_mass_channel,
                                 norm_to_ref = norm_to_ref,
                                 remove_beads = remove_beads,
                                 to_plot = to_plot,
                                 k = k,
                                 markers_to_keep = markers_to_keep)

  # save normalized FCS files
  flowCore::write.FCS(ff_norm, filename = file.path(out_dir,
                                                    gsub(".FCS","_beadNorm.fcs",
                                                         basename(file),
                                                         ignore.case = TRUE)))
}
#' Bead-based normalization
#'
#' @description Performs bead-based normalization using beads spiked in
#' the sample. It is based on functions from CATALYST package.
#'
#' @param files Character vector or list with the paths of the raw files.
#' @param cores Number of cores to be used. Works only for not-Widows users.
#' @param markers_to_keep Character vector, marker names to be kept after
#' the normalization, can be full marker name e.g. "CD45" or "CD".
#' If NULL (default) all markers will be normalized and kept in flowframe.
#' Selection of the markers will reduce file volume and speedup the analysis.
#' Non-mass channels like Time, Event_length, Gaussian parameters and in addition
#' palladium barcoding channels are kept if non_mass_ch set to NULL.
#' @param non_mass_channel Character vector, non-mass channels to keep for
#' further analysis. Can be full channel name like Eu151Di or 151.
#' By default "Time" and "event_length" will be always kept in the flow frame.
#' @param beads Character, as in CATALYST::normCytof, "dvs"
#' (for bead masses 140, 151, 153 ,165, 175)
#' or "beta" (for bead masses 139, 141, 159, 169, 175)
#' or a numeric vector of masses. Default is set to "dvs".
#' @param norm_to_ref flow frame, created by baseline_file function to which
#' input data will be normalized, default is set to NULL.
#' @param to_plot Logical if to plot bead gate and bead normalization lines
#' for each file.Defaults is set to TRUE.
#' @param out_dir Character, pathway to where the bead normalized fcs files
#' and plots should be saved, for plots only if argument to_plot = TRUE,
#' default is set to file.path(getwd(), BeadNorm).
#' @param k The same as in CATALYST::normCytof, integer width of the
#' median window used for bead smoothing (affects visualizations only).
#' @param remove_beads Logical, as in CATALYST::normCytof if beads should be
#' removed from fcs files. Default set to TRUE. Note, should be set to FALSE if
#' none of the channels is beads-specific.
#' @param ... Additional arguments to pass to normCytof.
#'
#' @return Save bead-normalized fcs files and plots to out_dir.
#'
#' @examples
#' # set input directory (pathway to the files that are going to be normalized)
#' raw_data_dir <- file.path(dir, "RawFiles")
#'
#' # set a directory where bead-normalized fcs files and plots will be saved
#' bead_norm_dir <- file.path(dir, "BeadNorm")
#'
#' # define full pathway to the files that you want to normalize
#' files <- list.files(raw_data_dir,
#'                     pattern = ".FCS$",
#'                     full.names = TRUE)
#'
#' # create baseline file to which all the files will be normalized
#' set.seed(2)
#' ref_sample <- baseline_file(fcs_files = files,
#'                             beads = "dvs",
#'                             out_dir = bead_norm_dir)
#'
#' # Normalize files
#' bead_normalize(files, cores = 1,
#'                out_dir = bead_norm_dir,
#'                non_mass_channel = NULL,
#'                norm_to_ref = ref_sample,
#'                to_plot = TRUE,
#'                remove_beads = TRUE,
#'                k = 80,
#'                markers_to_keep = c("CD", "HLA", "IgD", "TCR", "Ir",
#'                                  "Viability","IL", "IFNa",
#'                                    "TNF", "TGF", "MIP", "MCP", "Granz"))
#'
#' @export
bead_normalize <- function(files,
                           cores = 1,
                           markers_to_keep = NULL,
                           non_mass_channel = NULL,
                           beads = "dvs",
                           norm_to_ref = NULL,
                           remove_beads = TRUE,
                           to_plot = TRUE,
                           out_dir = NULL,
                           k = 80,
                           ...){
  # Check parameters
  if(!is(files, "character") & !is(files, "list")) {
    stop("files must be a character vector or a list")
  }

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  if (any(!is(cores, "numeric") | cores < 1)){
    stop("cores must be a positive number")
  }



  # Parallelized analysis
  BiocParallel::bplapply(files, function(x) {
    .save_bead_normalize(x,
                         markers_to_keep,
                         non_mass_channel,
                         beads,
                         norm_to_ref,
                         remove_beads,
                         to_plot,
                         out_dir,
                         k)},
    BPPARAM = BiocParallel::MulticoreParam(workers = cores))

}


#' gate_out_beads
#'
#' @description removes beads from the files that contain them e.g files before
#' bead normalization
#'
#' @param bead_channel character, the mass for bead channel that is exclusively used for
#' beads identification, no marker is present at this channel, default 140.
#' @param flow_frame Flow frame, if unstransformed arcsine_transform should be
#' kept as default, TRUE.
#'
#' @return flow frame with beads removed
.gate_out_beads <- function(bead_channel,
                           flow_frame){
  ch <- grep(pattern = bead_channel, x = flowCore::colnames(flow_frame), value = TRUE)
  ids <- flow_frame[,ch] > 0
  # calculate threshold
  th <- flowDensity::deGate(obj = flow_frame[ids,], channel = ch)
  #remove beads
  cells_to_remove <- flow_frame@exprs[, ch] < th
  flow_frame <- flow_frame[cells_to_remove,]
  return(flow_frame)
}

#' Plots quantiles for the markers
#'
#' @description Calculates quantiles (0.01, 0.25, 0.5, 0.75, 0.99) for
#' selected markers and plots them as diagnostic plots.
#'
#' @param files_before_norm Character, full path to the unnormalized fcs_files.
#' @param files_after_norm Character, full path to the normalized fcs_files.
#' @param batch_pattern Character, batch pattern to be match in the fcs file name
#' @param uncommon_prefix Character vector or string, uncommon prefix in
#' the basename of the fcs files. The file names need to match, so uncommon prefix
#' needs to be removed. If NULL (default) prefix like
#' "Norm|_CC_gated.fcs|_gated.fcs|_beadNorm.fcs|.FCS|.fcs" will be removed.
#' Default is set to NULL.
#' @param bead_channel character, the mass for bead channel that is exclusively
#' used for beads identification (no marker is assign to this channel),
#' Default 140.
#' @param remove_beads Logical, if beads needs to be removed. This needs to be
#' set to TRUE if files contain beads e.g before beads normalization,
#' default is set to TRUE. For the visualization purpose the beads will be
#' removed using channel set in bead_channel.
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5.
#' @param markers_to_plot character vector, marker names to be plotted, can be
#' full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted
#' @param manual_colors character, vector of the colors to be used,
#' the number of colors needs to be equal to the length of batch_pattern
#' @param out_dir Character, pathway to where the plots should be saved,
#' default is set to working directory.
#' @param transform_list Transformation list to pass to the flowCore
#' transform function, see flowCore::transformList(), if different transformation
#' than arcsine is needed. Only if arcsine_transform is FALSE. If NULL and
#' arcsine_transform = FALSE no transformation will be applied.Default set to NULL.
#'
#' @return Save the pdf with plots to out_dir.
#'
#' @examples
#' # Define files for visualization
#' # Before normalization
#' raw_data_dir <- file.path(dir, "RawFiles")
#' files_b <- list.files(raw_data_dir,
#'                       pattern = ".FCS$",
#'                       ignore.case = T,
#'                       full.names = TRUE)
#'
#' # After normalization
#' bead_norm_dir <- file.path(dir, "BeadNorm")
#' files_a <- list.files(bead_norm_dir,
#'                       pattern = "_beadNorm.fcs$",
#'                       ignore.case = T,
#'                       full.names = TRUE)
#'
#' # Define batch id and sample id for each file
#' batch_pattern <- stringr::str_match(basename(files_b), "(?i).*(day[0-9]*).*.FCS")[,2]
#'
#' plot_marker_quantiles(files_after_norm = files_a,
#'                      files_before_norm = files_b,
#'                      batch_pattern = batch_pattern,
#'                      arcsine_transform = TRUE,
#'                      remove_beads = TRUE,
#'                      bead_channel = "140",
#'                      uncommon_prefix = "_beadNorm.fcs|.FCS",
#'                      markers_to_plot = c("CD", "HLA", "IgD", "IL", "TNF",
#'                                          "TGF", "GR", "IFNa"),
#'                      manual_colors = c("darkorchid4", "darkorange", "darkgreen"),
#'                      out_dir = bead_norm_dir)
#'
#' @import ggplot2
#'
#' @importFrom (magrittr,"%>%")
#'
#' @export
plot_marker_quantiles <- function(files_before_norm,
                                  files_after_norm,
                                  batch_pattern = NULL,
                                  batch_labels = NULL,
                                  remove_beads = FALSE,
                                  bead_channel = "140",
                                  uncommon_prefix = NULL,
                                  arcsine_transform = TRUE,
                                  markers_to_plot = NULL,
                                  plot_name = "Marker_distribution_across_aliquots_and_batches",
                                  manual_colors = NULL,
                                  out_dir = NULL,
                                  transform_list = NULL){

  if(!(length(files_after_norm) == length(files_before_norm))){
    stop("files_before and after does not have the same length")
  }

  fcs_files <- c(files_after_norm, files_before_norm)
  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
                   basename(files_before_norm))

  if(is.null(batch_labels) & is.null(batch_pattern)){
    stop("define batch_labels or batch_pattern")
  } else if (!(is.null(batch_labels)) & !(is.null(batch_pattern))){
    stop("both batch_labels and batch_pattern are defined, desellect one option by  setting to NULL")
  }


  tmp <- c(paste0(files_after_norm, "_YES"), paste0(files_before_norm, "_NO"))

  if(!is.null(batch_labels)){
    if(length(tmp) != length(rep(batch_labels, 2))){
      stop("The lenght of batch labels is not equal to the lenght of files")
    }
  }

  ff_tmp <- flowCore::read.FCS(file.path(fcs_files[1]))

  if (!is.null(markers_to_plot)){

    if(!is.character(markers_to_plot)){
      stop ("markers are not a character vector")
    }

    matches <- paste(markers_to_plot, collapse="|")

    norm_markers <- grep(matches,
                         FlowSOM::GetMarkers(ff_tmp, find_mass_ch(ff_tmp,
                                                          value = TRUE)),
                         value = TRUE, ignore.case = FALSE)
  } else {
    norm_markers <- find_mass_ch(ff_tmp, value = TRUE)
    norm_markers <- FlowSOM::GetMarkers(ff_tmp, norm_markers)
  }


  quantile_values <-  c(0.01, 0.25, 0.5, 0.75, 0.99)
  quantiles <- expand.grid(File = tmp,
                           Marker = norm_markers,
                           Quantile = quantile_values,
                           Value = NA)

  if(!is.null(batch_pattern)){
    quantiles <- cbind(quantiles, "Batch" = stringr::str_match(
      basename(as.character(quantiles$File)), batch_pattern)[,1])
  }

  if(!is.null(batch_labels)){
    quantiles <- cbind(quantiles, "Batch" = rep(batch_labels, 2))
  }

  quantiles$Normalization <- gsub(".*.fcs_|.*.FCS_", "", quantiles$File)
  quantiles$File <- gsub("_YES|_NO", "", quantiles$File)

  for (file in fcs_files) {
    print(file)

    ff <- flowCore::read.FCS(file, transformation = FALSE)

    if(arcsine_transform){
      ff <- flowCore::transform(ff,
                                flowCore::transformList(grep("Di", flowCore::colnames(ff),
                                                             value = TRUE),
                                                        CytoNorm::cytofTransform))
    } else if (!is.null(transform_list)){
      ff <- flowCore::transform(ff, transform_list)
    } else {
      ff <- ff
    }

    norm <- quantiles$Normalization[(which(quantiles$File == file)[1])]

    if(norm == "NO" & remove_beads){
      ff <- .gate_out_beads(bead_channel = bead_channel, flow_frame = ff)
    }

    for (marker in names(norm_markers)) {
      quantiles_res <-stats::quantile(flowCore::exprs(ff)[, marker],
                                quantile_values)
      for (i in seq_along(quantiles_res)) {
        quantiles <- quantiles %>%
          dplyr::mutate(Value = replace(Value,
                                        File == file &
                                          Marker == norm_markers[marker] &
                                          Quantile == quantile_values[i],
                                        quantiles_res[i]))

      }
    }
  }

  if(is.null(uncommon_prefix)){
    quantiles$Sample <- gsub(pattern = "Norm_", replacement = "",
                             ignore.case = TRUE,
                             x = gsub(
                               pattern = "_CC_gated.fcs|_gated.fcs|_beadNorm.fcs|.FCS|.fcs",
                               replacement = "", ignore.case = TRUE,
                               x = basename(as.character(quantiles$File))))
  }
  else {
    uncommon_prefix <- paste(uncommon_prefix, collapse = ("|"))
    quantiles$Sample <- gsub(pattern = "Norm_", replacement = "",
                             ignore.case = TRUE,
                             x = gsub(pattern = uncommon_prefix,
                                      replacement = "",
                                      ignore.case = TRUE,
                                      x =  basename(as.character(
                                        quantiles$File))))
  }

  ncols <- length(unique(quantiles$Batch))
  p <- quantiles %>% dplyr::filter(Normalization == "YES") %>%
    ggplot2::ggplot(aes(x = Sample,
               y = Value,
               color = Batch)) +
    geom_point(data = quantiles %>% dplyr::filter(Normalization == "NO"),
               aes(alpha = ifelse(Quantile == "0.5", 2, 0)), color = "grey31") +
    geom_line(data = quantiles %>% dplyr::filter(Normalization == "NO"),
              aes(alpha = ifelse(Quantile != "0.5", 2, 0),
                  group = interaction(Batch, Quantile),
                  size = ifelse(Quantile %in% c("0.05", "0.95"), 1,
                                ifelse(Quantile == "0.5", 0, 2))),
              alpha = 0.5, color = "grey31") +
    ylab(label = levels(quantiles$Marker))+
    geom_point(aes(alpha = ifelse(Quantile == "0.5", 2, 0))) +
    geom_line(aes(alpha = ifelse(Quantile != "0.5", 2, 0),
                  group = interaction(Batch, Quantile),
                  size = ifelse(Quantile %in% c("0.05", "0.95"), 1,
                                ifelse(Quantile == "0.5", 0, 2))),
              alpha = 0.5) +
    scale_size_identity() +
    scale_alpha_identity() +
    facet_wrap(~ Marker + Batch, ncol = ncols, scales = "free_x") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")

  if (!is.null(manual_colors)){
    p <- p + ggplot2::scale_colour_manual(values = c(manual_colors))
  }

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- getwd()
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}


  ggplot2::ggsave(filename = paste0(plot_name, ".pdf"),
                  plot = p,
                  path = out_dir,
                  width = length(fcs_files)*0.25,
                  height = length(norm_markers)*4, limitsize = FALSE)
}

#' Prepares FlowSOM
#'
#' @description Builds FlowSOM tree
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param phenotyping_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted
#' @param nCells Numeric, the total number of cells, to use for FlowSOM clustering.
#' This number is determined by total number of fcs files, as a defult 1000 cells
#' is used per file
#' @param xdim Numeric, parameter to pass to FlowSOM, width of the SOM grid
#' @param ydim Numeric, parameter to pass to FlowSOM, geight of the SOM grid
#' @param nClus Numeric, exact number of clusters for metaclustering
#' @param out_dir Character, pathway to where the FlowSOM clustering plot should
#'  be saved, default is set to working directory.
#' @param batch Character, Character, aqusition batch for each fcs file..
#' Pass to FlowSOM plot name, defult is set to NULL
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5. Default set to TRUE.
#' @param seed numeric, set to obtain reproducible results, default 1
#' @param transform_list Transformation list to pass to the flowCore
#' transform function, see flowCore::transformList, if different transformation
#' than arcsine is needed. Only if arcsine_transform is FALSE. If NULL and
#' arcsine_transform = FALSE no transformation will be applied.
#' @param my_colors An array specifying colors to be used for the background
#' coloring of metaclusters in FlowSOM and t-SNE plot. Must have a length equal
#' to the nClus.
#'
#' @param to_plot Logical, if FlowSOM tree and t-SNE map should be plotted,
#' default set to TRUE.
#'
#' @return fsom object
#'
#' @export
fsom_aof <- function(fcs_files,
                     phenotyping_markers,
                     nCells = length(fcs_files)*10000,
                     xdim = 10,
                     ydim = 10,
                     nClus = 10,
                     out_dir = NULL,
                     batch = NULL,
                     arcsine_transform = TRUE,
                     transform_list = NULL,
                     my_colors = NULL,
                     seed = 1,
                     to_plot = TRUE){


  if(!exists("phenotyping_channels")){
    o <- capture.output(ff_tmp <- flowCore::read.FCS(file.path(files[1])))
    markers <- FlowSOM::GetMarkers(ff_tmp, flowCore::colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers,
                                       collapse = ("|")), markers, value = TRUE)
  }

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Quality_Control")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  if(arcsine_transform){
    trans <- flowCore::transformList(names(phenotyping_channels),
                                     CytoNorm::cytofTransform)
  }
  else {
    if(is.null(transform_list)){
      stop("transform_list must be defined")
    }
    trans <- transform_list
  }

  fsom <- CytoNorm::prepareFlowSOM(file = fcs_files,
                                   colsToUse = names(phenotyping_channels),
                                   seed = seed,
                                   nCells = nCells,
                                   transformList = trans,
                                   FlowSOM.params = list(xdim = xdim,
                                                         ydim = ydim,
                                                         nClus = nClus,
                                                         scale = FALSE))

  if(to_plot){
    if(is.null(my_colors)){
      backgroundColors <- NULL
    }
    else {

      if(max(as.numeric(fsom$metaclustering)) < length(my_colors)){
        warning("The number of colors is greater than the number of metaclusters only
             needed number of colors will be used")
        nmcl <- max(as.numeric(fsom$metaclustering))
        backgroundColors <- my_colors[1:nmcl]

      }

      if(max(as.numeric(fsom$metaclustering)) > length(my_colors)){
        warning("The number of colors is lower than the number of metaclusters
              default colors will be used")
        backgroundColors <- NULL

      }

      if(max(as.numeric(fsom$metaclustering)) == length(my_colors)){
        backgroundColors <- my_colors
      }
    }

    if(!is.null(batch)){
      filename <- paste0(batch, "_FlowSOM_clustering.pdf")
    }
    else {
      filename <- "FlowSOM_clustering.pdf"
    }

    fsomPlot <- FlowSOM::PlotStars(fsom = fsom,
                                   title = "FlowSOM clustering",
                                   backgroundValues = fsom$metaclustering,
                                   maxNodeSize = 3,
                                   backgroundColors = backgroundColors)
    fsomTsne <- FlowSOM::PlotDimRed(fsom = fsom, plotFile = NULL, seed = seed, cTotal = 20000,
                                    title = "tSNE visualization of FlowSOM metaclusters")

    figure <- ggpubr::ggarrange(fsomPlot, fsomTsne,
                                # labels = c("FlowSOM clustering", "tsne"),
                                ncol = 2, nrow = 1)

    ggplot2::ggsave(filename = filename, plot = figure, device = "pdf", path = out_dir,
                    width =24, height = 10)
  }

  return(fsom)
}

#' Calculates scaled aof scores
#'
#' @description Calculates scaled AOF scores and sample quality score.
#' Additionally plots heatmaps for both raw  aof scores and scaled aof scores.
#' @param aof_scores Matrix, array, Aof scores obtained using function
#' cytutils::greedyCytometryAof.
#' @param out_dir Character, pathway to where the plots and scores should be
#' saved, if NULL file.path(getwd(), "Quality_Control) will be used.
#' @param aof_channels Character vector with the markers and their corresponding
#' channel names used for aof_scoring. Used only for plotting markers instead
#' of channels. If NULL (default) the colnames(aof_scores) will be used.
#' @param batch Character, acquisition batch for each fcs file.
#' Pass to FlowSOM plot name, default is set to NULL
#' be saved, default is set to working directory.
#'
#' @return Returns data frame with sample scores.Saved heatmaps and data frames
#' for AOF scores and AOF scaled scores.
#'
#' @export
scaled_aof_score <- function(aof_scores, out_dir = NULL, aof_channels = NULL,
                             batch = NULL){
  aof_scores_scaled <- scale(aof_scores)
  aof_scores_scaled <- pmax(aof_scores_scaled, 0)^2
  sample_scores <- apply(aof_scores_scaled, 1, sum, na.rm = TRUE)

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Quality_Control")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  df <- as.data.frame(sample_scores)

  list_scores <- list("aof_scores_per_marker" = aof_scores,
                      "scaled_AOF" = aof_scores_scaled)

  for (name in names(list_scores)) {

    if(!is.null(batch)){
      filename <- file.path(out_dir, paste0(batch, "_", name, ".pdf"))
      main <- paste0(batch, "_", name)
    } else {
      filename <- file.path(out_dir, paste0(name, ".pdf"))
      main <- name
    }

    if(is.null(aof_channels)){
      phenotyping_channels <- colnames(list_scores[[name]])
    }
    else if (all(colnames(aof_scores) %in% names(aof_channels))){

      phenotyping_channels <- aof_channels[colnames(aof_scores)]
    } else {
      phenotyping_channels <- colnames(list_scores[[name]])
      warning("aof_channels do not correspond to the channels selected for
              AOF scoring, colnames from aof_scores will be used for plotting")
    }

    pheatmap::pheatmap(list_scores[[name]],
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       color = colorRampPalette(
                         RColorBrewer::brewer.pal(n = 9,
                                                  name = "YlGnBu"))(100),
                       display_numbers = TRUE,
                       labels_col = phenotyping_channels,
                       labels_row = basename(rownames(list_scores[[name]])),
                       filename = filename,
                       main = main,
                       number_format = "%.1f",
                       fontsize_number = 8,
                       number_color = "black",
                       width = 10)
  }

  if(is.null(batch)){
    saveRDS(list_scores, file.path(out_dir, "AOF_scores_and_Scaled_AOF_scores.RDS"))
  }
  else {
    saveRDS(list_scores, file.path(out_dir,
                                   paste0(batch, "_AOF_scores_and_Scaled_AOF_scores.RDS")))
  }
  return(df)
}

#' Calculates AOF scores and scaled AOF scores
#'
#' @description  Calculates AOF (Average Overlap Frequency) scores using flowSOM
#' object and greedy algorithm from cytutils package.
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param phenotyping_markers Character vector, marker names to be used for
#' clustering, can be full marker name e.g. "CD45" or "CD" if all CD-markers
#' needs to be plotted.
#' @param fsom FlowSOM object as generated by fsom_aof
#' @param out_dir Character, pathway to where the AOF scores and plots should
#' be saved, default is set to file.path(getwd(), "Quality_Control")
#' @param batch Character, acquisition batch for each fcs file.
#' This argument is passed to AOF plot names, default is set to NULL.
#'
#' @return Returns data frame with the scaled AOF scores and heatmap plots
#' representing AOF scores and scaled AOF scores.
#'
#' @export
aof_scoring <- function(fcs_files,
                        phenotyping_markers,
                        fsom,
                        out_dir = NULL,
                        batch = NULL){

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Quality_Control")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  if(!exists("phenotyping_channels")){

    ff_tmp <- flowCore::read.FCS(file.path(files[1]))
    markers <- FlowSOM::GetMarkers(ff_tmp, flowCore::colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers,
                                       collapse = ("|")), markers, value = TRUE)

    if(length(grep("Ir", phenotyping_channels)) > 1){
      phenotyping_channels <- phenotyping_channels[-(grep("Ir",
                                                          phenotyping_channels)[2])]
    }
  }

  aof_scores <- lapply(fcs_files, function(file) {
    print(paste("calculating AOF", file))
    File_ID <- which(fcs_files == file)
    idx <- which(fsom$data[,"File"] == File_ID)
    fcs_data <- fsom$data[idx,]
    MC <- fsom$metaclustering[fsom$map$mapping[idx, 1]]

    aof_tmp <- cytutils::greedyCytometryAof(fcs_data = fcs_data,
                                            y = MC,
                                            channel_names = names(phenotyping_channels),
                                            width = 0.05,
                                            cofactor = 5,
                                            verbose = TRUE)
    return(aof_tmp$Aof)
  })

  aof_scores <- do.call("rbind", aof_scores)
  rownames(aof_scores) <- fcs_files
  colnames(aof_scores) <- names(phenotyping_channels)

  scaled_aof_score(aof_scores = aof_scores,
                   out_dir = out_dir,
                   aof_channels = phenotyping_channels,
                   batch = batch)
}

#' Detects outliers based on sample quality scores
#'
#' @description Detects outlier files based on sample quality scores,
#' generates plot for outlier and .csv file which indicates which fcs files
#' could be discarded from further analysis.
#'
#' @param scores List of scaled scores per acquisition batch or data frame
#' of scaled scores, both generated by scaled_aof_score or aof_scoring function
#' @param out_dir Character, pathway to where the plot and .csv files with
#' quality scores should be saved, default is set to NULL, thus
#' file.path(getwd(), "Quality_Control" will be generated.
#' @param sd How many standard deviation should be use to detect outliers
#' default is set to 3.
#'
#' @return Save .RDS and .csv Quality scores for further analysis, and
#' plots and save .png for Quality AOF scores for all files.
#'
#' @export
file_outlier_detecion <- function(scores, out_dir = NULL, sd) {

  if(!inherits(scores, "data.frame") & !inherits(scores, "list")){
    stop("df scores are neither data frame nor list of the data frames")
  }

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Quality_Control")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  if(inherits(scores, "list")){
    df_scores <- do.call(rbind, scores)
  }
  else {
    df_scores <- scores
  }

  df_scores$file_names <- basename(rownames(df_scores))

  scores_median <- stats::median(df_scores$sample_scores)
  scores_MAD <- stats::mad(df_scores$sample_scores)

  df_scores$quality <- ifelse(df_scores$sample_scores >
                                (scores_median + sd * scores_MAD),"bad","good")

  bad_scores <- sum(df_scores$quality == "bad")

  colors <- c("bad" = "red", "good" = "darkblue", "threshold= " = "orange")

  max_score <- max(df_scores$sample_scores)
  max_pctgs <- max_score + (max_score * 0.1)

  p <- ggplot2::ggplot(df_scores, aes(x = file_names, y = sample_scores,
                                      color = quality)) +
    geom_point(size = 4) +
    scale_colour_manual(values = colors) +
    ylim(-0.5, max_pctgs) +
    annotate(geom="text", x = mean(as.numeric(as.factor(df_scores$file_names))),
             y= max_score - 0.05*max_score, label=paste("N bad = ", bad_scores),
             color="red", size = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text = element_text(size = 11),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          # panel.border = fill = "black",
          axis.line = element_line(colour = "black"),
          legend.position = "none") +
    scale_x_discrete(breaks = df_scores$file_names[df_scores$quality == "bad"])

  if(scores_median + sd * scores_MAD <= max_pctgs){
    p + geom_hline(yintercept = scores_median + sd * scores_MAD,
                   linetype = "dashed", color = "darkgreen", size = 1)
  }

  ggplot2::ggsave(filename = "Quality_AOF_score.png", plot = p,
                  path = file.path(out_dir))

  saveRDS(df_scores, file.path(out_dir, "Quality_AOF_score.RDS"))
  write.csv(df_scores, file = file.path(out_dir, "Quality_AOF_score.csv"))
}

#' Check the quality of acquired files
#'
#' @description Wrapper function to perform sample quality scoring.
#' First, it clusters the data per each batch (if batch argument is defined) and
#' calculates the AOF scores and Quality scores per batch using AOF algorithm.
#' Next, based on Quality scores, it detects outliers across all the files,
#' regarding the batch.
#'
#' @param fcs_files Character, full path to fcs files.
#' @param file_batch_id Character vector with batch label for each fcs_file,
#' the order and the length needs to be the same as in fcs_files. if only batch
#' is processed can be prepared as e.g. file_batch_id <- rep("batch", length(files))
#' @param out_dir Character, pathway to where the plots should be saved,
#' default is set to NULL, which means that the following path will be created
#' file.path(getwd(), "Quality_Control").
#' @param phenotyping_markers Character vector, marker names to be used for
#' flowsom clustering including DNA marker Iridium and viability staining
#' if available. Can be full marker name e.g. "CD45" or pattern "CD" if
#' all CD-markers needs to be plotted. Default is set to NULL, thus all the mass
#' channels will be used.
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5. Default is set to TRUE.
#' If FALSE, the transform list to pass to the flowCore transform function must
#' be defined and pass as an additional argument to fsom_aof function.
#' @param sd Numeric, number of standard deviation allowed for file outlier
#' detection, default = 3.
#' @param nClus Numeric, as in FlowSOM, number of metaclusters to be obtained
#' @param ... Arguments to be passed to fsom_aof function for FlowSOM parameter
#' adjustment and plotting: xdim, ydim, transform_list, my_colors, seed, to_plot.
#'
#' @return Plots Quality AOF scores for all files and save .RDS and .csv Quality
#' scores for further analysis, files are saved in out_dir.
#'
#' @importFrom flowCore transformList
#' @import ggplot2
#'
#' @export
file_quality_check <- function(fcs_files,
                               file_batch_id = NULL,
                               out_dir = NULL,
                               phenotyping_markers = NULL,
                               arcsine_transform = TRUE,
                               sd = 3,
                               nClus = 10,
                               ...){

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Quality_Control")
  }

  if(length(file_batch_id) != length(fcs_files)){
    stop("the lenght of the file_batch_id is not equal to the lenght of fcs_files")
  }

  if(!dir.exists(out_dir)){dir.create(out_dir)}

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  if (!is.null(file_batch_id)) {
    scores <- lapply(unique(file_batch_id), function(batch) {
      print(batch)

      files <- fcs_files[file_batch_id == batch]
      fsom <- fsom_aof(fcs_files = files,
                       phenotyping_markers = phenotyping_markers,
                       out_dir = out_dir,
                       arcsine_transform = arcsine_transform,
                       nClus = nClus,
                       batch = batch, ...)
      return(aof_scoring(fcs_files = files,
                         phenotyping_markers = phenotyping_markers,
                         fsom = fsom, out_dir = out_dir, batch = batch))
    })
    names(scores) <- unique(file_batch_id)
  }
  else {
    files <- fcs_files
    fsom <- fsom_aof(fcs_files = files, phenotyping_markers = phenotyping_markers,
                     out_dir = out_dir, arcsine_transform = arcsine_transform,
                     nClus = nClus,
                     batch = NULL)

    scores <- aof_scoring(fcs_files = files,
                          phenotyping_markers = phenotyping_markers,
                          fsom = fsom, out_dir = out_dir, batch = NULL)
  }

  file_outlier_detecion(scores = scores, out_dir = out_dir,
                        sd = sd)
}


.debarcode_ind <- function(file,
                           fcs_files,
                           file_batch_id,
                           file_score,
                           out_dir,
                           min_threshold,
                           threshold,
                           to_plot,
                           barcodes_used,
                           less_than_th,
                           barcode_key) {

  print(paste0("   ", Sys.time()))
  print(paste0("   Debarcoding ", file))
  ff <- flowCore::read.FCS(file, transformation = FALSE)

  file_id <- which(file == fcs_files)
  batch_id <- file_batch_id[file_id]

  if(!is.null(barcodes_used)){
    if(is.list(barcodes_used)){
      s_key <- barcode_key[rownames(barcode_key) %in% barcodes_list[[batch_id]],]
    } else {
      s_key <- barcode_key[rownames(barcode_key) %in% barcodes_used,]
    }

  } else {
    s_key <- barcode_key
  }

  dat <- CATALYST::prepData(ff)
  dat <- CATALYST::assignPrelim(dat, bc_key = s_key)
  rownames(dat)[SummarizedExperiment::rowData(dat)$is_bc]
  # table(colData(dat)$bc_id)
  dat <- CATALYST::estCutoffs(dat)

  if (min_threshold){
    if(any(S4Vectors::metadata(dat)$sep_cutoffs < threshold)){
      warning(paste0("cutoff lower than 0.18 has been detected for ", basename(file),
                     ", cutoff will be set to 0.18"))
      fileOut <- basename(file)
    }

    id <- S4Vectors::metadata(dat)$sep_cutoffs < threshold
    S4Vectors::metadata(dat)$sep_cutoffs[id] <- threshold

  } else {
    if(any(S4Vectors::metadata(dat)$sep_cutoffs < threshold)){
      warning(paste0("cutoff lower than ", threshold, " detected for ", basename(file)))
      fileOut <- basename(file)
    }
  }

  id <- is.na(S4Vectors::metadata(dat)$sep_cutoffs)
  S4Vectors::metadata(dat)$sep_cutoffs[id] <- 1

  if (to_plot){
    p <- CATALYST::plotYields(dat, which = rownames(s_key))

    pdf(file.path(out_dir, paste(gsub(".fcs", "_yields.pdf", basename(file)))))
    for (name in names(p)){
      print(p[[name]])
    }
    dev.off()
  }

  dat <- CATALYST::applyCutoffs(dat)

  if (to_plot){
    p <- CATALYST::plotEvents(dat, n = 500)

    pdf(file.path(out_dir, paste(gsub(".fcs", "_debarcode_quality.pdf",
                                      basename(file)))))
    for (name in names(p)){
      print(p[[name]])
    }
    dev.off()
  }

  dat <- dat[, dat$bc_id !=0]
  fs <- CATALYST::sce2fcs(dat, split_by = "bc_id")

  tmp_dir <- file.path(out_dir, batch_id)
  if(!dir.exists(tmp_dir)) dir.create(tmp_dir)

  file_name <- gsub("_cleaned.fcs|.fcs", "", basename(file))

  flowCore::write.flowSet(fs, outdir = tmp_dir,
                          filename = paste0(rownames(fs@phenoData), "_", file_name,
                                            "_debarcoded.fcs"))

  filePaths <- file.path(tmp_dir, paste0(rownames(fs@phenoData), "_", file_name,
                                        "_debarcoded.fcs"))

  return(list(fileOut, filePaths))
}

#' Debarcodes files
#'
#' @description performs sample debarcoding
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param cores Number of cores to be used
#' @param file_batch_id Character vector with batch label for each fcs_file,
#' the order and the length needs to be the same as in fcs_files. If only batch
#' is processed can be prepared as e.g. file_batch_id <- rep("batch", length(files))
#' @param file_score Data frame with quality scores obtained from
#' file_quality_check.Default set to NULL.
#' @param out_dir Character, pathway to where the plots should be saved,
#' only if argument to_plot = TRUE, default is set to working directory
#' @param min_threshold Logical, if the minimal threshold for barcoding
#' should be applied.
#' @param threshold Numeric, value for the minimum threshold for debarcoding,
#' default is set to 0.18, only if min_threshold set to TRUE.
#' @param to_plot Logical, if plots for yields and debarcoding quality should
#' be plotted.
#' @param barcodes_used Character vector with the names of the barcodes that were used, eg.
#' barcode 1 is the same as A1. Or a list with the barcodes name per batch.
#' If NULL (default) all the barcodes contained in sample_key will be used,
#' regarding the batch.
#' @param less_than_th Logical, if the name of the files for which lower threshold
#' than set in parameter threshold was detected. Default is set to FALSE.
#' @param barcode_key matrix as in CATALYST::assignPrelim, the debarcoding scheme.
#' A binary matrix with sample names as row names and numeric masses as column names
#' OR a vector of numeric masses corresponding to barcode channels.
#' When the latter is supplied, 'assignPrelim' will create a scheme of the
#' appropriate format internally.
#'
#' @return save debarcoded fcs files in out_dir. If parameter to_plot set
#' to TRUE save plots for yields and debarcode_quality in out_dir. If less_than_th
#' set to TRUE save file names with lower than the value set in threshold parameter
#' are saved into file called "files_with_lower_debarcoding_threshold.RDS" in out_dir
#'
#' @examples
#'
#' # Set input directory
#'clean_dir <- file.path(dir, "Cleaned")
#'
#'# Define files for debarcoding
#' files <- list.files(clean_dir,
#'                    pattern = "_cleaned.fcs$",
#'                    full.names = TRUE)
#'
#'# Read in file scores if calculated
#'file_scores <- readRDS(list.files(path = dir,
#'                                  recursive = TRUE,
#'                                  full.names = TRUE,
#'                                  pattern = "Quality_AOF_score.RDS"))
#'
#'# Define file batch ID for each file
#'file_batch_id <- stringr::str_match(basename(files),
#'                                    "(day[0-9]*).*.fcs")[,2]
#'
#'# Read in metadata
#'md <- utils::read.csv(file.path(dir, "RawFiles", "meta_data.csv"))
#'
#'# read in barcode key
#'sample_key <- CATALYST::sample_key

#'# Extract information about barcodes used in each batch
#'barcodes_list <- list()
#'for (batch in unique(file_batch_id)){
#'  idx <- md[md[,"BATCH"] == batch, "BARCODE"]
#'  barcodes_list[[batch]] <- rownames(sample_key)[idx]
#'}
#'
#'# Debarcode files
#'debarcode_files(fcs_files = files,
#'                out_dir = NULL,
#'                file_score = file_scores,
#'                min_threshold = TRUE,
#'                barcodes_used = barcodes_list,
#'                file_batch_id = file_batch_id,
#'                less_than_th = TRUE,
#'                barcode_key = sample_key)
#'
#' @export
debarcode_files <- function(fcs_files,
                            cores = 1,
                            file_batch_id,
                            file_score = NULL,
                            out_dir = NULL,
                            min_threshold = TRUE,
                            threshold = 0.18,
                            to_plot = TRUE,
                            barcodes_used = NULL,
                            less_than_th = FALSE,
                            barcode_key = NULL){

  if(anyDuplicated(fcs_files) != 0){
    stop("names of fcs files are duplicated")
  }

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  if(!is.null(file_score)){
    if(!inherits(file_score, "data.frame")) {
      stop("file_scores is not a data frame")
    } else {
      # Select good quality files
      good_files <- file_scores$file_names[file_scores$quality == "good"]
      fcs_files_clean <- fcs_files[basename(fcs_files) %in% good_files]
      fcs_files <- fcs_files_clean
    }
  }

   if(length(file_batch_id) != length(fcs_files)){
    stop("the lenght of the file_batch_id is not equal to the lenght of fcs_files")
  }

  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Debarcoded")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  outFiles <-  BiocParallel::bplapply(fcs_files, function(file) {
    .debarcode_ind(file,
                   fcs_files = fcs_files,
                   file_batch_id = file_batch_id,
                   file_score = file_score,
                   out_dir = out_dir,
                   min_threshold = min_threshold,
                   threshold = threshold,
                   to_plot = to_plot,
                   barcodes_used = barcodes_used,
                   less_than_th = less_than_th,
                   barcode_key = barcode_key)
  },
  BPPARAM = BiocParallel::MulticoreParam(workers = cores))

  lessFiles <- lapply(outFiles, function(x) {return(x[[1]])})
  debarcodedFiles <- lapply(outFiles, function(x) {return(x[[2]])})

  if(less_than_th){
    saveRDS(unlist(lessFiles), file.path(out_dir, "files_with_lower_debarcoding_threshold.RDS"))
  }
  return(unlist(debarcodedFiles))
}


.aggregate_ind <- function(fcs_files,
                           cores = 1,
                           channels_to_keep = NULL,
                           outputFile = "aggregate.fcs",
                           maxcells = NULL,
                           write_agg_file = FALSE,
                           out_dir = getwd()) {
  nFiles <- length(fcs_files)
  flowFrame <- NULL

  if(!dir.exists(aggregate_dir))(dir.create(aggregate_dir))

  for (i in seq_len(nFiles)) {
    f <- flowCore::read.FCS(fcs_files[i])

    if(!is.null(maxcells)){
      c <- sample(seq_len(nrow(f)), min(nrow(f), maxcells))
      f <- f[c,]
    }

    m <- matrix(rep(i, nrow(f)))
    m2 <- m + stats::rnorm(length(m), 0, 0.1)
    m <- cbind(m, m2)
    colnames(m) <- c("File", "File_scattered")
    prev_agg <- length(grep("File[0-9]*$", colnames(f)))
    if (prev_agg > 0) {
      colnames(m) <- paste0(colnames(m), prev_agg + 1)
    }
    if(is.null(channels_to_keep)){
      f <- flowCore::fr_append_cols(f, m)
    } else {
      f <- flowCore::fr_append_cols(f[ , channels_to_keep], m)
    }
    if (is.null(flowFrame)) {
      flowFrame <- f
      flowFrame@description$`$FIL` <- gsub(".*/", "",
                                           outputFile)
      flowFrame@description$FILENAME <- gsub(".*/", "",
                                             outputFile)
    }
    else {
      f@exprs[, "Time"] <- f@exprs[, "Time"] + max(flowFrame@exprs[,"Time"]) + 1000
      flowCore::exprs(flowFrame) <- rbind(flowCore::exprs(flowFrame),
                                          flowCore::exprs(f))
    }
  }
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame) -
                                 1, "Rmin", sep = "")]] <- 0
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame) -
                                 1, "Rmax", sep = "")]] <- nFiles + 1
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame),
                               "Rmin", sep = "")]] <- 0
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame),
                               "Rmax", sep = "")]] <- nFiles + 1
  flowFrame@description[[paste("$P", ncol(flowFrame) - 1,
                               "B", sep = "")]] <- 32
  flowFrame@description[[paste("$P", ncol(flowFrame), "B",
                               sep = "")]] <- 32
  flowFrame@description$FIL <- gsub(".*/", "", outputFile)

  if(write_agg_file == TRUE){

    flowCore::write.FCS(x = flowFrame, filename = file.path(out_dir, outputFile), endian = "big")
  }

  return(file.path(out_dir, outputFile))

  # return(flowFrame)
}

#' Deconvolute and aggregate debarcoded files
#'
#' @description Performs aggregation of debarcoded files, assigning user defined
#' name to each file.
#'
#' @param fcs_files Character, full path to the fcs_files.
#' @param md Metadata. Must contain the following columns:
#' batch_column: defines to which batch each file belongs;
#' barcode_name: defines to which barcode each file belongs
#' fcs_new_name: a name for the fcs file that will be given after deconvolution
#' and aggregation.
#' @param cores Number of cores to be used.
#' @param channels_to_keep Character vector with channel names to be kept.
#' Default NULL.
#' @param maxcells Numeric, maximum cells to randomly aggregate from each file,
#' default is set to NULL, which means that all the cells will be aggregated.
#' @param write_agg_file Logical, if the fcs files should be saved, if TRUE
#' files will be saved in out_dir. Default set to TRUE.
#' @param out_dir Character, pathway to where the files should be saved,
#' if NULL (default) files will be saved to file.path(getwd(), Aggregated).
#'
#' @return List of the pathways to aggregated files.
#'
#' @examples
#'
#' # Set input directory
#' debarcode_dir <- file.path(dir, "Debarcoded")

#' # Define files for debarcoding
#' files <- list.files(debarcode_dir,
#'                     pattern = "_debarcoded.fcs$",
#'                     full.names = TRUE, recursive = T)
#'
#' # Define out_dir for aggregated files
#' aggregate_dir <- file.path(dir, "Aggregated")
#'
#' # Bring metadata
#' md <- utils::read.csv(file.path(dir, "RawFiles", "meta_data.csv"))
#'
#' # Assign barcodes names
#' md$barcode_name <- paste0(rownames(CATALYST::sample_key)[md$BARCODE])
#'
#' # Assign new sample names specifying patient id and its batch name
#' md$fcs_new_name <- paste0(md$ID, "_", md$STIM, "_", md$BATCH, ".fcs")
#'
#' # Aggregate and deconvolute file names
#' aggregate_files(fcs_files = files,
#'                 md,
#'                 barcode_column = "barcode_name",
#'                 batch_column = "BATCH",
#'                 cores = 1,
#'
#'                 out_dir = aggregate_dir,
#'                 write_agg_file = TRUE)
#'
#' @export
aggregate_files <- function(fcs_files,
                            md,
                            barcode_column,
                            batch_column,
                            cores = 1,
                            channels_to_keep = NULL,
                            maxcells = NULL,
                            write_agg_file = TRUE,
                            out_dir = NULL){

  # Check parameters
  if(!is(files, "character") & !is(files, "list")) {
    stop("files must be a character vector or a list")
  }

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  if (any(!is(cores, "numeric") | cores < 1)){
    stop("cores must be a positive number")
  }

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Aggregated")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  # Parallelized analysis
  aggregatedFiles <- BiocParallel::bplapply(seq_len(nrow(md)), function(i) {
    patterns <- as.character(md[i, c(barcode_column, batch_column)])

    files_to_agg <- grep(pattern = patterns[2],
                         grep(pattern = patterns[1],
                              fcs_files, value = TRUE),
                         value = TRUE)

    print(paste0("Creating ", md[[i, "fcs_new_name"]]))

    outputFile = md[[i, "fcs_new_name"]]

    .aggregate_ind(fcs_files = files_to_agg,
                   channels_to_keep = channels_to_keep,
                   outputFile = outputFile,
                   maxcells = maxcells,
                   write_agg_file = write_agg_file,
                   out_dir = out_dir)
  },
    BPPARAM = BiocParallel::MulticoreParam(workers = cores))

  return(unlist(aggregatedFiles))

}


#' Gate intact cells
#'
#' @description Performs gating of intact cells using flowDensity package.
#'
#' @param flow_frame A flowframe that contains cytometry data.
#' @param file_name Character, the file name used for saving the flow frame
#' (if save_gated_flow_frame = TRUE) and for plotting, if NULL (default)
#' the file name stored in keyword FIL will be used.
#' @param tinypeak_removal_head Numeric from 0-1, as in deGate to exclude/include
#' tiny peaks in the head of the density distribution curve for both Iridium
#' channels.
#' @param tinypeak_removal_tail The same as tinypeak_removal1 but for the tail
#' in the density distribution curve.
#' @param alpha_head Numeric, 0-1, as in deGate specify the significance of change
#' in the slope being detected at the head of the density distribution curve.
#' @param alpha_tail The same as in alpha1 but for the tail of the density
#' distribution curve.
#' @param arcsine_transform Logical, if the data should be transformed
#' with arcsine transformation and cofactor 5. If FALSE the data won't be
#' transformed, thus transformed flow frame should be used if needed.
#' Default TRUE.
#' @param save_gated_flow_frame Logical, if gated flow frame should be saved.
#' Only cells falling into intact cell region will be saved. Default set to FALSE.
#' @param suffix Character, suffix placed in the name of saved fcs file, only
#' if save_gated_flow_frame = TRUE.Defult is "_intact_gated".
#' @param out_dir Character, pathway to where the files should be saved,
#' if NULL (default) files will be saved to file.path(getwd(), Gated).
#' @param ... Additional parameters to pass to flowDensity::deGate()
#'
#' @return An untransformed flow frame with intact cells only
gate_intact_cells <- function(flow_frame,
                              file_name = NULL,
                              tinypeak_removal_head = 0.8,
                              tinypeak_removal_tail = 0.8,
                              alpha_head = 0.05,
                              alpha_tail = 0.1,
                              arcsine_transform = TRUE,
                              save_gated_flow_frame = FALSE,
                              out_dir = NULL,
                              suffix = "_intact_gated",
                              ...){

  # Check parameters
  if(!is(flow_frame, "flowFrame")) {
    stop("flow_frame must be a flow frame" )
  }

  if (is.null(file_name)){
    file_name <- ff@description$FIL
  } else {
    file_name
  }

  if(arcsine_transform == TRUE){

    ff_t <- flowCore::transform(ff,
                                flowCore::transformList(
                                  flowCore::colnames(ff)[grep("Di", flowCore::colnames(ff))],
                                  CytoNorm::cytofTransform))
  } else {
    ff_t <- ff
  }

  selection <- matrix(TRUE,
                      nrow = nrow(ff),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("intact")))

  tr <- list()
  for(m in c("Ir193Di", "Ir191Di")){

    tr[[m]] <- c(flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal_head,
                                     upper = FALSE, use.upper = TRUE,
                                     alpha = alpha_head, verbose = F, count.lim = 3, ...),
                 flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal_tail,
                                     upper = TRUE, use.upper = TRUE,
                                     alpha = alpha_tail, verbose = F, count.lim = 3, ...))
  }

  for(m in c("Ir193Di", "Ir191Di")){
    selection[ff_t@exprs[,m] < tr[[m]][1], "intact"] <- FALSE
    selection[ff_t@exprs[,m] > tr[[m]][2], "intact"] <- FALSE
  }

  percentage <- (sum(selection)/length(selection))*100
  flowDensity::plotDens(ff_t, c("Ir193Di", "Ir191Di"),
                        main = paste0(basename(file_name)," ( ", format(round(percentage, 2),
                                                                        nsmall = 2), "% )"))

  graphics::abline(h = c(tr[["Ir191Di"]]))
  graphics::abline(v = c(tr[["Ir193Di"]]))
  graphics::points(ff_t@exprs[!selection[,"intact"], c("Ir193Di", "Ir191Di")], pch = ".")

  ff <- ff[selection[,"intact"], ]

  if(save_gated_flow_frame){

    .save_flowframe(ff, out_dir, suffix, file_name)
  }

  return(ff)
}


.save_flowframe <- function(flow_frame,
                            out_dir,
                            suffix,
                            file_name){
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Gated")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  fname <- gsub(pattern = ".fcs", replacement = paste0(suffix, ".fcs"),
                x = file_name, ignore.case = TRUE)

  flowCore::write.FCS(x = flow_frame, filename = file.path(out_dir, fname))
}

#' remove_mad_outliers
#'
#' @description detects outliers in the selected channel(s) using MAD
#' (mean absolute deviation).
#'
#' @param flow_frame
#' @param channels character, channel names used for gating, default is set to
#' "Event_length"
#' @param n_mad numeric, how many MAD should be use to detect outliers
#' @param mad_f function used to compute deviation, default set to "mad"
#' @param plot logicle, if to plot the data, default TRUE
#' @param center
#' @param main character, title of the plot, default set to ""
#' @param ... other arguments to pass plotDens
#'
#' @return matrix with the selected cells
.remove_mad_outliers <- function(flow_frame,
                                channels = "Event_length",
                                n_mad = 2,
                                mad_f = mad,
                                plot = TRUE,
                                center = "center",
                                main = "",
                                ...){
  boundaries <- matrix(NA,
                       nrow = 5,
                       ncol = length(channels),
                       dimnames = list(c("median", "center", "mad", "l_lim", "u_lim"),
                                       channels))
  for (channel in channels) {
    x <- flow_frame@exprs[, channel]
    boundaries["median", channel] <- stats::median(x)
    boundaries["center", channel] <- stats::density(x)$x[which.max(stats::density(x)$y)]
    boundaries["mad", channel] <- mad_f(x,
                                        center = boundaries[center, channel] )
    boundaries["l_lim", channel] <- boundaries[center, channel] - n_mad * boundaries["mad", channel]
    boundaries["u_lim", channel] <- boundaries[center, channel] + n_mad * boundaries["mad", channel]
  }

  selection <- rep(TRUE, nrow(flow_frame))
  for (channel in channels) {
    selection <- selection & (flow_frame@exprs[, channel] > boundaries["l_lim", channel])
    selection <- selection & (flow_frame@exprs[, channel] < boundaries["u_lim", channel])
  }
  percentage <- (sum(selection)/length(selection))*100
  if (plot) {
    flowDensity::plotDens(flow_frame,
                          c(channels, "Ir191Di"),
                          main = paste0(main, " ( ", format(round(percentage, 2),
                                                            nsmall = 2), "% )"),
                          ...)
    if(length(channels) == 2) {
      graphics::points(flow_frame@exprs[!selection, channels], col = "red", pch = ".")
      graphics::abline(v = boundaries[c("l_lim", "u_lim"), channels[1]], col = "grey")
      graphics::abline(h = boundaries[c("l_lim", "u_lim"), channels[2]], col = "grey")
    } else if(length(channels) == 1) {
      graphics::points(flow_frame@exprs[!selection, c(channels, "Ir191Di")], pch = ".")
      graphics::abline(v = boundaries[c("l_lim", "u_lim"), channels[1]], col = "grey")
    }
  }

  return(selection)
}


#' Gate singlet cells
#'
#' @param flow_frame A flowframe that contains cytometry data.
#' @param file_name Character, the file name used for saving the flow frame
#' (if save_gated_flow_frame = TRUE) and for plotting, if NULL
#' the file name stored in keyword FIL will be used,
#' default is set to NULL.
#' @param channels character, channels name to be used for gating, default is
#' to Event_length
#' @param n_mad numeric, number of MADs to detect outliers
#' @param arcsine_transform Logical, if the data should be transformed
#' with arcsine transformation and cofactor 5. If FALSE the data won't be
#' transformed, thus transformed flow frame should be used if needed.
#' Default TRUE.
#' @param save_gated_flow_frame Logical, if gated flow frame should be saved.
#' Only cells falling into intact cell region will be saved. Default set to FALSE.
#' @param suffix Character, suffix placed in the name of saved fcs file, only
#' if save_gated_flow_frame = TRUE. Defult is "_singlets_gated".
#' @param out_dir Character, pathway to where the files should be saved,
#' if NULL (default) files will be saved to file.path(getwd(), Gated).
#' @param ... arguments to pass to plotDens
#'
#' @return An untransformed flow frame with singlets only
gate_singlet_cells <- function(flow_frame,
                               file_name = NULL,
                               channels = "Event_length",
                               arcsine_transform = TRUE,
                               save_gated_flow_frame = NULL,
                               suffix = "_singlets_gated",
                               n_mad = 2,
                               out_dir = NULL,
                               ...){

  # Check parameters
  if(!is(flow_frame, "flowFrame")) {
    stop("flow_frame must be a flow frame" )
  }

  if (is.null(file_name)){
    file_name <- flow_frame@description$FIL
  } else {
    file_name
  }

  if(arcsine_transform == TRUE){

    flow_frame_t <- flowCore::transform(flow_frame,
                                        flowCore::transformList(
                                          flowCore::colnames(flow_frame)[grep("Di", flowCore::colnames(flow_frame))],
                                                                CytoNorm::cytofTransform))
  } else {
    flow_frame_t <- flow_frame
  }

  selection <- matrix(TRUE,
                      nrow = nrow(flow_frame),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("singlets")))

  selection[, "singlets"] <- .remove_mad_outliers(flow_frame = flow_frame_t,
                                                 channels = channels,
                                                 main = paste("Singlets", file_name),
                                                 n_mad = n_mad,
                                                 xlim = c(0, 100), ylim = c(0, 8), ...)

  flow_frame <- flow_frame[selection[,"singlets"], ]

  if(save_gated_flow_frame){

    .save_flowframe(flow_frame, out_dir, suffix, file_name)
  }

  return(flow_frame)

}

#' Gate live cells
#'
#' @description Performs gating of live cells using flowDensity package
#'
#' @param flow_frame A flowframe that contains cytometry data.
#' @param file_name Character, the file name used for saving the flow frame
#' (if save_gated_flow_frame = TRUE) and for plotting, if NULL
#' the file name stored in keyword FIL will be used,
#' default is set to NULL.
#' @param viability_channel Character, the channel name used for viability staining
#' @param tinypeak_removal_viability, numeric from 0-1, as in deGate to exclude/include
#' tiny peaks in the tail of the density ditribution curve for both viability channel
#' @param tinypeak_removal_Iridium the same as tinypeak_removal_viablity but for
#' the head and tail of the density ditribution curve in Iridium channel
#' @param alpha_viability numeric, 0-1, as in deGate specify the significance of change
#' in the slope of viability channel
#' @param alpha_Iridium the same as in alpha_viability but for the Iridium
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5.
#' @param save_gated_flow_frame Logical, if gated flow frame should be saved.
#' Only cells falling into intact cell region will be saved. Default set to FALSE.
#' @param suffix Character, suffix placed in the name of saved fcs file, only
#' if save_gated_flow_frame = TRUE.Defult is "_intact_gated".
#' @param out_dir Character, pathway to where the files should be saved,
#' if NULL (default) files will be saved to file.path(getwd(), Gated).
#' @param ... arguments to pass to flowDensity::plotDens()
#'
#' @return An untransformed flow frame with live cells only

gate_live_cells <- function(flow_frame,
                            file_name = NULL,
                            viability_channel,
                            tinypeak_removal_viability = 0.8,
                            alpha_viability = 0.1,
                            tinypeak_removal_Iridium = 0.8,
                            alpha_Iridium = 0.05,
                            arcsine_transform = TRUE,
                            save_gated_flow_frame = FALSE,
                            suffix = "_live_gated",
                            out_dir = NULL,... ){

  # Check parameters
  if(!is(flow_frame, "flowFrame")) {
    stop("flow_frame must be a flow frame" )
  }

  ff <- flow_frame

  if (is.null(file_name)){
    file_name <- ff@description$FIL
  } else {
    file_name
  }

  if(arcsine_transform == TRUE){

    ff_t <- flowCore::transform(ff,
                                flowCore::transformList( flowCore::colnames(ff)[grep("Di",
                                                                                     flowCore::colnames(ff))],
                                              CytoNorm::cytofTransform))
  } else {
    ff_t <- ff
  }

  selection <- matrix(TRUE,
                      nrow = nrow(ff),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("live")))


  v_ch <- grep(viability_channel,  flowCore::colnames(ff), value = T)

  tr <- list()
  for(m in c("Ir191Di", v_ch)){
    if (m == v_ch) {
      upper = TRUE
      alpha = alpha_viability
      tr[[m]] <- flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal_viability,
                                     upper = upper, use.upper = TRUE,
                                     alpha = alpha, verbose = F, count.lim = 3)

    } else {
      alpha = alpha_Iridium
      tr[[m]] <- c(flowDensity::deGate(ff_t, m,
                                       tinypeak.removal = tinypeak_removal_Iridium,
                                       upper = FALSE, use.upper = TRUE,
                                       alpha = alpha,  verbose = F, count.lim = 3),
                   flowDensity::deGate(ff_t, m,
                                       tinypeak.removal = tinypeak_removal_Iridium,
                                       upper = TRUE, use.upper = TRUE,
                                       alpha = alpha, verbose = F, count.lim = 3))

    }
  }

  for(m in c(v_ch, "Ir191Di")){
    if (m == v_ch) {
      selection[ff_t@exprs[,m] > tr[[m]][1], "live"] <- FALSE
    } else {
      selection[ff_t@exprs[,m] < tr[[m]][1], "live"] <- FALSE
      selection[ff_t@exprs[,m] > tr[[m]][2], "live"] <- FALSE
    }
  }
  percentage <- (sum(selection)/length(selection))*100
  flowDensity::plotDens(ff_t, c(v_ch, "Ir191Di"),
                        main = paste0(file_name," ( ", format(round(percentage, 2),
                                                              nsmall = 2), "% )"),
                        xlim = c(0, 8), ylim = c(0, 8), ...)

  graphics::abline(h = tr[["Ir191Di"]])
  graphics::abline(v = tr[[v_ch]])

  graphics::points(ff_t@exprs[!selection[,"live"], c(v_ch, "Ir191Di")], pch = ".")

  ff <- ff[selection[,"live"], ]

  if(save_gated_flow_frame){

    .save_flowframe(flow_frame, out_dir, suffix, file_name)
  }

  return(ff)

}


















plot_gate <- function(live_out,
                      filename = "gating.png",
                      n_plots = 3) {

  n_plots <- 3

  allPlots <- lapply(live_out, function(x) {
    return(list(x[[2]][[1]], x[[2]][[2]], x[[2]][[3]]))
  })

  allPlots <- unlist(allPlots, recursive = FALSE)

  nRows <- length(live_out)

  # par(mar=c(1, 1, 1, 1))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  # par(mfcol=c(5,3),mai=c(0.5,0.5,0.5,0))
  grDevices::png(filename,
                 width = n_plots * 3,
                 height = nRows * 3)

  ggpubr::ggarrange(plotlist = allPlots, nrow = nRows, ncol = n_plots)

  dev.off()

}

.plot_batch_ind <- function(name,
                            files,
                            batch_labels,
                            batch_pattern,
                            out_dir,
                            clustering_markers,
                            arcsine_transform,
                            manual_colors,
                            cells_total,
                            transform_list,
                            n_neighbors) {

  ff_agg <- FlowSOM::AggregateFlowFrames(fileNames = files,
                                         cTotal = length(files) * cells_total,
                                         verbose = TRUE,
                                         writeMeta = FALSE,
                                         writeOutput = FALSE,
                                         outputFile = file.path(out_dir,
                                                                paste0("aggregated_for_batch_plotting.fcs")))

  if(arcsine_transform){
    ff_agg <- flowCore::transform(ff_agg,
                              flowCore::transformList(grep("Di", flowCore::colnames(ff_agg),
                                                           value = TRUE),
                                                      CytoNorm::cytofTransform))
  } else if (!is.null(transform_list)){
    ff_agg <- flowCore::transform(ff_agg, transform_list)
  } else {
    ff_agg <- ff_agg
  }

  markers <- FlowSOM::GetMarkers(ff_agg, flowCore::colnames(ff_agg))

  cl_markers <- paste(clustering_markers, collapse="|")
  cl_markers <- grep(cl_markers, markers, value = T)

  ff_agg@exprs[, names(cl_markers)] <- apply(ff_agg@exprs[, names(cl_markers)],
                                             2, function(x){
                                               q <- stats::quantile(x, 0.9999)
                                               x[x > q] <- q
                                               x
                                             })

  samp <- length(files)
  ff_samp <- ff_agg@exprs[sample(nrow(ff_agg@exprs), samp*cells_total), ]

  dimred_res <- uwot::umap(X = ff_samp[, names(cl_markers)],
                           n_neighbors = n_neighbors, scale = TRUE)

  dimred_df <- data.frame(dim1 = dimred_res[,1], dim2= dimred_res[,2],
                          ff_samp[, names(cl_markers)])

  dimred_df$file_id <- ff_samp[,"File2"]

  if(!(is.null(batch_pattern))){
    dimred_df$batch <- sapply(files[dimred_df$file_id], function(file) {
      stringr::str_match(file, batch_pattern)[,1]})
  }

  if(!(is.null(batch_labels))){
    dimred_df$batch <- batch_labels
  }


  p <- ggplot2::ggplot(dimred_df,  aes_string(x = "dim1", y = "dim2", color = "batch")) +
    geom_point(aes(color = batch), size = 3, position="jitter") +
    ggtitle(name)+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    theme(panel.background = element_rect(fill = "white", colour = "black",
                                          size = 2, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(color="black", size=26,
                                       hjust = 0.95, face = "bold"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.text.x = element_text(size = 23, color = "black"),
          strip.background = element_rect(fill = "white"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 22),
          legend.position = "bottom",
          # legend.key.size = unit(3,"point"),
          legend.key = element_blank())

  if (!is.null(manual_colors)){
    p <- p+scale_color_manual(values = manual_colors)
  }

  return(p)

}

#' Visualize batch using umap dimensional reduction
#'
#' @description Plots batch effect using UMAP and clustering markers
#'
#' @param files_before_norm Character, full path to the unnormalized fcs_files.
#' @param files_after_norm Character, full path to the normalized fcs_files.
#' @param cores Number of cores to be used
#' @param out_dir Character, pathway to where the files should be saved,
#' if NULL (default) files will be saved to file.path(getwd(), CytoNormed).
#' @param clustering_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted.
#' These markers are used for building and plotting UMAP.
#' @param arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5, default is set to TRUE.
#' @param batch_pattern Character, batch pattern to be match in the fcs file name.
#' @param manual_colors Character, vector of the colors to be used,
#' the number of colors needs to be equal to the length of batch_patter.
#' @param cells_total Number of cells to plot per each file.
#' @param transform_list Transformation list to pass to the flowCore
#' transform function, see flowCore::transformList(), if different transformation
#' than arcsine is needed. Only if arcsine_transform is FALSE. If NULL and
#' arcsine_transform = FALSE no transformation will be applied.
#' @param n_neighbors The size of local neighborhood in UMAP analysis, default
#' set to 15, as in uwot::umap().
#' It is recommended to set it to the number of files in each batch.
#'
#' @import ggplot2
#'
#' @examples
#' # Define files before normalization
#' gate_dir <- file.path(dir, "Gated")
#' files_before_norm <- list.files(gate_dir,
#'                                 pattern = ".fcs",
#'                                 full.names = T)
#'
#' # Define files after normalization
#' norm_dir <- file.path(dir, "CytoNormed")
#' files_after_norm <- list.files(norm_dir,
#'                                pattern = ".fcs",
#'                                full.names = T)
#'
#' # files needs to be in the same order, check and order if needed
#' test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
#'                  basename(files_before_norm))
#'
#' batch_labels <- stringr::str_match(basename(files_before_norm), "day[0-9]*")[,1]
#'
#' # Plot batch effect
#' set.seed(789)
#' plot_batch(files_before_norm = files_before_norm,
#'            files_after_norm = files_after_norm,
#'            batch_labels = batch_labels,
#'            cores = 1,
#'            out_dir = norm_dir,
#'            clustering_markers = c("CD", "IgD", "HLA"),
#'            manual_colors = c("darkorchid4", "darkorange", "chartreuse4"))
#'
#' @return save plots for batch effect in the out_dir

plot_batch <- function(files_before_norm,
                       files_after_norm,
                       batch_labels = NULL,
                       batch_pattern = NULL,
                       cores = 1,
                       out_dir = NULL,
                       clustering_markers = "CD|HLA|IgD|PD|BAFF|TCR",
                       arcsine_transform = TRUE,
                       manual_colors = NULL,
                       cells_total = 1000,
                       transform_list = NULL,
                       n_neighbors = length(files_before_norm)){

  if(!(length(files_after_norm) == length(files_before_norm))){
    stop("files_before and after does not have the same length")
  }

  fcs_files <- c(files_after_norm, files_before_norm)

  if (!all(file.exists(fcs_files))){
    stop("incorrect file path, the fcs file does not exist")
  }

  test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
                   basename(files_before_norm))

  files_list <- list("files_before_norm" = files_before_norm,
                     "files_after_norm" = files_after_norm)


  if(is.null(batch_labels) & is.null(batch_pattern)){
    stop("define batch_labels or batch_pattern")
  } else if (!(is.null(batch_labels)) & !(is.null(batch_pattern))){
    stop("both batch_labels and batch_pattern are defined, desellect one option by  setting to NULL")
  }

  if(!is.null(batch_labels)){
    if(length(files_after_norm) != length(batch_labels)){
      stop("The lenght of batch labels is not equal to the lenght of files")
    }
  }

    if(is.null(out_dir)){
      out_dir <- file.path(getwd(), "CytoNormed")
    }
    if(!dir.exists(out_dir)){dir.create(out_dir)}


  # Parallelized analysis
  plots <- BiocParallel::bplapply(names(files_list), function(x) {
    .plot_batch_ind(name = x,
                    files = files_list[[x]],
                    batch_labels = batch_labels,
                    batch_pattern = batch_pattern,
                    out_dir = out_dir,
                    clustering_markers = clustering_markers,
                    arcsine_transform = arcsine_transform,
                    manual_colors = manual_colors,
                    cells_total = cells_total,
                    transform_list = transform_list,
                    n_neighbors = n_neighbors)},
    BPPARAM = BiocParallel::MulticoreParam(workers = cores))



  png(file.path(norm_dir, "batch_effect.png"),
      width = length(plots)*1500,
      height = 1500, res = 300)
  gridExtra::grid.arrange(grobs = plots, ncol = 2)
  dev.off()
}

test_match_order <- function(x,y) {

  if (isTRUE(all.equal(x,y))) print('Files are ordered')

  if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y))))
  warning('Perfect match but wrong order. Please order the files')

  if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y))))
    warning('No match, please make sure that files are in the same order')
}


#' Prepares data for plotting cell frequency and MSI
#'
#' @description Performs dimensional reduction and constructs
#' data frame for plotting cell frequencies and MSI
#' per clusters and metaclusters obtained from extract_pctgs_msi_per_flowsom function.
#'
#' @param frequency_msi_list List containing matrices with cell frequency and msi
#' obtained in step extract_pctgs_msi_per_flowsom.
#' @param matrix_type The name of the matrix to be plotted.
#' @param seed Numeric set to obtain reproducible results, default NULL.
#' @param n_neighbours The size of local neighborhood in UMAP analysis, default
#' set to 15, as in uwot::umap().
#' It is recommended to set it to the number of files in each batch.
#'
#' @return data frame for plotting
#'
#' @examples
#' # Define files before normalization
#' gate_dir <- file.path(dir, "Gated")
#' files_before_norm <- list.files(gate_dir,
#'                                 pattern = ".fcs",
#'                                 full.names = T)
#'
#' # Define files after normalization
#' norm_dir <- file.path(dir, "CytoNormed")
#' files_after_norm <- list.files(norm_dir,
#'                                pattern = ".fcs",
#'                                full.names = T)
#'
#' # files needs to be in the same order, check and order if needed
#' test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
#'                  basename(files_before_norm))
#'
#' batch_labels <- stringr::str_match(basename(files_before_norm), "day[0-9]*")[,1]
#'
#' mx <- extract_pctgs_msi_per_flowsom(files_after_norm = files_after_norm,
#'                                     files_before_norm = files_before_norm,
#'                                     nCells = 50000,
#'                                     phenotyping_markers = c("CD", "HLA", "IgD"),
#'                                     functional_markers = c("MIP", "MCP", "IL",
#'                                                            "IFNa", "TNF", "TGF",
#'                                                            "Gr"),
#'                                     xdim = 10,
#'                                     ydim = 10,
#'                                     n_metaclusters = 35,
#'                                     out_dir = norm_dir,
#'                                     arcsine_transform = TRUE,
#'                                     save_matrix = TRUE,
#'                                     seed = 343)
#' # create the list to store the plots
#'  plots <- list()
#'  for (name in names(mx[[1]])){
#'  df_plot <- prepare_data_for_plotting(frequency_msi_list = mx,
#'                                        matrix_type = name,
#'                                       n_neighbours = 11, seed = 1)
#'
#'
#'  batch <- stringr::str_match(rownames(df_plot), "day[0-9]*")[,1]
#'  samples_id <- ifelse(grepl("p1", rownames(df_plot)),"p1",
#'                       ifelse(grepl("p2", rownames(df_plot)), "p2", "ref"))
#'  stimulation <- stringr::str_match(rownames(df_plot), "UNS|RSQ|IMQ|LPS|CPG")[,1]
#'
#'  plots[[name]] <- plot_batch_using_freq_msi(df_plot = df_plot, fill = batch,
#'                                             shape = samples_id, color = batch,
#'                                             split_by_normalization = TRUE, title = name)
#'
#'}
#'
#' gg_a <- grid_arrange_common_legend(plot_lists = plots, nrow = 2, ncol = 2,
#'                                   position = "right")
#'
#'ggplot2::ggsave(filename ="batch_effect_frequency_MSI.png",
#'                device = "png",
#'               path = norm_dir,
#'                plot = gg_a,
#'                units = "cm",
#'                width = 22,
#'                height = 14, dpi = 300)
#'
prepare_data_for_plotting <- function(frequency_msi_list,
                                      matrix_type,
                                      n_neighbours = 15,
                                      seed = NULL){
  print(matrix_type)

  # set the number of closest neighborhoods
  n <- n_neighbours

  # process files before normalization
  df_b <- frequency_msi_list[["before"]][[matrix_type]]

  # select the columns for which MSI is higher than 0.2 SD
  if(grepl("mfi", matrix_type)){
    id_cols <-  which(apply(df_b, 2, mean) > 1)
    df_b <- df_b[,id_cols]
  }

  # build UMAP for files before the normalization
  if(!is.null(seed)){
    set.seed(seed)
  }
  df_b_umap <- data.frame(uwot::umap(df_b, n_neighbors = n, scale = T))

  # process files after normalization
  df_a <- frequency_msi_list[["after"]][[matrix_type]]

  # select the columns for which MSI is higher than 0.2 SD
  if(grepl("mfi", matrix_type)){
    id_cols <-  which(apply(df_a, 2, mean) > 1)
    df_a <- df_a[,id_cols]
  }

  # build UMAP for files after the normalization
  if(!is.null(seed)){
    set.seed(seed)
  }
  df_a_umap <- data.frame(uwot::umap(df_a, n_neighbors = n, scale = T))

  # extract rownames to use the for ggplot annotation
  rnmes <- c(rownames(df_b), rownames(df_a))

  #join two UMAP data frames
  dr <- data.frame(rbind(df_b_umap, df_a_umap), check.names = F)
  colnames(dr) <- c("dim1", "dim2")

  dr$normalization <- c(rep("Raw", length(rownames(df_b_umap))),
                        rep("Normalized", length(rownames(df_a_umap))))

  return(dr)
}


#' Plot the distribution of the features across batches in two dimensional space
#'
#' @description Uses ggplot draw the distribution of the samples across batches
#'
#' @param df_plot Data frame that contains dim1 and dim2 obtained upon dimensional
#' reduction. This data frame is generated by prepare_data_for_plotting.
#' @param fill As in ggplot2, defines by which variable the points are colored.
#' @param shape As in ggplot2, defines by which variable the shapes are plotted.
#' @param color As in ggplot2, defines by which variable the borderd of the
#' points are colored.
#' @param split_by_batch Logical, if TRUE (defult) the plots are split
#' (facet_wrap) by normalization.
#' @param fill_legend_name Character, specifies the name of the legend for fill.
#' @param color_legend_name Character, specifies the name of the legend for fill.
#' @param shape_legend_name Character, specifies the name of the legend for fill.
#' @param title Character, the title of the plot.
#'
#' @return ggplot
#' @export
#'
#' @examples
#' # Define files before normalization
#' gate_dir <- file.path(dir, "Gated")
#' files_before_norm <- list.files(gate_dir,
#'                                 pattern = ".fcs",
#'                                 full.names = T)
#'
#' # Define files after normalization
#' norm_dir <- file.path(dir, "CytoNormed")
#' files_after_norm <- list.files(norm_dir,
#'                                pattern = ".fcs",
#'                                full.names = T)
#'
#' # files needs to be in the same order, check and order if needed
#' test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
#'                  basename(files_before_norm))
#'
#' batch_labels <- stringr::str_match(basename(files_before_norm), "day[0-9]*")[,1]
#'
#' mx <- extract_pctgs_msi_per_flowsom(files_after_norm = files_after_norm,
#'                                     files_before_norm = files_before_norm,
#'                                     nCells = 50000,
#'                                     phenotyping_markers = c("CD", "HLA", "IgD"),
#'                                     functional_markers = c("MIP", "MCP", "IL",
#'                                                            "IFNa", "TNF", "TGF",
#'                                                            "Gr"),
#'                                     xdim = 10,
#'                                     ydim = 10,
#'                                     n_metaclusters = 35,
#'                                     out_dir = norm_dir,
#'                                     arcsine_transform = TRUE,
#'                                     save_matrix = TRUE,
#'                                     seed = 343)
#' # create the list to store the plots
#'  plots <- list()
#'  for (name in names(mx[[1]])){
#'  df_plot <- prepare_data_for_plotting(frequency_msi_list = mx,
#'                                        matrix_type = name,
#'                                       n_neighbours = 11, seed = 1)
#'
#'
#'  batch <- stringr::str_match(rownames(df_plot), "day[0-9]*")[,1]
#'  samples_id <- ifelse(grepl("p1", rownames(df_plot)),"p1",
#'                       ifelse(grepl("p2", rownames(df_plot)), "p2", "ref"))
#'  stimulation <- stringr::str_match(rownames(df_plot), "UNS|RSQ|IMQ|LPS|CPG")[,1]
#'
#'  plots[[name]] <- plot_batch_using_freq_msi(df_plot = df_plot, fill = batch,
#'                                             shape = samples_id, color = batch,
#'                                             split_by_normalization = TRUE, title = name)
#'
#'}
#'
#' gg_a <- grid_arrange_common_legend(plot_lists = plots, nrow = 2, ncol = 2,
#'                                   position = "right")
#'
#'ggplot2::ggsave(filename ="batch_effect_frequency_MSI.png",
#'                device = "png",
#'               path = norm_dir,
#'                plot = gg_a,
#'                units = "cm",
#'                width = 22,
#'                height = 14, dpi = 300)
#'
plot_batch_using_freq_msi <- function(df_plot,
                                      fill = NULL,
                                      shape = NULL,
                                      color = NULL,
                                      split_by_normalization = TRUE,
                                      fill_legend_name = NULL,
                                      color_legend_name = NULL,
                                      shape_legend_name = NULL,
                                      title = NULL){


  p <- ggplot(df_plot, aes(x = dim1, y = dim2))+
    geom_point(data=df_plot, aes(x=dim1, y=dim2, fill = fill,
                                        shape = shape, color = color),
               size = 3)+
    ggtitle(title)

  if(split_by_normalization){
    df_plot$normalization <- factor(df_plot$normalization,
                                    levels = c("Raw", "Normalized"))
    p <- p + facet_wrap(~normalization)
  }

  p <- p + theme(panel.background = element_rect(fill = "white", colour = "black",
                                                 size = 1, linetype = "solid"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 # axis.text = element_blank(),
                 # axis.ticks = element_blank(),
                 # axis.title.y = element_blank(),
                 # axis.title.x = element_blank(),
                 legend.position = "right",
                 legend.key=element_blank(),
                 # title = element_text(size = 10),
                 #strip.text = element_blank(),
                 strip.background = element_rect(fill = "white", colour = "black"))


  p <- p + labs(fill = fill_legend_name, color = color_legend_name, shape = shape_legend_name)
  return(p)
}


#' UMAP
#'
#' @description Builds UMAP on aggregated flow frame
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param clustering_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted
#' @param functional_markers Character vector, marker names for functional markers
#' e.g cytokines, phosphorylated proteins etc. The functional markers will be plotted
#' in the file "Marker distribution across aliquots and batches"
#' @param out_dir Character, pathway to where the plots should be saved,
#' default is set to working directory.
#' @param batch_pattern Character, bach pattern to be match in the fcs file name
#' @param arcsine_transform arcsine_transform Logical, if the data should be transformed with
#' arcsine transformation and cofactor 5, default is set to TRUE
#' @param cells_total numeric, number of cells taken from each file to buil UMAP
#' @param seed numeric set to obtain reproducible results, default is set to 1
#'
#' @return data frame with UMAP coordinates

UMAP <- function(fcs_files,
                 clustering_markers = c("CD", "HLA", "IgD"),
                 functional_markers = c("IL", "TNF", "TGF", "Gr", "IF"),
                 out_dir = getwd(),
                 batch_pattern = "day[0-9]*",
                 arcsine_transform = TRUE,
                 cells_total = 1000){

  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  ff_agg <- FlowSOM::AggregateFlowFrames(fileNames = fcs_files,
                                         cTotal = length(fcs_files) * cells_total,
                                         verbose = TRUE,
                                         writeMeta = TRUE,
                                         writeOutput = TRUE,
                                         outputFile = file.path(out_dir, paste0("aggregated_for_UMAP_analysis.fcs")))

  if (arcsine_transform == TRUE){
    ff_agg <- flowCore::transform(ff_agg,
                                  transformList(grep("Di", colnames(ff_agg),
                                                     value = TRUE),
                                                CytoNorm::cytofTransform))
  }

  markers <- FlowSOM::GetMarkers(ff_agg, colnames(ff_agg))

  cl_markers <- paste(clustering_markers, collapse="|")
  cl_markers <- grep(cl_markers, markers, value = T)

  ff_agg@exprs[, names(cl_markers)] <- apply(ff_agg@exprs[, names(cl_markers)], 2, function(x){
    q <- quantile(x, 0.9999)
    x[x > q] <- q
    x
  })

  dimred_res <- uwot::umap(X = ff_agg@exprs[, names(cl_markers)],
                           n_neighbors = 15, scale = TRUE)

  if (!is.null(functional_markers)){
    f_markers <- paste(functional_markers, collapse="|")
    f_markers <- grep(f_markers, markers, value = T)
    all_markers <- c(cl_markers, f_markers)
  } else {
    all_markers <- cl_markers
  }

  df <- data.frame(dim1 = dimred_res[,1], dim2= dimred_res[,2],
                   ff_agg@exprs[, names(all_markers)])

  colnames(df)[3:ncol(df)] <- gsub("_|-", "",
                                   sub("^[0-9]*[A-Za-z]*|^[0-9]*[A-Za-z]*_","",
                                       all_markers))

  df$file_id <- ff_agg@exprs[,"File2"]
  df$batch2 <- sapply(df$file_id, function(id) {
    stringr::str_match(fcs_files[id], batch_pattern)[,1]
  })
  df$sample_name2 <-  sapply(df$file_id, function(id) {
    gsub(".fcs|_gated.fcs|_CC_gated.fcs","", basename(fcs_files[id]))  })

  return(df)

}

#' manual_labels
#'
#' @description Get the vector of manual labels for each cell
#'
#' @param manual_matrix matrix with TRUE and FALSE values for each cell and
#' each population name
#' @param cell_types character, the cells of interest
#'
#' @return vector of manual labels for each cells

manual_labels <- function (manual_matrix, cell_types)
{
  if (is.list(manual_matrix)) {
    manual_matrix <- do.call(rbind, manual_matrix)
  }
  manual <- rep("Unknown", nrow(manual_matrix))
  for (cellType in cell_types) {
    manual[manual_matrix[, cellType]] <- cellType
  }
  manual <- factor(manual, levels = c("Unknown", cell_types))
  return(manual)
}

#'@description Calculates breaks for flow frame splitting
.make_breaks <- function(event_per_flow_frame, total_events){
  breaks <- .split_flowFrames(seq_len(total_events),
                             event_per_flow_frame)

  names(breaks) <- seq_along(breaks)

  return(list("breaks"=breaks, "events_per_flowframe"=event_per_flow_frame))
}

#'@description Calculates beginning and end of each flow frame.
#'Code from Annelies Emmaneel (2021). PeacoQC:
#'Peak-based selection of high quality cytometry data. R package
#'version 1.4.0.
#'
.split_flowFrames <- function(vec, seg.length) {
  starts=seq(1, length(vec), by=seg.length)
  ends  = starts + seg.length - 1
  ends[ends > length(vec)]=length(vec)

  lapply(seq_along(starts), function(i) vec[starts[i]:ends[i]])
}

#' Split big flow frames into smaller fcs files
#'
#' @description Split the big flow frames into smaller ones.
#'
#' @param flow_frame Flow frame containing cytometry data.
#' @param event_per_flow_frame Numeric, the number of events to be split to
#' small flow frames, default is set to 500000.
#' @param min_cell_per_fcs Numeric, minimal number of cells in flow frame
#' to save fcs file, default 20000.
#' @param out_dir Character, pathway to where the files should be saved,
#' default is set to file.path(getwd(), Splitted).
#'
#' @return Save splitted fcs files in out_dir.
split_big_flowFrames <- function(flow_frame,
                                 event_per_flow_frame = 500000,
                                 out_dir = NULL,
                                 min_cell_per_fcs = 20000){

  total_events <- nrow(flow_frame)
  res_breaks <- .make_breaks(event_per_flow_frame, total_events)

  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "Splitted")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}

  for (i in as.numeric((names(res_breaks$breaks)))) {
    id <- res_breaks$breaks[[i]]

    if (i < 10) {
      num <- paste0("0", i)
    }

    if(nrow(flow_frame[id, ]) > min_cell_per_fcs){
      flowCore::write.FCS(flow_frame[id, ],
                          file.path(out_dir,
                                    gsub(".fcs|.FCS",
                                         paste0("_", num, ".fcs"), flow_frame@description$ORIGINALGUID)))
    }

  }
}


#' Extracts percentages and MSI for cell populations
#'
#' @description Performs FlowSOM clustering and extracts cluster and metacluster
#' frequency and MSI. It is imputing 0 values when NAs are detected in MSI for
#' clusters and metaclusters.
#'
#' @param file_list List, pathway to the files before and after normalization
#' @param nCells Numeric, number of cells to be cluster per each file,
#' default is set to 50 000.
#' @param phenotyping_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted.
#' @param functional_markers Character vector, marker names to be used for
#' functional markers, can be full marker name
#' e.g. "IL-6" or "IL" if all IL-markers needs to be plotted.
#' @param xdim Numeric, parameter to pass to FlowSOM, width of the SOM grid,
#' default is set to 10.
#' @param ydim Numeric, parameter to pass to FlowSOM, geight of the SOM grid,
#' default is set to 10.
#' @param n_metaclusters Numeric, exact number of clusters for metaclustering
#' in FlowSOM, default is set to 35.
#' @param out_dir Character, pathway to where the FlowSOM clustering plot should
#' be saved, default is set to file.path(getwd(), "CytoNormed").
#' @param seed Numeric, set to obtain reproducible results, default is set to NULL.
#' @param arcsine_transform arcsine_transform Logical, if the data should
#' be transformed with arcsine transformation and cofactor 5, default is set to TRUE
#' @param save_matrix Logical, if the results should be saved, if TRUE (default)
#' list of matrices will be saved in out_dir.
#' @param transform_list Transformation list to pass to the flowCore
#' transform function, see flowCore::transformList, if different transformation
#' than arcsine is needed. Only if arcsine_transform is FALSE. If NULL and
#' arcsine_transform = FALSE no transformation will be applied.Default set to NULL.
#' @param save_flowsom_result Logical, if FlowSOM and FlowSOM plots should be
#' saved. If TRUE (default) files will be saved in out_dir.
#'
#' @import ggplot2
#'
#' @examples
#' #' # Define files before normalization
#' gate_dir <- file.path(dir, "Gated")
#' files_before_norm <- list.files(gate_dir,
#'                                 pattern = ".fcs",
#'                                 full.names = T)
#'
#' # Define files after normalization
#' norm_dir <- file.path(dir, "CytoNormed")
#' files_after_norm <- list.files(norm_dir,
#'                                pattern = ".fcs",
#'                                full.names = T)
#'
#' # files needs to be in the same order, check and order if needed
#' test_match_order(x = basename(gsub("Norm_","",files_after_norm)),
#'                  basename(files_before_norm))
#'
#' batch_labels <- stringr::str_match(basename(files_before_norm), "day[0-9]*")[,1]
#'
#' mx <- extract_pctgs_msi_per_flowsom(files_after_norm = files_after_norm,
#'                                     files_before_norm = files_before_norm,
#'                                     nCells = 50000,
#'                                     phenotyping_markers = c("CD", "HLA", "IgD"),
#'                                     functional_markers = c("MIP", "MCP", "IL",
#'                                                            "IFNa", "TNF", "TGF",
#'                                                            "Gr"),
#'                                     xdim = 10,
#'                                     ydim = 10,
#'                                     n_metaclusters = 35,
#'                                     out_dir = norm_dir,
#'                                     arcsine_transform = TRUE,
#'                                     save_matrix = TRUE,
#'                                     seed = 343)
#'
#' @return list of matrices that contain calculation for
#' cl_pctgs (cluster percentages), mcl_pctgs (metaclusters percentages),
#' cl_msi (cluster MSIs for selected markers), mcl_msi (metaclusters MSI
#' for selected markers). If save_matrix = TRUE, saves this matrices in out_dir.
#' FlowSOM objects for normalized and unnormalized data,
#' if save_flowsom_result set to TRUE.
#'
#' @export
extract_pctgs_msi_per_flowsom <- function(files_before_norm,
                                          files_after_norm,
                                          nCells = 50000,
                                          phenotyping_markers = c("CD", "HLA", "IgD"),
                                          functional_markers = NULL,
                                          xdim = 10,
                                          ydim = 10,
                                          n_metaclusters = 35,
                                          out_dir = NULL,
                                          seed = NULL,
                                          arcsine_transform = TRUE,
                                          save_matrix = TRUE,
                                          transform_list = NULL,
                                          save_flowsom_result = TRUE) {

  if (!all(file.exists(c(files_after_norm, files_before_norm)))){
    stop("incorrect file path, the fcs file does not exist")
  }

  file_list <- list("before" = files_before_norm,
                    "after" = files_after_norm)

    if(is.null(out_dir)){
      out_dir <- file.path(getwd(), "CytoNormed")
    }
    if(!dir.exists(out_dir)){dir.create(out_dir)}

  res <- list()
  for (f in names(file_list)){

    nCells <- length(file_list[[f]]) * 50000
    print(paste("aggregating files for", f, "normalization"))
    ff_agg <- FlowSOM::AggregateFlowFrames(fileNames = file_list[[f]],
                                  cTotal = nCells,
                                  writeOutput = F,
                                  outputFile = file.path(out_dir, paste0(f, "_flowsom_agg.fcs")))


    if(arcsine_transform){
      ff_aggt <- flowCore::transform(ff_agg,
                                flowCore::transformList(grep("Di", flowCore::colnames(ff_agg),
                                                             value = TRUE),
                                                        CytoNorm::cytofTransform))
    } else if (!is.null(transform_list)){
      ff_aggt <- flowCore::transform(ff_agg, transform_list)
    } else {
      ff_aggt <- ff_agg
    }

    markers <- FlowSOM::GetMarkers(ff_agg, flowCore::colnames(ff_agg))
    phenotyping_channels <- grep(paste(phenotyping_markers,
                                       collapse = ("|")), markers, value = TRUE)
    functional_channels <- grep(paste(functional_markers,
                                      collapse = ("|")), markers, value = TRUE)

    # Define parameters for FlowSOM analysis
    xdim <- xdim
    ydim <- ydim
    nClus <- n_metaclusters
    s <- seed

    print(paste("building FlowSOM for", f, "normalization"))
    fsom <- FlowSOM::FlowSOM(ff_aggt,
                             colsToUse = names(phenotyping_channels),
                             scale = FALSE,
                             nClus = nClus,
                             seed = s,
                             xdim = xdim,
                             ydim = ydim)

    if(save_flowsom_result){
      fsomPlot <- FlowSOM::PlotStars(fsom, backgroundValues = fsom$metaclustering)
      fsomTsne <- FlowSOM::PlotDimRed(fsom = fsom, cTotal = 5000, seed = s)

      figure <- suppressWarnings(ggpubr::ggarrange(fsomPlot, fsomTsne,
                                                   # labels = c("FlowSOM clustering", "tsne"),
                                                   ncol = 2, nrow = 1))

      ggplot2::ggsave(filename = paste0(f, "_FlowSOM.pdf"), plot = figure, device = "pdf", path = out_dir,
                      width =24, height = 10)

      saveRDS(object = fsom, file = file.path(out_dir, paste0(f, "_flowsom.RDS")))
    }


    # Define matrices for frequency (pctgs) calculation and MSI (msi). These calculation is performed
    # for clusters (cl) and metaclusters (mcl)
    cl_pctgs <- matrix(data = NA, nrow = length(file_list[[f]]),
                       ncol = xdim * ydim,
                       dimnames = list(basename(file_list[[f]]), 1:(xdim*ydim)))

    mcl_pctgs <- matrix(data = NA, nrow = length(file_list[[f]]),
                        ncol = nClus,
                        dimnames = list(basename(file_list[[f]]), 1:nClus))
    mfi_cl_names <- apply(expand.grid(paste0("Cl", seq_len(fsom$map$nNodes)),
                                      FlowSOM::GetMarkers(ff_agg,
                                                          unique(c(phenotyping_channels,functional_channels)))),
                          1, paste, collapse = "_")
    mfi_mc_names <- apply(expand.grid(paste0("MC", 1:nClus),
                                      FlowSOM::GetMarkers(ff_agg,
                                                          unique(c(phenotyping_channels,functional_channels)))),
                          1, paste, collapse = "_")
    cl_msi <- matrix(NA,
                     nrow = length(file_list[[f]]),
                     ncol = fsom$map$nNodes * length(unique(names(c(phenotyping_channels,
                                                             functional_channels)))),
                     dimnames = list(basename(file_list[[f]]), mfi_cl_names))
    mcl_msi <- matrix(NA,
                      nrow = length(file_list[[f]]),
                      ncol =  length(mfi_mc_names),
                      dimnames = list(basename(file_list[[f]]), mfi_mc_names))

    print(paste("calculating frequency and msi for:", f, "normalization"))

    for (i in unique(fsom$data[,"File2"])){

      file <- basename(file_list[[f]][i])

      id <- which(fsom$data[,"File2"] == i)
      fsom_subset <- FlowSOM::FlowSOMSubset(fsom = fsom, ids = id)

      cl_counts <- rep(0, xdim * ydim)
      counts_tmp <- table(FlowSOM::GetClusters(fsom_subset))
      cl_counts[as.numeric(names(counts_tmp))] <- counts_tmp

      cl_pctgs[file,] <- (cl_counts/sum(cl_counts, na.rm = T))*100

      mcl_counts <- tapply(cl_counts, fsom$metaclustering, sum)
      mcl_pctgs[file,] <- tapply(cl_pctgs[file,], fsom$metaclustering, sum)

      cluster_mfis <- FlowSOM::GetClusterMFIs(fsom_subset)
      cl_msi[file,] <- as.numeric(cluster_mfis[,unique(names(c(phenotyping_channels,functional_channels)))])
      mcluster_mfis <- as.matrix(FlowSOM::GetMetaclusterMFIs(list(FlowSOM = fsom_subset,
                                                     metaclustering = fsom$metaclustering)))
      mcl_msi[file,] <- as.numeric(mcluster_mfis[,unique(names(c(phenotyping_channels,functional_channels)))])

    }

    # impute 0 values for NAs
    mfi_cl_imp <- apply(cl_msi, 2,
                        function(x){
                          missing <- which(is.na(x))
                          x[missing] <- 0
                          x
                        })

    mfi_mcl_imp <- apply(mcl_msi, 2,
                         function(x){
                           missing <- which(is.na(x))
                           x[missing] <- 0
                           x
                         })

    # store the matrices in the list for convenient plotting
    all_mx <- list("Cluster_frequencies" = cl_pctgs,
                   "Metacluster_frequencies" = mcl_pctgs,
                   "Cluster_MSIs" = mfi_cl_imp,
                   "Metacluster_MSIs" = mfi_mcl_imp)

    res[[f]] <- all_mx
  }

  if(save_matrix){
    saveRDS(object = res, file = file.path(out_dir, "cell_frequency_and_msi_list_using_FlowSOM.RDS"))
  }

  return(res)

}


##### flowAI internal functions ### taken from the package


#' #' @references this code uses internal functions from flowAI package
#' #' Monaco, G., Chen, H., Poidinger, M., Chen, J., de Magalhães, J.P.,
#' #' and Larbi, A. (2016). flowAI: automatic and interactive anomaly discerning
#' #' tools for flow cytometry data. Bioinformatics 32, 2473–2480.
.flow_rate_check_adapted <- function (x, FlowRateData, alpha = alpha,
                                      use_decomp = use_decomp)
{
  fr_frequences <- FlowRateData$frequencies
  fr_cellBinID <- FlowRateData$cellBinID
  second_fraction <- FlowRateData$info["second_fraction"]
  if (length(unique(fr_frequences[, 2])) == 1) {
    fr_autoqc <- NULL
  }
  else {
    fr_autoqc <- .anomaly_detection_addapted(fr_frequences[, "tbCounts"],
                                             alpha = alpha,
                                             use_decomp = use_decomp)
  }
  if (is.null(fr_autoqc) || is.null(fr_autoqc$anoms)) {
    badPerc <- 0
    newx <- x
    goodCellIDs <- fr_cellBinID$cellID
    badCellIDs <- NULL
  }
  else {
    goodCellIDs <- fr_cellBinID$cellID[!(fr_cellBinID$binID %in%
                                           fr_autoqc$anoms$index)]
    badCellIDs <- setdiff(fr_cellBinID$cellID, goodCellIDs)
    badPerc <- round(1 - (length(goodCellIDs)/nrow(fr_cellBinID)),
                     4)
    params <- flowCore::parameters(x)
    keyval <- flowCore::keyword(x)
    sub_exprs <- flowCore::exprs(x)
    sub_exprs <- sub_exprs[goodCellIDs, ]
    newx <- flowCore::flowFrame(exprs = sub_exprs, parameters = params,
                                description = keyval)
  }
  cat(paste0(100 * badPerc, "% of anomalous cells detected in the flow rate check. \n"))
  return(list(anoms = fr_autoqc$anoms, frequencies = fr_frequences,
              FRnewFCS = newx, goodCellIDs = goodCellIDs, badCellIDs = badCellIDs,
              res_fr_QC = data.frame(second_fraction = second_fraction,
                                     num_obs = fr_autoqc$num_obs, badPerc = badPerc)))
}


.anomaly_detection_addapted <- function (x, max_anoms = 0.49, direction = "both", alpha = 0.01,
                                         use_decomp = TRUE, period = 1, verbose = FALSE)
{
  if (is.vector(x) && is.numeric(x)) {
    x <- ts(x, frequency = period)
  }
  else if (is.ts(x)) {
  }
  else {
    stop("data must be a time series object or a vector that holds numeric values.")
  }
  if (length(rle(is.na(c(NA, x, NA)))$values) > 3) {
    stop("Data contains non-leading NAs. We suggest replacing NAs with interpolated values (see na.approx in Zoo package).")
  }
  else {
    x <- na.omit(x)
  }
  if (max_anoms > 0.49) {
    stop(paste("max_anoms must be less than 50% of the data points (max_anoms =",
               round(max_anoms * length(x), 0), " data_points =",
               length(x), ")."))
  }
  if (!direction %in% c("pos", "neg", "both")) {
    stop("direction options are: pos | neg | both.")
  }
  if (!(0.01 <= alpha || alpha <= 0.1)) {
    print("Warning: alpha is the statistical significance level, and is usually between 0.01 and 0.1")
  }
  if (is.null(period)) {
    stop("Period must be set to the number of data points in a single period")
  }
  if (use_decomp) {
    x_cf <- .cffilter_adapted(x)
    med_t <- trunc(median(x_cf$trend))
    sign_n <- sign(x_cf$trend - med_t)
    sign_n[which(sign_n == 0)] <- 1
    x_2 <- as.vector(trunc(abs(x - med_t) + abs(x_cf$cycle)) *
                       sign_n)
    trend <- x_cf$trend
  }
  else {
    x_2 <- as.vector(x - median(x))
    trend <- x
  }
  anomaly_direction = switch(direction, pos = data.frame(one_tail = TRUE,
                                                         upper_tail = TRUE), neg = data.frame(one_tail = TRUE,
                                                                                              upper_tail = FALSE), both = data.frame(one_tail = FALSE,
                                                                                                                                     upper_tail = TRUE))
  n <- length(x_2)
  data_det <- data.frame(index = 1:length(x), values = x_2,
                         or_values = trend)
  max_outliers <- trunc(n * max_anoms)
  func_ma <- match.fun(median)
  func_sigma <- match.fun(IQR)
  R_idx <- 1L:max_outliers
  num_anoms <- 0L
  one_tail <- anomaly_direction$one_tail
  upper_tail <- anomaly_direction$upper_tail
  for (i in 1L:max_outliers) {
    if (verbose)
      message(paste(i, "/", max_outliers, "completed"))
    if (one_tail) {
      if (upper_tail) {
        ares <- data_det[[2L]] - func_ma(data_det[[2L]])
      }
      else {
        ares <- func_ma(data_det[[2L]]) - data_det[[2L]]
      }
    }
    else {
      ares = abs(data_det[[2L]] - func_ma(data_det[[2L]]))
    }
    data_sigma <- func_sigma(ares)
    if (data_sigma == 0)
      break
    ares <- ares/data_sigma
    R <- max(ares)
    temp_max_idx <- which(ares == R)[1L]
    R_idx[i] <- data_det[[1L]][temp_max_idx]
    data_det <- data_det[-which(data_det[[1L]] == R_idx[i]),
    ]
    if (one_tail) {
      p <- 1 - alpha/(n - i + 1)
    }
    else {
      p <- 1 - alpha/(2 * (n - i + 1))
    }
    t <- qt(p, (n - i - 1L))
    lam <- t * (n - i)/sqrt((n - i - 1 + t^2) * (n - i +
                                                   1))
    if (R > lam)
      num_anoms <- i
  }
  if (num_anoms > 0) {
    R_idx <- R_idx[1L:num_anoms]
    all_data <- data.frame(index = 1:length(x), anoms = x)
    anoms_data <- subset(all_data, (all_data[[1]] %in% R_idx))
  }
  else {
    anoms_data <- NULL
  }
  return(list(anoms = anoms_data, num_obs = n))
}

.cffilter_adapted <- function (x, pl = NULL, pu = NULL, root = FALSE, drift = FALSE,
                               type = c("asymmetric", "symmetric", "fixed", "baxter-king",
                                        "trigonometric"), nfix = NULL, theta = 1)
{
  type = match.arg(type)
  if (is.null(root))
    root <- FALSE
  if (is.null(drift))
    drift <- FALSE
  if (is.null(theta))
    theta <- 1
  if (is.null(type))
    type <- "asymmetric"
  if (is.ts(x))
    freq = frequency(x)
  else freq = 1
  if (is.null(pl)) {
    if (freq > 1)
      pl = trunc(freq * 1.5)
    else pl = 2
  }
  if (is.null(pu))
    pu = trunc(freq * 8)
  if (is.null(nfix))
    nfix = freq * 3
  nq = length(theta) - 1
  b = 2 * pi/pl
  a = 2 * pi/pu
  xname = deparse(substitute(x))
  xo = x
  x = as.matrix(x)
  n = nrow(x)
  nvars = ncol(x)
  if (n < 5)
    warning("# of observations < 5")
  if (n < (2 * nq + 1))
    stop("# of observations must be at least 2*q+1")
  if (pu <= pl)
    stop("pu must be larger than pl")
  if (pl < 2) {
    warning("pl less than 2 , reset to 2")
    pl = 2
  }
  if (root != 0 && root != 1)
    stop("root must be 0 or 1")
  if (drift < 0 || drift > 1)
    stop("drift must be 0 or 1")
  if ((type == "fixed" || type == "baxter-king") && nfix <
      1)
    stop("fixed lag length must be >= 1")
  if (type == "fixed" & nfix < nq)
    stop("fixed lag length must be >= q")
  if ((type == "fixed" || type == "baxter-king") && nfix >=
      n/2)
    stop("fixed lag length must be < n/2")
  if (type == "trigonometric" && (n - 2 * floor(n/2)) != 0)
    stop("trigonometric regressions only available for even n")
  theta = as.matrix(theta)
  m1 = nrow(theta)
  m2 = ncol(theta)
  if (m1 > m2)
    th = theta
  else th = t(theta)
  g = convolve(th, th, type = "open")
  cc = g[(nq + 1):(2 * nq + 1)]
  j = 1:(2 * n)
  B = as.matrix(c((b - a)/pi, (sin(j * b) - sin(j * a))/(j *
                                                           pi)))
  R = matrix(0, n, 1)
  if (nq > 0) {
    R0 = B[1] * cc[1] + 2 * t(B[2:(nq + 1)]) * cc[2:(nq +
                                                       1)]
    R[1] = pi * R0
    for (i in 2:n) {
      dj = Bge(i - 2, nq, B, cc)
      R[i] = R[i - 1] - dj
    }
  }
  else {
    R0 = B[1] * cc[1]
    R[1] = pi * R0
    for (j in 2:n) {
      dj = 2 * pi * B[j - 1] * cc[1]
      R[j] = R[j - 1] - dj
    }
  }
  AA = matrix(0, n, n)
  if (type == "asymmetric") {
    if (nq == 0) {
      for (i in 1:n) {
        AA[i, i:n] = t(B[1:(n - i + 1)])
        if (root)
          AA[i, n] = R[n + 1 - i]/(2 * pi)
      }
      AA[1, 1] = AA[n, n]
      AAu = AA
      AAu[!upper.tri(AAu)] <- 0
      AA = AA + .flipud_adapted(.fliplr_adapted(AAu))
    }
    else {
      A = Abuild(n, nq, g, root)
      Ainv = solve(A)
      for (np in 0:ceiling(n/2 - 1)) {
        d = matrix(0, n, 1)
        ii = 0
        for (jj in (np - root):(np + 1 + root - n)) {
          ii = ii + 1
          d[ii] = Bge(jj, nq, B, cc)
        }
        if (root == 1)
          d[n - 1] = R[n - np]
        Bhat = Ainv %*% d
        AA[np + 1, ] = t(Bhat)
      }
      AA[(ceiling(n/2) + 1):n, ] = .flipud_adapted(.fliplr_adapted(AA[1:floor(n/2),
      ]))
    }
  }
  if (type == "symmetric") {
    if (nq == 0) {
      for (i in 2:ceiling(n/2)) {
        np = i - 1
        AA[i, i:(i + np)] = t(B[1:(1 + np)])
        if (root)
          AA[i, i + np] = R[np + 1]/(2 * pi)
        AA[i, (i - 1):(i - np)] = AA[i, (i + 1):(i +
                                                   np)]
      }
      AA[(ceiling(n/2) + 1):n, ] = .flipud_adapted(.fliplr_adapted(AA[1:floor(n/2),
      ]))
    }
    else {
      for (np in nq:ceiling(n/2 - 1)) {
        nf = np
        nn = 2 * np + 1
        A = Abuild(nn, nq, g, root)
        Ainv = solve(A)
        d = matrix(0, nn, 1)
        ii = 0
        for (jj in (np - root):(-nf + root)) {
          ii = ii + 1
          d[ii] = Bge(jj, nq, B, cc)
        }
        if (root)
          d[nn - 1] = R[nf + 1]
        Bhat = Ainv %*% d
        AA[np + 1, 1:(2 * np + 1)] = t(Bhat)
      }
      AA[(ceiling(n/2) + 1):n, ] = .flipud_adapted(.fliplr_adapted(AA[1:floor(n/2),
      ]))
    }
  }
  if (type == "fixed") {
    if (nq == 0) {
      bb = matrix(0, 2 * nfix + 1, 1)
      bb[(nfix + 1):(2 * nfix + 1)] = B[1:(nfix + 1)]
      bb[nfix:1] = B[2:(nfix + 1)]
      if (root) {
        bb[2 * nfix + 1] = R[nfix + 1]/(2 * pi)
        bb[1] = R[nfix + 1]/(2 * pi)
      }
      for (i in (nfix + 1):(n - nfix)) AA[i, (i - nfix):(i +
                                                           nfix)] = t(bb)
    }
    else {
      nn = 2 * nfix + 1
      A = Abuild(nn, nq, g, root)
      Ainv = solve(A)
      d = matrix(0, nn, 1)
      ii = 0
      for (jj in (nfix - root):(-nfix + root)) {
        ii = ii + 1
        d[ii] = Bge(jj, nq, B, cc)
      }
      if (root)
        d[nn - 1] = R[nn - nfix]
      Bhat = Ainv %*% d
      for (ii in (nfix + 1):(n - nfix)) AA[ii, (ii - nfix):(ii +
                                                              nfix)] = t(Bhat)
    }
  }
  if (type == "baxter-king") {
    bb = matrix(0, 2 * nfix + 1, 1)
    bb[(nfix + 1):(2 * nfix + 1)] = B[1:(nfix + 1)]
    bb[nfix:1] = B[2:(nfix + 1)]
    bb = bb - sum(bb)/(2 * nfix + 1)
    for (i in (nfix + 1):(n - nfix)) AA[i, (i - nfix):(i +
                                                         nfix)] = t(bb)
  }
  if (type == "trigonometric") {
    jj = 1:(n/2)
    jj = jj[((n/pu) <= jj & jj <= (n/pl) & jj < (n/2))]
    if (!any(jj))
      stop("frequency band is empty in trigonometric regression")
    om = 2 * pi * jj/n
    if (pl > 2) {
      for (t in 1:n) {
        for (k in n:1) {
          l = t - k
          tmp = sum(cos(om * l))
          AA[t, k] = tmp
        }
      }
    }
    else {
      for (t in 1:n) {
        for (k in n:1) {
          l = t - k
          tmp = sum(cos(om * l))
          tmp2 = (cos(pi * (t - l)) * cos(pi * t))/2
          AA[t, k] = tmp + tmp2
        }
      }
    }
    AA = AA * 2/n
  }
  if (root) {
    tst = max(abs(c(apply(AA, 1, sum))))
    if ((tst > 1e-09) && root) {
      warning("Bhat does not sum to 0 ")
      cat("test =", tst, "\n")
    }
  }
  if (drift)
    x = undrift(x)
  x.cycle = AA %*% as.matrix(x)
  if (type == "fixed" || type == "symmetric" || type == "baxter-king") {
    if (nfix > 0)
      x.cycle[c(1:nfix, (n - nfix + 1):n)] = NA
  }
  x.trend = x - x.cycle
  if (is.ts(xo)) {
    tsp.x = tsp(xo)
    x.cycle = ts(x.cycle, start = tsp.x[1], frequency = tsp.x[3])
    x.trend = ts(x.trend, start = tsp.x[1], frequency = tsp.x[3])
    x = ts(x, start = tsp.x[1], frequency = tsp.x[3])
  }
  if (type == "asymmetric")
    title = "Chiristiano-Fitzgerald Asymmetric Filter"
  if (type == "symmetric")
    title = "Chiristiano-Fitzgerald Symmetric Filter"
  if (type == "fixed")
    title = "Chiristiano-Fitzgerald Fixed Length Filter"
  if (type == "baxter-king")
    title = "Baxter-King Fixed Length Filter"
  if (type == "trigonometric")
    title = "Trigonometric Regression Filter"
  res <- list(cycle = x.cycle, trend = x.trend, fmatrix = AA,
              title = title, xname = xname, call = as.call(match.call()),
              type = type, pl = pl, pu = pu, nfix = nfix, root = root,
              drift = drift, theta = theta, method = "cffilter", x = x)
  return(structure(res, class = "mFilter"))
}

.flipud_adapted <- function (x)
{
  apply(as.matrix(x), 2, rev)
}

.fliplr_adapted <- function (x)
{
  t(apply(as.matrix(x), 1, rev))
}






#' Plot listed plots with one common legend
#'
#' @param plot_lists The list of ggplots
#' @param nrow Numeric, number of rows to arrange the plots
#' @param ncol Numeric, number of columns to arrange the plots.
#' @param position Character, Position of the legend, either "bottom" or "right".
#'
#' @return combained ggplot
#' @export
#'
#' @import gridExtra, ggplot2
#'
#' @examples
#' gg_a <- grid_arrange_common_legend(plot_lists = plots, nrow = 2, ncol = 2,
#' position = "right")

grid_arrange_common_legend <- function(plot_lists,
                                       nrow = 1,
                                       ncol = length(plot_lists),
                                       position = c("bottom", "right")) {

  position <- match.arg(position)
  g <- ggplotGrob(plot_lists[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plot_lists, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = grid::unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = grid::unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid::grid.newpage()
  grid::grid.draw(combined)
  return(combined)

}



#' Plot 2D scatter plots
#'
#' @description Plots biaxial plots
#'
#' @param fcs_files Character, pathway to fcs files.
#' @param markers_to_plot Character, pattern of the markers to be plotted e.g.
#' "CD" (all CD markers will be plotted), "CD41$" (only CD41 will be plotted).
#' @param y_marker Character, The marker to be plotted on the y-axis.
#' @param out_dir Character, path where fill are saved, if NULL (default)
#' files are saved in getwd()
#' @param out_put_name The name of the file that is saved in out_dir, default is
#' "marker_plots.png".
#'
#' @return Save plots in out_dir.
#' @export
#'
#' @examples
#' plot_2D_scatter_plots(fcs_files = fcs_files,
#'                       markers_to_plot = c("CD", "HLA"),
#'                       y_marker = "CD45$",
#'                       out_dir = getwd(),
#'                       out_put_name = "marker_plots.png")
#'
#'
plot_2D_scatter_plots <- function(fcs_files,
                                  markers_to_plot,
                                  y_marker,
                                  out_dir = NULL,
                                  out_put_name = "marker_plots.png"){


  if(!all(file.exists(fcs_files))){
    stop("fcs files do not exist")
  }

  ff <- flowCore::read.FCS(fcs_files[1], transformation = FALSE)

  if (!is.null(markers_to_plot)){

    if(!is.character(markers_to_plot)){
      stop ("markers are not a character vector")
    }

    matches <- paste(markers_to_plot, collapse="|")

    norm_markers <- grep(matches,
                         FlowSOM::GetMarkers(ff, find_mass_ch(ff,
                                                              value = TRUE)),
                         value = TRUE, ignore.case = FALSE)

  } else {
    norm_markers <- find_mass_ch(ff, value = TRUE)
    norm_markers <- FlowSOM::GetMarkers(ff, norm_markers)
  }

  matches_y <- paste(y_marker, collapse="|")
  ym <- grep(matches_y,
             FlowSOM::GetMarkers(ff, find_mass_ch(ff,
                                                  value = TRUE)),
             value = TRUE, ignore.case = FALSE)

  n_plots <- length(norm_markers) - 1
  png(file.path(out_dir, paste0(out_put_name)),
      width = n_plots * 300, height = length(fcs_files) * 300)
  layout(matrix(1:(length(fcs_files) * n_plots), ncol = n_plots, byrow = TRUE))

  for(file in fcs_files){
    ff <- read.FCS(file, transformation = FALSE)

    cols_to_trans <- grep(pattern = "Di", x = colnames(ff), value = TRUE)

    ff_t <- flowCore::transform(ff, flowCore::transformList(cols_to_trans, cytofTransform))

    for(m in names(norm_markers)){

      if(m != names(ym)){
        plotDens(obj = ff_t, channels = c(m,names(ym)), main = basename(file))
        print(paste("plotting", basename(file), m))
      }
    }
  }
  dev.off()
}


#' Title
#' 
#' @description Perform normalization using reference files. Takes advantage of
#' CytoNorm package. 
#'
#' @param df Data frame containing following columns:
#' file_paths (the full path to the files to be normalized), 
#' batch_label (batch label for each file), ref_ids (logical defining TRUE values
#' for reference sample).
#' @param markers_to_normalize Character vector, marker names to be normalized, 
#' can be full marker name e.g. "CD45$" (only CD45 marker will be picked) or 
#' "CD" (all markers containig "CD" will be used).
#' If NULL (default) all non-mass markers will be normalized.
#' @param transformList 
#' @param arcsine_transform 
#' @param nQ 
#' @param limit 
#' @param quantileValues 
#' @param goal 
#' @param to_plot 
#' @param norm_with_clustering 
#' @param seed 
#' @param nCells 
#' @param xdim 
#' @param ydim 
#' @param nClus 
#' @param clustering_markers 
#' @param out_dir 
#' @param save_model 
#'
#' @return
#' @export
#'
#' @examples
train_REF_model <- function(df,
                            markers_to_normalize = NULL,
                            transformList = NULL,
                            arcsine_transform = TRUE,
                            nQ = 101,
                            limit = NULL,
                            quantileValues = NULL,
                            goal = "mean",
                            to_plot = TRUE,
                            norm_with_clustering = FALSE, 
                            seed = NULL, 
                            nCells = 10000,
                            xdim = 10,
                            ydim = 10,
                            nClus = 10,
                            clustering_markers = NULL, 
                            out_dir = NULL, 
                            save_model = FALSE){
  
  if(!is(df, "data.frame")){
    stop("df is not a data frame")
  }
  
  if(!all(colnames(df) == c("file_paths", "batch_labels", "ref_ids"))){
    stop("colnames in df does not match: file_paths, batch_ids, ref_ids 
        please correct colnames")
  }
  
  if(!all(file.exists(df$file_paths))){
    id <- !file.exists(df$file_paths)
    print(df$file_paths[id])
    stop("above files have incorrect path")
  }
 
  
  flow_frame <- flowCore::read.FCS(df$file_paths[1])
  if (!is.null(markers_to_normalize)){
    
    matches <- paste(markers_to_normalize, collapse="|")
    
    m_to_keep <- names(grep(matches, FlowSOM::GetMarkers(flow_frame, flowCore::colnames(flow_frame)),
                            ignore.case = TRUE, value = TRUE))
    m_to_keep <- grep(pattern = "File", x = m_to_keep, ignore.case = TRUE,
                      invert = TRUE, value = TRUE)
    
  } else {
    m_to_keep <- grep(pattern = "Time|length|Center|Offset|Width|Residual|File|File_scattered",
                      x = flowCore::colnames(flow_frame),
                      ignore.case = TRUE, value = TRUE, invert = TRUE)
  }
  
  # create out_dir if does not exist
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "CytoNormed")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  
  files_ref <- df$file_paths[df$ref_ids]
  message(paste("The following reference files were found"))
  print(basename(files_ref))
  
  labels_ref <- df$batch_labels[df$ref_ids]
  
  if(arcsine_transform & !is.null(transformList)){
    stop("define either arcsine_transform or transformList parameter")
  } else if(arcsine_transform==FALSE & is.null(transformList)){
    stop("define either arcsine_transform or transformList parameter")
  }
  
  if(arcsine_transform){
    trans <- flowCore::transformList(m_to_keep,
                                            CytoNorm::cytofTransform)
   
  } else {
    trans <- transformList
  }
  
  if(to_plot & !norm_with_clustering){
    
    png(file.path(out_dir, "REF_normalization_taining_model.png"),
        width = length(m_to_keep) * 300,
        height = (length(files_ref) * 2 + 1) * 300)
  }
  
  
  if(norm_with_clustering){
    print("Clustering the data using FlowSOM")
    
    if(!is.null(clustering_markers)){
      markers_all <- FlowSOM::GetMarkers(flow_frame, flowCore::colnames(flow_frame))
      clustering_pattern <- paste(clustering_markers, collapse = "|")
      clustering_pattern <- names(grep(clustering_pattern, markers_all, value = TRUE))
    } else {
      clustering_pattern <- m_to_keep
    }

    
    nC <- length(files_ref)*nCells
    model <- CytoNorm::CytoNorm.train(files = files_ref, 
                                      labels = labels_ref,
                                      channels = m_to_keep, 
                                      transformList = trans, 
                                      plot = to_plot, 
                                      seed = seed, 
                                      normParams = list(nQ = nQ,
                                                        limit = limit,
                                                        quantileValues = quantileValues,
                                                        goal = goal),
                                      FlowSOM.params = list(nCells = nC,
                                                            xdim = xdim,
                                                            ydim = ydim,
                                                            nClus = nClus,
                                                            scale = FALSE,
                                                            colsToUse = clustering_pattern), 
                                      outputDir = out_dir)
    
  } else {
    model <- CytoNorm::QuantileNorm.train(files = files_ref,
                                          labels = labels_ref,
                                          channels = m_to_keep,
                                          transformList = trans,
                                          nQ = 2,
                                          limit = c(0,8),
                                          quantileValues = c(0.05, 0.95),
                                          goal = "mean",
                                          plot = to_plot)
  }
  
  
  
  if(to_plot & !norm_with_clustering){
    dev.off()
  }
  print(out_dir)
  if (save_model){
    saveRDS(model, file.path(out_dir, "REF_normalization_model.RSD"))
  }
  return(model)
}


normalize_REF <- function(model, 
                          df,
                          arcsine_transform = TRUE,
                          transformList = NULL,
                          transformList.reverse = NULL,
                          out_dir = NULL, 
                          norm_with_clustering = FALSE){
  
  files <- df$file_paths
  
  if(!all(file.exists(files))){
    stop("the files does not exist, please specify the corect pathway")
  }
  
  labels <- df$batch_labels
  
  if(arcsine_transform & !is.null(transformList)){
    stop("define either arcsine_transform or transformList parameter")
  } else if(arcsine_transform==FALSE & is.null(transformList)){
    stop("define either arcsine_transform or transformList parameter")
  }
  
  if(arcsine_transform){
    flow_frame <- flowCore::read.FCS(df$file_paths[1])
    
    
    trans <- flowCore::transformList(flowCore::colnames(flow_frame),
                                     CytoNorm::cytofTransform)
    trans_rev <- flowCore::transformList(flowCore::colnames(flow_frame),
                                         CytoNorm::cytofTransform.reverse)
    
  } else {
    if(is.null(transformList)){
      stop("define transformList")
    }
    
    trans <- transformList
    
    if(is.null(transformList.reverse)){
      stop("define transformList")
    }
    trans_rev <- transformList.reverse
  }
  
  if(is.null(out_dir)){
    out_dir <- file.path(getwd(), "CytoNormed")
  }
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  
  if(norm_with_clustering){
    CytoNorm::CytoNorm.normalize(model = model,
                                 files = files, 
                                 labels = labels, 
                                 transformList = trans, 
                                 transformList.reverse = trans_rev, 
                                 outputDir = out_dir)
    
  } else {
    CytoNorm::QuantileNorm.normalize(model = model,
                                     files = files,
                                     labels = labels,
                                     transformList = trans,
                                     transformList.reverse = trans_rev,
                                     outputDir = out_dir)
  }
}
