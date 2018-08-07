# loadAdmixture.R
# Functions for parsing output from admixture

#' Read Admixture Output
#' @param  qfile a valid Q file from ADMIXTURE
#' @param pfile a corresponding P file from ADMXIXTURE
#' @param logfile optional logfile from corresponding ADMIXTURE run
#' @param bootstrap_bias optional file containing bias estimates for Q matrix
#' @param bootstrap_se optional file containing standard error estiamtes for Q matrix
#' @return an \link{admix} object containing the output of of an admixture run
#' @importFrom data.table fread
#' @importFrom stringr str_split_fixed
#' @export
#' @examples
#' qfin <- system.file("extdata/hapmap3_files", "hapmap3.2.Q", package = "starmie")
#' pfin <- system.file("extdata/hapmap3_files", "hapmap3.2.P", package = "starmie")
#' my_admix <- loadAdmixture(qfin, pfin)
#' # add log file
#' logfin <- system.file("extdata/hapmap3_files", "log2.out", package = "starmie")
#' my_admix <- loadAdmixture(qfin, pfin, logfile = logfin)
loadAdmixture <- function(qfile, pfile,logfile = NULL, bootstrap_bias = NULL, bootstrap_se = NULL) {
  # i/o checks
  if ( !(is.character(qfile) & length(qfile) == 1) )
    stop("qfile must be a character vector of length 1.")
  if ( !(is.character(pfile) & length(pfile) == 1) )
    stop("pfile must be a character vector of length 1.")
  if ( !(is.null(bootstrap_bias)) ) {
    if( !(is.character(bootstrap_bias) & length(bootstrap_bias) == 1))
      stop("bootstrap_bias must be a character vector of length 1.")
  }

  if ( !(is.null(bootstrap_bias)) ) {
    if( !(is.character(bootstrap_bias) & length(bootstrap_bias) == 1))
      stop("bootstrap_bias must be a character vector of length 1.")
  }

  if ( !(is.null(bootstrap_se)) ) {
    if ( !(is.character(bootstrap_se) & length(bootstrap_se) == 1))
      stop("bootsrap_se must be a character vector of length 1.")
  }

  if ( !is.null(logfile) ) {
    if ( !(is.character(logfile) & length(logfile) == 1) )
      stop("logfile must be a character vector of length 1.")
  }

  # create new admixture object
  admix_obj <- admix()
  # qfile data reader
  q_df <- fread(qfile, header = FALSE, data.table = FALSE)
  K <- ncol(q_df)
  nsamples <- nrow(q_df)
  colnames(q_df) <- paste("Cluster", seq(1, K))
  # pfile data reader
  p_df <- fread(pfile, header = FALSE, data.table = FALSE)
  nmarkers <- nrow(p_df)
  if(ncol(p_df) != K) {
    stop("Number of populations does not match between Q and P files")
  }

  if ( !is.null(bootstrap_bias) ) {
    Q_bias_df <- fread(bootstrap_bias, header = FALSE, data.table = FALSE)
    if (ncol(Q_bias_df) != K)
      stop("Number of populations does not match between admixture inputs.")
    colnames(Q_bias_df) <- colnames(q_df)
  } else {
    Q_bias_df <- NULL
  }

  if ( !is.null(bootstrap_se) ) {
    Q_se_df <- fread(bootstrap_se, header = FALSE, data.table = FALSE)
    if (ncol(Q_se_df) != K)
      stop("Number of populations does not match between admixture inputs.")
    colnames(Q_se_df) <- colnames(q_df)
  } else {
    Q_se_df <- NULL
  }

  # read in logfiles
  if ( !is.null(logfile) ) {
    run_summary <- read_logfiles(logfile)
    if(run_summary$log_info$K != K) {
      stop("logfile does not match input P and Q files")
    }
    log_info <- run_summary$log_info
    fst_matrix <- run_summary$fst_matrix
  } else {
    log_info <- NULL
    fst_matrix <- NULL
  }

  # output list containing q_info, p_info and log_info

  admix_obj$K <- K
  admix_obj$nsamples <- nsamples
  admix_obj$nmarkers <- nmarkers
  admix_obj$Q_df <- q_df
  admix_obj$P_df <- p_df
  admix_obj$Q_bias_df <- Q_bias_df
  admix_obj$Q_se_df <- Q_se_df
  admix_obj$fst_matrix <- fst_matrix
  admix_obj$log_info <- log_info
  return(admix_obj)
}

read_logfiles <- function(logfile) {
  # grab final logliklihood and CV error (if available) and number of bootstrap resamples
  # if applicable
  logfin <- readLines(logfile)
  seed <- as.integer(gsub("Random seed:", "",
                          logfin[grep("^Random seed", logfin)]))
  loglik <- as.numeric(gsub("Loglikelihood: ", "",
                            logfin[grep("^Loglikelihood:", logfin)]))

  # check for CV error
  match_cverror <- grep("^CV error", logfin)
  if (length(match_cverror) == 1) {
    cverror <- as.numeric(gsub(".*: ", "\\1", logfin[match_cverror]))
  } else {
    cverror <- NA
  }
  # check for bootstrap resamples
  match_bootstraps <- grep("\\bbootstrap resamplings\\b", logfin)
  if (length(match_bootstraps) == 1) {
    bootstrap_replicates <- as.integer(sub(".*\\s+([[:digit:]]+)\\s+.*",
                                "\\1", logfin[match_bootstraps]))
  } else {
    bootstrap_replicates <- NA
  }

  # read Fst info
  fst_lines <- logfin[grep("^Pop[0-9]{1,}",logfin)]

  if (length(fst_lines) < 1) {
    message("No Fst lines found, assuming K=1")
    K <- 1
    fst_matrix <- matrix(NA, 1,1)
    rownames(fst_matrix) <- "Cluster 1"
  } else {
    K <- length(fst_lines)
    fst_matches <- stringr::str_split_fixed(fst_lines, "\\t", n = K + 1)
    fst_matrix <- fst_matches[, 2:(K+1)]
    diag(fst_matrix) <- 0
    class(fst_matrix) <- "numeric"
    rownames(fst_matrix) <- paste("Cluster", seq(1,K))

  }


  run_summary <- list(log_info = data.frame(K = K,
                                            logL = loglik,
                                            CVerror = cverror,
                                            bootstrap_replicates = bootstrap_replicates,
                                            stringsAsFactors = FALSE),
                      fst_matrix = fst_matrix)
  run_summary
}
