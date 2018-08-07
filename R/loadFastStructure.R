# loadFastStructure.R
# Functions for parsing output from fastSTRUCTURE

#' Read fastStructure Output
#' @param  qfile a valid Q file from fastSTRUCTURE
#' @param pfile a corresponding P file from fastSTRUCTURE
#' @param logfile optional logfile from corresponding fastSTRUCTURE run
#' @return an \link{fast} object containing the output of of an fastSTRUCTURE run
#' @importFrom data.table fread
#' @importFrom stringr str_split_fixed
#' @export
#' @examples
#' qfin <- system.file("extdata/fastStructure_files", "hapmap3.2.meanQ", package = "starmie")
#' pfin <- system.file("extdata/fastStructure_files", "hapmap3.2.meanP", package = "starmie")
#' my_fast <- loadFastStructure(qfin, pfin)
#' # add log file
#' logfin <- system.file("extdata/fastStructure_files", "hapmap3.2.log", package = "starmie")
#' my_fast <- loadFastStructure(qfin, pfin, logfile = logfin)
loadFastStructure <- function(qfile, pfile, logfile = NULL) {
  # i/o checks
  if ( !(is.character(qfile) & length(qfile) == 1) )
    stop("qfile must be a character vector of length 1.")
  if ( !(is.character(pfile) & length(pfile) == 1) )
    stop("pfile must be a character vector of length 1.")

  if ( !is.null(logfile) ) {
    if ( !(is.character(logfile) & length(logfile) == 1) )
      stop("logfile must be a character vector of length 1.")
  }

  # create new fastStructure object
  fast_obj <- fast()
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

  # read in logfiles
  if ( !is.null(logfile) ) {
    run_summary <- read_fast_logfiles(logfile)
    log_info <- run_summary$log_info
  } else {
    log_info <- NULL
  }

  # output list containing q_info, p_info and log_info

  fast_obj$K <- K
  fast_obj$nsamples <- nsamples
  fast_obj$nmarkers <- nmarkers
  fast_obj$Q_df <- q_df
  fast_obj$P_df <- p_df
  fast_obj$log_info <- log_info
  return(fast_obj)
}

read_fast_logfiles <- function(logfile) {
  # grab final logliklihood and CV error (if available) and number of bootstrap resamples
  # if applicable
  logfin <- readLines(logfile)

  loglik <- as.numeric(gsub("Marginal Likelihood = ", "",
                            logfin[grep("^Marginal Likelihood", logfin)]))

  # check for CV error
  match_cverror <- grep("^CV error", logfin)
  if (length(match_cverror) == 1) {
    cverror <- gsub("CV error = ", "", logfin[match_cverror])
    cv.sd <- as.numeric(gsub(".*, ", "", cverror))
    cverror <- as.numeric(gsub(", .*", "", cverror))
  } else {
    cverror <- NA
    cv.sd <- NA
  }

  run_summary <- list(log_info = data.frame(logL = loglik,
                                            cverror = cverror,
                                            cv.sd = cv.sd,
                                            stringsAsFactors = FALSE))
  run_summary
}
