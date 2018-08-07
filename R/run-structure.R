# runStructure.R

#' Run STRUCTURE in current path
#'
#' @param path_to_structure path to structure binary executable (ie. /usr/bin/structure)
#' @param input_file file name of input data
#' @param main_params file name of mainparams file for STRUCTURE
#' @param extra_params file name of extraparams file for STRUCTURE
#' @param out_prefix prefix path/name for logging
#' @param n_K number of assumed populations to try
#' @param n_replicates number of replicates
#' @param n_cores number of cores
#' @note Set RANDOMIZE = 0 in main params file to avoid using same seed. Haven't tested
#' on Windows.
#' @importFrom stringr str_pad str_extract
#' @importFrom parallel mcmapply detectCores
#' @importFrom stats runif
#' @export
#' @examples
#'\dontrun{
#'input_file <- system.file("inst/extdata/microsat_testfiles", "locprior.str", package = "starmie")
#'main_params <- system.file("inst/extdata/microsat_testfiles", "mainparams", package = "starmie")
#'extra_params <-  system.file("inst/extdata/microsat_testfiles", "extraparams", package = "starmie")
#'runStructure("structure", input_file, main_params, extra_params, "test", 5, 2, 2)
#'}
runStructure <- function(path_to_structure, input_file, main_params, extra_params,
                         out_prefix, n_K, n_replicates, n_cores) {
  # i/o checks
  # check R can find files
  stopifnot(file.exists(path_to_structure))
  stopifnot(file.exists(input_file))
  stopifnot(file.exists(main_params))
  stopifnot(file.exists(extra_params))

  if (!is.integer(n_K) & n_K < 1 ) {
    stop("Assumed populations must be greater than 1")
  }

  if (!is.integer(n_replicates) & n_replicates < 1 ) {
    stop("Number of replicates must be greater than 1")
  }

  # set up vectors for running replictes
  K <- 1L:n_K
  replicates <- 1L:n_replicates

  # file names created by STRUCTURE output
  out_files <- outer(replicates, K,
                     function(x,y) paste0(out_prefix,
                                          stringr::str_pad(x, width = 2, pad = 0),
                                          "_K_",
                                          stringr::str_pad(y, width = 2, pad = 0),
                                          ".out"))
  log_files <- gsub("out", "log", out_files)

  run_structure_single <- function(out_file, log_file) {
    k <- as.integer(stringr::str_extract(out_file, "[0-9]{2}\\b"))
    seed <- round(runif(1) * 1e+08)

    cmd <- paste(path_to_structure, "-K", k,
                 "-i", input_file,
                 "-m", main_params,
                 "-e", extra_params,
                 "-D", seed , "-o", out_file, '&>', log_file)

    # save log as character file
    system(cmd)
    return(cmd)
  }

  message(paste("Running STRUCTURE on", n_cores, "core with",
                n_K, "populations with", n_replicates, "replicates."))
  if(requireNamespace("parallel", quietly = TRUE)) {

    if (n_cores > parallel::detectCores()) {
      stop("Number of cores greater than number available on machine.")
    }

    cmds_run <- parallel::mcmapply(run_structure_single, out_files, log_files,
                                   mc.cores = n_cores, mc.set.seed = TRUE)
    message(paste("Commands run\n", paste(cmds_run, collapse = "\n")))

  } else {
    if (n_cores > 1L) {
      stop("parallel package needed to run multiple calls. Please install it.",
           call. = FALSE)
    } else {
      cmds_run <- mcmapply(run_structure_single, out_files, log_files)
      message(paste("Commands run\n", paste(cmds_run, collapse = "\n")))

    }
  }
}
