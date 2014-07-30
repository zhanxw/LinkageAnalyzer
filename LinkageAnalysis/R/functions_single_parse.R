get_run <- function(dir, run, filenames) {
  # load file
  results_file <- file.path(dir, run, "results.RData")
  suppressWarnings(rm("analysis"))
  tryCatch(load(results_file), error = function(e) report("e", "Record doesn't exist!",
    filenames$log_file))
  if (length(get("analysis")) != 8) {
    report("e", paste("Incomplete record for ", run, "!", sep = ""), filenames$log_file)
  }
  return(get("analysis"))
}

get_all_runs <- function(dir, runs, filenames, silent) {
  # get each run
  all_results <- c()
  for (run in runs) {
    if (silent == FALSE) {
      report("m", paste("Parsing:", run), filenames$log_file)
    }
    all_results[[run]] <- get_run(dir, run, filenames)
  }

  # check consistency
  same_file <- all_results[[run]]$results[, 1:3]  # the last run

  # check each run
  for (run in runs) {
    each_file <- all_results[[run]]$results[, 1:3]
    if (dim(each_file)[1] != dim(same_file)[1] || !all(same_file == each_file)) {
      report("e", "The analyses are not conducted on the same genes!", filenames$log_file)
    }
  }

  return(all_results)
}
