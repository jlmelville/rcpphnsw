stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
# appears only if called from an environment where a logical verbose = TRUE
# OR force = TRUE
tsmessage <- function(..., domain = NULL, appendLF = TRUE, force = FALSE,
                      time_stamp = TRUE) {
  verbose <- get0("verbose", envir = sys.parent())

  if (force || (!is.null(verbose) && verbose)) {
    msg <- ""
    if (time_stamp) {
      msg <- paste0(stime(), " ")
    }
    message(msg, ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}

check_random_seed <- function(random_seed) {
  seed_max <- .Machine$integer.max

  if (!is.numeric(random_seed) ||
    length(random_seed) != 1 ||
    is.na(random_seed) ||
    random_seed < 0 ||
    random_seed > seed_max ||
    random_seed != floor(random_seed)) {
    stop(
      "random_seed must be an integer between 0 and ",
      seed_max,
      call. = FALSE
    )
  }

  as.integer(random_seed)
}
