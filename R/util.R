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

# Do work and update a progress bar
progress_for <- function(n, nchunks, fun) {
  message("0%   10   20   30   40   50   60   70   80   90   100%")
  message("[----|----|----|----|----|----|----|----|----|----|")
  remaining <- n
  chunk_end <- 0
  for (i in 1:nchunks) {
    chunk_start <- chunk_end + 1
    chunk_end <- chunk_start + round(remaining / (nchunks - i + 1)) - 1
    remaining <- remaining - (chunk_end - chunk_start + 1)

    fun(chunk_start, chunk_end)

    message("*", appendLF = FALSE)
    utils::flush.console()
  }
  message("|")
}
