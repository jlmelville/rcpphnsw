stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
# appears only if called from an environment where a logical verbose = TRUE
# OR force = TRUE
tsmessage <- function(
  ...,
  domain = NULL,
  appendLF = TRUE,
  force = FALSE,
  time_stamp = TRUE
) {
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
  check_whole_number(random_seed, "random_seed", lower = 0)
}

check_whole_number <- function(
  value,
  name,
  lower,
  upper = .Machine$integer.max
) {
  if (
    !is.numeric(value) ||
      length(value) != 1 ||
      is.na(value) ||
      !is.finite(value) ||
      value != floor(value) ||
      value < lower ||
      value > upper
  ) {
    stop(
      name,
      " cannot be outside the whole-number range ",
      lower,
      " to ",
      upper,
      call. = FALSE
    )
  }

  as.integer(value)
}

check_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  }

  value
}

check_input_matrix <- function(X, byrow, allow_empty = FALSE) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix", call. = FALSE)
  }

  if (byrow) {
    nitems <- nrow(X)
    ndim <- ncol(X)
  } else {
    nitems <- ncol(X)
    ndim <- nrow(X)
  }

  if (ndim == 0) {
    stop("X must have at least one dimension", call. = FALSE)
  }
  if (!allow_empty && nitems == 0) {
    stop("X must contain at least one item", call. = FALSE)
  }

  list(nitems = nitems, ndim = ndim)
}
