make_progress <- function(max) {
  message("0%   10   20   30   40   50   60   70   80   90   100%")
  message("[----|----|----|----|----|----|----|----|----|----|")

  max_stars <- 51 # length of the progress bar
  value <- 0
  curr_stars <- 0
  list(
    value = value,
    curr_stars = curr_stars,
    msm = max_stars / max,
    max_stars = max_stars
  )
}

increment_progress <- function(progress) {
  if (progress$curr_stars >= progress$max_stars) {
    return(progress)
  }
  progress$value <- progress$value + 1
  num_stars <- round(progress$msm * progress$value)
  if (num_stars > progress$curr_stars) {
    # Number of new stars to print
    num_new_stars <- num_stars - progress$curr_stars

    # If we are going to reach the end of the progress bar
    # save space for the terminal "|"
    if (num_stars >= progress$max_stars) {
      num_new_stars <- num_new_stars - 1
    }
    new_stars <- paste(rep("*", num_new_stars), collapse = "")

    message(new_stars, appendLF = FALSE)
    utils::flush.console()
    progress$curr_stars <- num_stars
  }
  if (progress$curr_stars >= progress$max_stars) {
    # The terminal "|" character that appears instead of a *
    message("|")
  }
  progress
}
