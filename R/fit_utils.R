sum_by_group <- function(
    x,
    group,
    n_groups = max(group)
) {

  out <- numeric(n_groups)

  tmp <- rowsum(
    as.numeric(x),
    group = group,
    reorder = FALSE
  )

  index <- as.integer(
    rownames(tmp)
  )

  out[index] <- as.numeric(
    tmp[, 1]
  )

  out
}


prepare_frailty_fit <- function(
    sim
) {

  required <- c(
    "long",
    "subjects",
    "n_pairs_observed",
    "pair_zig"
  )

  missing <- setdiff(
    required,
    names(sim)
  )

  if (length(missing) > 0L) {
    stop(
      "simulation object is missing: ",
      paste(missing, collapse = ", ")
    )
  }

  long <- sim$long
  subjects <- sim$subjects

  required_long <- c(
    "pair_index",
    "id",
    "zig",
    "start",
    "stop",
    "status",
    "y"
  )

  required_subjects <- c(
    "pair_index",
    "id",
    "zig",
    "entry",
    "stop",
    "status",
    "y_entry"
  )

  missing_long <- setdiff(
    required_long,
    names(long)
  )

  missing_subjects <- setdiff(
    required_subjects,
    names(subjects)
  )

  if (length(missing_long) > 0L) {
    stop(
      "long data are missing: ",
      paste(missing_long, collapse = ", ")
    )
  }

  if (length(missing_subjects) > 0L) {
    stop(
      "subject data are missing: ",
      paste(missing_subjects, collapse = ", ")
    )
  }

  if (any(long$stop <= long$start)) {
    stop("all start-stop intervals must have positive length")
  }

  if (any(!long$status %in% c(0L, 1L))) {
    stop("status must contain only 0 and 1")
  }

  n_pairs <- sim$n_pairs_observed
  n_subjects <- nrow(subjects)

  if (length(sim$pair_zig) != n_pairs) {
    stop("pair_zig is incompatible with n_pairs_observed")
  }

  long_subject_index <- match(
    long$id,
    subjects$id
  )

  if (anyNA(long_subject_index)) {
    stop("some long-format subjects are missing from subject data")
  }

  pair_index_long <- long$pair_index
  pair_index_subject <- subjects$pair_index

  d_pair <- sum_by_group(
    long$status,
    pair_index_long,
    n_pairs
  )

  d_subject <- sum_by_group(
    long$status,
    long_subject_index,
    n_subjects
  )

  list(
    long = long,
    subjects = subjects,
    n_pairs = n_pairs,
    n_subjects = n_subjects,
    pair_zig = sim$pair_zig,
    dz_index = which(
      sim$pair_zig == "DZ"
    ),
    mz_index = which(
      sim$pair_zig == "MZ"
    ),
    pair_index_long = pair_index_long,
    pair_index_subject = pair_index_subject,
    long_subject_index = long_subject_index,
    d_pair = d_pair,
    d_subject = d_subject,
    start = long$start,
    stop = long$stop,
    log_stop = log(long$stop),
    status = long$status,
    y = long$y,
    entry = subjects$entry,
    y_entry = subjects$y_entry
  )
}
