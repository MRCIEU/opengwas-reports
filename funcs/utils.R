verify_path <- function(path, how = c("fail", "warning")) {
  #' Verify `path` and fail or raise warning
  how <- match.arg(how)
  if (how == "fail") {
    if (!file_exists(path)) {
      stop(glue("File or directory not exists: {path}"))
      quit("no")
    }
  } else {
    if (!file_exists(path)) {
      logwarn(glue("File or directory not exists: {path}"))
    }
  }
}

reuse_or_process <- function(path, func, args, reuse = FALSE) {
  #' If `file` exists and `reuse`, reuse it,
  #' otherwise generate this file according to (apply func args)
  loginfo(glue("Processing {path}..."))
  if (file_exists(path) && reuse) {
    loginfo(glue("Reuse {path}"))
  } else {
    loginfo(glue("Generating {path}..."))
    do.call(what = func, args = args)
  }
  loginfo(glue("Processing {path} finished"))
}

get_build_version <- function() {
  #' get build version in the format of "<branch>-<hash>-<date>"
  branch <- system("git branch | grep \\* | cut -d ' ' -f2", intern = TRUE)
  commit <- system("git rev-parse --verify --short HEAD", intern = TRUE)
  date <- system(glue("git show -s --format=%ci {commit}"), intern = TRUE) %>%
    lubridate::as_date()
  glue("{branch}-{commit}-{date}")
}

split_by_chunk <- function(vec, n_chunks, chunk_idx) {
  chunks <- function(vec, n_chunks) {
    chunk_size <- ceiling(length(vec) / n_chunks)
    1:n_chunks %>%
      map(function(idx, vec, chunk_size) {
        vec[(1 + (idx - 1) * chunk_size):
        min((idx - 1) * chunk_size + chunk_size, length(vec))]
      },
      vec = vec, chunk_size = chunk_size
      )
  }
  vec_split <- chunks(vec, n_chunks)
  vec_idx <- vec_split[[chunk_idx]]
  vec_idx
}

neg_log10 <- function(pval, is_neg_log10 = FALSE) {
  #' If not `is_neg_log10`, transform by -log10(pval),
  #' otherwise return as is.
  if (!is_neg_log10) {
    res <- -log10(pval)
  } else {
    res <- pval
  }
  res
}

unlog <- function(val, is_log = FALSE,
                  func = function(x) 10^-x) {
  #' If `is_log`, restore from log transformed value, by `func`
  if (is_log) {
    val <- func(val)
  }
  val
}
