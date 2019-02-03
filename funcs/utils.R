verify_path <- function(path, how = c("fail", "warning")) {
  #' Verify `path` and fail or raise warning
  how = match.arg(how)
  if (how == "fail") {
    if (!file_exists(path)) {
      stop(glue("File or directory not exists: {path}"))
      quit("no")
    }
  } else {
    if (!file_exists(path)) {
      warning(glue("File or directory not exists: {path}"))
    }
  }
}

reuse_or_process <- function(path, func, args, reuse = FALSE) {
  #' If `file` exists and `reuse`, reuse it,
  #' otherwise generate this file according to (apply func args)
  message(glue("Processing {path}..."))
  if (file_exists(path) && reuse) {
    message(glue("Reuse {path}"))
  } else {
    message(glue("Generating {path}..."))
    do.call(what = func, args = args)
  }
  message(glue("Processing {path} finished"))
}

get_build_version <- function() {
  # TODO: get build version in the format of "<branch>-<hash>-<date>"
}
