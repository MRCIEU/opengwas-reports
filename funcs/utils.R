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
  date <- system(glue('git show -s --format=%ci {commit}'), intern = TRUE) %>%
    lubridate::as_date()
  glue("{branch}-{commit}-{date}")
}
