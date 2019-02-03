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
