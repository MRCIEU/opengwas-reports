linters <- lintr::with_defaults(
  absolute_paths_linter = NULL,
  absolute_path_linter = NULL,
  infix_spaces_linter = NULL,
  camel_case_linter = NULL,
  object_usage_linter = NULL,
  spaces_left_parentheses_linter = NULL,
  open_curly_linter = NULL,
  commented_code_linter = NULL
)

for (file in list.files(
  pattern = "[.]R(md)?$",
  ignore.case = TRUE,
  recursive = TRUE
)) {
  print(lintr::lint(filename = file, linters = linters))
}
