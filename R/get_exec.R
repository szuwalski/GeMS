#' Get Executable
#' Get SimAss binary/executable location in package
#'
#' @param exe_name Name of the SimAss binary
#'
#' @return The path to a GeMS Age-Structured Model binary.
#'
#' @export
#' @examples
#' \dontrun{
#' get_exec("simass.exe")
#' }

get_exec <- function(exe_name=SimAssExec) {
  exec.path <- system.file("src", package = "GeMS")
  full.path <- file.path(exec.path,exe_name)
  if(!file.exists(full.path)) stop("Executable file not found. Please download and install package from https://github.com/szuwalski/GeMS")
  return(full.path)
}
