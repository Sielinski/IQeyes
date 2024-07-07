#######################
## read_pentacam_csv
#######################

#' Read a Pentacam CSV file
#' @description
#' Reads a Pentacam CSV file into a data frame. The function also converts the
#' exam date and time columns to lubridate date and period fields, respectively,
#' converts the patient's date of birth to a date field (and calculates their
#' age at the time of the exam), converts the exam eye to a factor, and converts
#' the file's column names to an R-friendly format.
#' @param file_name
#' A character string containing the name of the file name to read.
#' @param delimiter
#' A character identifying the delimiter that the Pentacam file uses to separate
#' fields/columns.
#' @param keep_ok_only
#' A Boolean. \code{TRUE} to return only records/rows with a Penatacam status
#' of "OK".
#' @param keep_dup_rows
#' A Boolean. \code{TRUE} to return all records/rows, including duplicates.
#' @param keep_empty_cols
#' A Boolean. \code{TRUE} to return all columns, even when none of the
#' records/rows contain values.
#' @return
#' A data frame containing the file that has been read.
#'
#' @family Utilities
#'
#' @importFrom readr read_delim
#' @importFrom readr locale
#' @importFrom dplyr select_if
#' @importFrom lubridate mdy
#' @importFrom lubridate hms
#' @importFrom stringr str_to_lower
#' @importFrom stringr str_trim
#'
#' @export
read_pentacam_csv <- function(file_name,
                              delimiter = ';',
                              keep_ok_only = T,
                              keep_dup_rows = F,
                              keep_empty_cols = F) {
  csv_dat <-
    readr::read_delim(
      file_name,
      delim = delimiter,
      locale = readr::locale(encoding = "windows-1252"),
      col_names = T,
      trim_ws = T,
      show_col_types = F
    )

  # get rid of empty columns
  if (!keep_empty_cols) {
    csv_dat <- csv_dat |>
      dplyr::select_if(function(x) !(all(is.na(x)) | all(x == '')))
  }

  # get rid of duplicates
  if (!keep_dup_rows) {
    csv_dat <- csv_dat |>
      unique()
  }

  # hang onto original column names (in case there's an issue)
  orig_colnames <- colnames(csv_dat)

  # Convert file's column names to a format more suitable for variable names
  colnames(csv_dat) <- colnames(csv_dat) |>
    convert_name()

  # debug column name conversion
  #if (dim(csv_dat)[2] != length(unique(colnames(csv_dat)))) {
  #  name_matches <- colnames(csv_dat) == unique(colnames(csv_dat))
  #  first_mismatch <- min(which(name_matches == F))
  #  # identify the column names are get duplicated
  #  which(colnames(csv_dat) == colnames(csv_dat)[first_mismatch])
  #  # take a look at what's happening to the column name when it's converted
  #  convert_name(orig_colnames[first_mismatch])
  #}

  # convert date/time types
  csv_dat$d_o_birth <- lubridate::mdy(csv_dat$d_o_birth, quiet = T)
  csv_dat$exam_date <- lubridate::mdy(csv_dat$exam_date)
  csv_dat$exam_time <- lubridate::hms(csv_dat$exam_time)

  # establish factors
  csv_dat$exam_eye <- csv_dat$exam_eye |>
    stringr::str_to_lower() |>
    stringr::str_trim(side = 'right') |>
    factor(levels = c('left', 'right'))

  # calculate age at time of exam
  csv_dat$age <- ((csv_dat$exam_date - csv_dat$d_o_birth) / dyears(1)) |>
    round(0)

  if (keep_ok_only & "status" %in% colnames(csv_dat)) {
    csv_dat <- csv_dat |>
      filter(status == 'OK')
  }

  return(csv_dat)
}


##################
## convert_name
##################

#' Convert a Pentacam name to an R-friendly format
#' @description
#' Converts all characters to lower case, replaces decimals and whitespace
#' characters with an underscore ("_"), replaces symbols with text (e.g., "°"
#' with "deg"), and more.
#' @param name
#' A character string containing the file or column name to be converted.
#' @return
#' A character string containing the converted name.
#' @examples
#' convert_name('r4.0 324°')
#'
#' @family Utilities
#'
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#'
#' @export
convert_name <- function(name) {
  name |>
    stringr::str_to_lower() |>
    stringr::str_replace(':$', '') |>  # ':' at the end of column name
    stringr::str_replace_all('(?<=[:digit:])\\.(?=[:digit:])', '_') |>  # decimal point in numbers
    stringr::str_replace_all('°', 'deg') |> # degree character
    stringr::str_replace_all('-', '_') |> # dash to underscore
    stringr::str_replace_all(' ', '_') |> # space character to underscore
    stringr::str_replace_all('\\,', '_') |> # commas to underscore
    stringr::str_replace_all('\\.', '_') |> # period (abbreviation) to underscore
    stringr::str_replace_all('\\(', '_') |> # open parenthesis to underscore
    stringr::str_replace_all('\\)', '_') |> # close parenthesis to underscore
    stringr::str_replace_all('__', '_') |> # consecutive underscores to just one
    stringr::str_replace_all('__', '_') |> # consecutive underscores again (might have just created more)
    stringr::str_replace_all('â', '') |> # introduced by anonymization
    stringr::str_replace('_$', '')  # underscore at the end of column name
}
