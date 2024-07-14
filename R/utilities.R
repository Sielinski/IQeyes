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
#' @details
#' In some CSVs, Pentacam uses uses a semicolon as the delimiter between
#' columns, regardless of a user's configuration settings in the Pentacam
#' software. This function automatically converts semicolons to the delimiter
#' character passed in as a parameter.
#'
#' @family Utilities
#'
#' @importFrom readr read_lines
#' @importFrom readr locale
#' @importFrom readr write_lines
#' @importFrom readr read_delim
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

  # Read the file into a character vector
  file_lines <- readr::read_lines(file_name,
                                  locale = readr::locale(encoding = "windows-1252"))


  # Replace semicolons with commas
  file_lines <- gsub(";",
                        delimiter,
                        file_lines)

  # Remove trailing commas
  file_lines <- gsub(paste0(delimiter, "$"),
                     "",
                     file_lines)

  # Write the cleaned lines to a temporary file
  temp_file <- tempfile()
  readr::write_lines(file_lines, temp_file)

  # Read the cleansed file
  csv_dat <-
    readr::read_delim(
      temp_file,
      delim = delimiter,
      #locale = readr::locale(encoding = "windows-1252"),
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
  if ('exam_date' %in% colnames(csv_dat)) csv_dat$exam_date <- lubridate::mdy(csv_dat$exam_date)
  if ('exam_time' %in% colnames(csv_dat)) csv_dat$exam_time <- lubridate::hms(csv_dat$exam_time)

  if ('d_o_birth' %in% colnames(csv_dat)) {
    csv_dat$d_o_birth <- lubridate::mdy(csv_dat$d_o_birth, quiet = T)
    # calculate age at time of exam
    csv_dat$age <- ((csv_dat$exam_date - csv_dat$d_o_birth) / dyears(1)) |>
      round(0)
  }

  # establish exam_eye as a factor
  if ('exam_eye' %in% colnames(csv_dat)) {
    csv_dat$exam_eye <- csv_dat$exam_eye |>
      stringr::str_to_lower() |>
      stringr::str_trim(side = 'right') |>
      factor(levels = c('left', 'right'))
  }

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


#######################
## parse_csv_to_list
#######################

#' Read a detailed Pentacam CSV into a list
#' @description
#' Reads one the Pentacam's exported CSV files into a list object. The list will
#' contain one element for each section of the file.
#'
#' @param file_name
#' A character string containing the name of the file name to read.
#' @param delimiter
#' A character identifying the delimiter that the Pentacam file uses to separate
#' fields/columns.
#' @param pachy_flag
#' A Boolean. \code{TRUE} to indicate that the file contains pachymetry data.

#' @return
#' A list containing the file that has been read.
#'
#' @details
#' The date and time fields in the \code{patient} and \code{examination}
#' sections are \emph{not} converted to lubridate date and time types.
#'
#' The section \code{cornea edge} contains a bug.
#'
#' @family Utilities
#'
#' @importFrom readr read_lines
#' @importFrom readr locale
#' @importFrom readr write_lines
#' @importFrom readr read_delim
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_to_lower
#'
#' @export
# Function to parse the CSV file and create a list object
parse_csv_to_list <- function(file_name, delimiter = ',', pachy_flag = F) {

  # Read the file into a character vector
  file_lines <- readr::read_lines(file_name,
                                  locale = readr::locale(encoding = "windows-1252"))

  # Remove trailing commas
  file_lines <- gsub(paste0(delimiter, "$"),
                     "", file_lines)

  # Write the cleansed lines to a temporary file
  temp_file <- tempfile()
  readr::write_lines(file_lines, temp_file)

  # Read the cleansed file
  csv_dat <-
    readr::read_delim(
      temp_file,
      delim = delimiter,
      col_names = F,
      na = c("", "NA", "NaN"),
      show_col_types = F
    )

  # Initialize an empty list to store the results
  result_list <- list()
  # Variables to store current section's metadata and data get set in helper
  current_section <- NULL
  in_matrix <- in_XY <- in_metadata <- F

  # this helper function will be called at the end of every section to package
  # the just completed section into the appropriate data structure, which will
  # become a list element
  helper_end_section <- function(in_metadata, in_matrix, in_XY,
                                 matrix_rows, x_values, y_values, current_meta){

    # Save previous section data, if any
    if (in_matrix) {
      # measurement map as a numeric matrix
      current_matrix <- do.call(rbind, matrix_rows)
      colnames(current_matrix) <- x_values
      rownames(current_matrix) <- y_values
      return(current_matrix)
    } else if (in_XY) {
      # X-Y metadata as a matrix
      if (length(matrix_rows) == 0)  {
        return(NULL)
      } else {
        current_matrix <- do.call(rbind, matrix_rows) |>
          matrix(ncol = 2, byrow = T, dimnames = list(NULL, c('x', 'y')))
        return(current_matrix)
      }
    } else if (in_metadata) {
      # metadata as a data frame
      return(as.data.frame(current_meta, stringsAsFactors = FALSE))
    } else {
      return(NULL)
    }
  }

  if (pachy_flag) csv_dat[1, 1] <- 'pachymetry'

  # Loop through each row to parse sections
  for (i in 1:nrow(csv_dat)) {
    # i <- 1 + i
    line <- csv_dat[i, ]

    if (all(is.na(line))) next

    # Check if the line identifies a new section
    if (grepl("^\\[.*\\]$", line[1])) {
      # Save previous section data, if any
      list_item <- helper_end_section(in_metadata, in_matrix, in_XY,
                                      matrix_rows, x_values, y_values, current_meta)

      if (!is.null(list_item)) result_list[[current_section]] <- list_item

      current_section <- line[1] |>
        stringr::str_replace('\\[', '') |>
        stringr::str_replace('\\]', '') |>
        stringr::str_replace_all(' ', '_') |>
        stringr::str_to_lower()
      in_matrix <- F
      in_XY <- F
      in_metadata <- T

      current_meta <- list()

    } else if (!in_metadata && stringr::str_detect(line[1], '[:alpha:]')) {
      # Handle matrix initialization
      # Save previous section data, if any
      list_item <- helper_end_section(in_metadata, in_matrix, in_XY,
                                      matrix_rows, x_values, y_values, current_meta)

      if (!is.null(list_item)) result_list[[current_section]] <- list_item

      current_section <- gsub(' ', '_', line[1]) |>
        tolower() |>
        paste0('_matrix')

      in_matrix <- T
      in_XY <- F
      in_metadata <- F

      matrix_rows <- list()
      x_values <- line[-1]
      y_values <- NULL

    } else if (in_metadata &&
               (line[1] == 'X' && line[2] == 'Y') || (stringr::str_detect(line[1], 'Cornea-Edge X'))) {
      # Handle the start of an X, Y section within metadata
      # Save previous section data, if any
      list_item <- helper_end_section(in_metadata, in_matrix, in_XY,
                                      matrix_rows, x_values, y_values, current_meta)

      if (!is.null(list_item)) result_list[[current_section]] <- list_item

      if (line[1] == 'X' && line[2] == 'Y') {
        current_section <- paste0(current_section, '_xy')
      } else {
        current_section <- paste0(current_section,
                                  sub(" .*", "", line[1]) |> tolower(),
                                  '_xy')
      }

      in_matrix <- F
      in_XY <- T
      in_metadata <- F

      matrix_rows <- list()

    } else if (in_matrix) {
      # Handle matrix data lines
      y_values <- c(y_values, as.numeric(line[1]))
      matrix_rows <- append(matrix_rows, list(as.numeric(line[-1])))

    } else if (in_XY) {
      # Handle matrix data lines
      matrix_rows <- append(matrix_rows, list(x = as.numeric(line[1]),
                                              y = as.numeric(line[2])))

    } else {
      # Handle typical metadata lines
      key <- line[1] |>
        unlist() |>
        stringr::str_replace('°', '_deg') |>
        stringr::str_replace('\\[mm³\\]', '') |>
        stringr::str_replace('@', 'at') |>
        stringr::str_replace_all(' ', '_') |>
        stringr::str_to_lower()
      value <- line[2]
      names(value) <- key

      current_meta <- append(current_meta, value)
    }

  }

  # Save the last section data
  list_item <- helper_end_section(in_metadata, in_matrix, in_XY,
                                  matrix_rows, x_values, y_values, current_meta)

  if (!is.null(list_item)) result_list[[current_section]] <- list_item

  return(result_list)
}

