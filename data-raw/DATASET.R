library(IQeyes)
library(readxl)

###############################
## Define reference datasets
###############################

## join_fields ##

# establish the fields that are common across files
join_fields <-
  c('last_name',
    'first_name',
    'pat_id',
    'exam_date',
    'exam_time',
    'exam_eye')


## ring_diameters ##

# To find the steepest and flattest meridians on the 3, 5, and 7 mm rings around
# the apex of cornea, create a data frame (ring_diameters) with radial positions
# for every integer angle from 0 - 359° on each of the target rings
ring_diameters <- data.frame(ring_diam = 3, angle = seq(0, 359, 1)) |>
  dplyr::bind_rows(data.frame(ring_diam = 5, angle = seq(0, 359, 1))) |>
  dplyr::bind_rows(data.frame(ring_diam = 7, angle = seq(0, 359, 1)))

# Add Cartesian coordinates to ring_diameters
ring_diameters <- ring_diameters |>
  dplyr::bind_cols(polar_to_cartesian(ring_diameters$ring_diam / 2, ring_diameters$angle)) |>
  dplyr::mutate(ring_diam = round(ring_diam, 5), angle = round(angle, 5))


## reference stats ##

reference_stats <- readxl::read_xlsx('data-raw/Indices.xlsx', sheet = 'Stats')

colnames(reference_stats) <- stringr::str_to_lower(colnames(reference_stats))
reference_stats$index <- stringr::str_to_lower(reference_stats$index)
reference_stats$population <- stringr::str_to_lower(reference_stats$population)

iq_reference_stats <- reference_stats |>
  dplyr::filter(iqeyes == T) |>
  dplyr::select(index, mean, sd, min, max, population, doi)

## Store datasets
usethis::use_data(join_fields, overwrite = T, internal = F)
usethis::use_data(ring_diameters, overwrite = T, internal = F)
usethis::use_data(iq_reference_stats, overwrite = T, internal = F)


#############################
## Read reference datasets
#############################

# a reference data frame, mapping the columns in COR-PWR to polar and Cartesian
# coordinates
radii_degrees <- read.csv('data-raw/radii_degrees.csv')

# measure points from COR-PWR that are adjacent to each other
adjacent_points_dat<- read.csv('data-raw/adjacent_points_dat.csv')

## Store datasets
usethis::use_data(radii_degrees, overwrite = T, internal = F)
usethis::use_data(adjacent_points_dat, overwrite = T, internal = F)


########################
## Sample datasets
########################

# function to convert date, times, and factors to their appropriate format
read_package_data <- function(file_name) {

  #file_name <- 'data-raw/plot_dat.csv'
  dat <- read.csv(file_name)

  dat$exam_date <- lubridate::ymd(dat$exam_date)
  dat$exam_time <- lubridate::hms(dat$exam_time)

  if ('d_o_birth' %in% colnames(dat)) dat$d_o_birth <- lubridate::ymd(dat$d_o_birth)

  # establish factors
  dat$exam_eye <- dat$exam_eye |>
    stringr::str_to_lower() |>
    stringr::str_trim(side = 'right') |>
    factor(levels = c('left', 'right'))

  return(dat)
}

# curvature data for the canonical shapes
canonical_curvature <- read_package_data('data-raw/canonical_shapes_curvature.csv')

# A sample exam's curvature data (from COR-PWR)
sample_curvature <- read_package_data('data-raw/plot_dat.csv')

# A sample exam's astigmatism data (from CHAMBER)
sample_astig <- read_package_data('data-raw/astig_dat.csv')

# A sample exam's characteristic contour
sample_contour <- read_package_data('data-raw/contours_dat.csv')

## Store datasets
usethis::use_data(canonical_curvature, overwrite = T, internal = F)
usethis::use_data(sample_curvature, overwrite = T, internal = F)
usethis::use_data(sample_astig, overwrite = T, internal = F)
usethis::use_data(sample_contour, overwrite = T, internal = F)
