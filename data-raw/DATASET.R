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


## Store datasets
usethis::use_data(join_fields, overwrite = T, internal = F)
usethis::use_data(radii_degrees, overwrite = T, internal = F)


#############################
## Read reference datasets
#############################

# a reference data frame, mapping the columns in COR-PWR to polar and Cartesian
# coordinates
radii_degrees <- read.csv('data-raw/radii_degrees.csv')

# measure points from COR-PWR that are adjascent to each other
adjacent_points_dat<- read.csv('data-raw/adjacent_points_dat.csv')

# curvature data for the canonical shapes
canonical_curvature <- read.csv('data-raw/canonical_shapes_curvature.csv')


## Store datasets
usethis::use_data(ring_diameters, overwrite = T, internal = T)
usethis::use_data(adjacent_points_dat, overwrite = T, internal = T)
usethis::use_data(canonical_curvature, overwrite = T, internal = F)


########################
## Sample datasets
########################

# A sample exam's curvature data (from COR-PWR)
sample_curvature <- read.csv('data-raw/plot_dat.csv')

# A sample exam's astigmatism data (from CHAMBER)
sample_astig <- read.csv('data-raw/astig_dat.csv')

# A sample exam's characteristic contour
sample_contour <- read.csv('data-raw/contours_dat.csv')

## Store datasets
usethis::use_data(sample_curvature, overwrite = T, internal = F)
usethis::use_data(sample_astig, overwrite = T, internal = F)
usethis::use_data(sample_contour, overwrite = T, internal = F)
