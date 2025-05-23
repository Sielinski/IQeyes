library(IQeyes)
library(readxl)
#library(dplyr)

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


## Absolute Scale ##

# the absolute power scale that is used to determine the contour levels for
# curvature maps. The scale is not linear
absolute_scale <- c(seq(10, 30, by = 2.5),
                    seq(30, 46, by = 0.5),
                    seq(46, 50, by = 1),
                    seq(50, 90, by = 2.5)) |> unique()


## reference stats ##

# mean and standard deviation for Pentacam index values
reference_stats <- readxl::read_xlsx('data-raw/Indices.xlsx', sheet = 'Stats')

colnames(reference_stats) <- stringr::str_to_lower(colnames(reference_stats))
reference_stats$index <- stringr::str_to_lower(reference_stats$index)
reference_stats$population <- stringr::str_to_lower(reference_stats$population)

iq_reference_stats <- reference_stats |>
  dplyr::filter(iqeyes == T) |>
  dplyr::select(index, mean, sd, min, max, population, doi)

## Store datasets
usethis::use_data(join_fields, overwrite = T, internal = F)
usethis::use_data(absolute_scale, overwrite = T, internal = F)
usethis::use_data(iq_reference_stats, overwrite = T, internal = F)


#############################
## Read reference datasets
#############################

# a reference data frame, mapping the columns in COR-PWR to polar and Cartesian
# coordinates
radii_degrees <- readRDS('data-raw/radii_degrees.RDS')

# measure points from COR-PWR that are adjacent to each other
adjacent_points_dat<- readRDS('data-raw/adjacent_points_dat.RDS')

risk_factors <- read.csv('data-raw/risk_factors.CSV')

## Store datasets
usethis::use_data(radii_degrees, overwrite = T, internal = F)
usethis::use_data(adjacent_points_dat, overwrite = T, internal = F)
usethis::use_data(risk_factors, overwrite = T, internal = F)


##########################
## Read sample datasets
##########################

## Reference Shapes

# cluster IDs and shape names for the reference shapes
reference_shapes <- readRDS('data-raw/reference_shapes.RDS')

# curvature data for the reference contours
reference_contours <- readRDS('data-raw/reference_contours.RDS')

# polygons for the reference shapes
reference_polygons <- readRDS('data-raw/reference_polygons.RDS')


## Sample Exam

# A sample exam's curvature data (from COR-PWR)
sample_curvature <- readRDS('data-raw/sample_curvature.RDS')

# A sample exam's astigmatism data (from CHAMBER)
sample_astig <- readRDS('data-raw/sample_astig.RDS')

# A sample exam's characteristic contour
sample_contour <- readRDS('data-raw/sample_contour.RDS')

## Store datasets
usethis::use_data(reference_shapes, overwrite = T, internal = F)
usethis::use_data(reference_contours, overwrite = T, internal = F)
usethis::use_data(reference_polygons, overwrite = T, internal = F)

usethis::use_data(sample_curvature, overwrite = T, internal = F)
usethis::use_data(sample_astig, overwrite = T, internal = F)
usethis::use_data(sample_contour, overwrite = T, internal = F)


########################
## Read fitted models
########################

# identify the characteristic contour
char_fit <- readRDS('data-raw/characteristic_contour_fit.RDS')

# identify the characteristic contour
bowtie_fit <- readRDS('data-raw/bowtie_fit.RDS')
bowtie_roc_coords <- readRDS('data-raw/bowtie_roc_coords.RDS')
bowtie_features <- readRDS('data-raw/bowtie_features.RDS')

## Store datasets
usethis::use_data(char_fit, overwrite = T, internal = F)
usethis::use_data(bowtie_fit, overwrite = T, internal = F)
usethis::use_data(bowtie_roc_coords, overwrite = T, internal = F)
usethis::use_data(bowtie_features, overwrite = T, internal = F)


#################
## Color scales
#################

source('color_scales.R')
