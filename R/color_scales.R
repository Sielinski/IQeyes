######################
## Elevation Colors
######################

elevation_colors <- matrix(c(Inf, 160, 40, 80,
                             35, 160, 40, 80,
                             30, 200, 0, 35,
                             25, 255, 40, 0,
                             20, 255, 90, 0,
                             15, 255, 130, 0,
                             10, 255, 180, 0,
                             5, 255, 220, 0,
                             0, 0, 220, 0,
                             -5, 0, 220, 255,
                             -10, 0, 180, 255,
                             -15, 0, 140, 255,
                             -20, 0, 100, 255,
                             -25, 0, 60, 255,
                             -30, 60, 0, 220,
                             -35, 130, 0, 180,
                             -Inf, 130, 0, 180),
                           byrow = T, ncol = 4)

elevation_colors <- data.frame(breaks = elevation_colors[, 1],
                               hex_color = rgb(elevation_colors[, 2] / 255,
                                               elevation_colors[, 3] / 255,
                                               elevation_colors[, 4] / 255)) |>
  dplyr::mutate(breaks = breaks - 2.5) |>
  dplyr::arrange(breaks)


#######################
## Pachymetry Colors
#######################

pachymetry_colors <- matrix(c(300, 200, 100, 100,
                              340, 160, 50, 50,
                              380, 180, 0, 0,
                              420, 220, 80, 0,
                              460, 220, 150, 0,
                              500, 220, 200, 0,
                              540, 255, 255, 0,
                              580, 90, 230, 0,
                              620, 0, 240, 252,
                              660, 0, 185, 252,
                              700, 20, 120, 240,
                              740, 40, 60, 220,
                              780, 40, 0, 140,
                              820, 100, 40, 150,
                              860, 140, 40, 180,
                              900, 180, 40, 220),
                            byrow = T, ncol = 4)

pachymetry_colors <- data.frame(breaks = pachymetry_colors[, 1],
                                hex_color = rgb(pachymetry_colors[, 2] / 255,
                                                pachymetry_colors[, 3] / 255,
                                                pachymetry_colors[, 4] / 255))

pachymetry_colors <- data.frame(hex_color = colorRampPalette(pachymetry_colors$hex_color)(61),
                                breaks = seq(300, 900, 10)) |>
  dplyr::mutate(breaks = breaks - 5) |>
  dplyr::arrange(breaks)


######################
## Curvature Colors
######################

curvature_colors <- data.frame(
  breaks = c(10, 20, 30, 32, 34, 36, 38, 40, 42, 44, 46, 50, 60, 70, 80, 90),
  #radius = c(34.5, 22.5,  10.5,  10.0,  9.6,  9.2,  8.8,  8.4,  8.0,  7.6,  7.2,  6.8,  6.2,  5.4,  4.6,  3.8),
  hex_color = c(
    '#B428DC',  # 10
    '#8C28B4',  # 20
    '#642896',  # 30
    '#28008C',  # 32
    '#283CDC',  # 34
    '#1478F0',  # 36
    '#00B9FC',  # 38
    '#00F0FC',  # 40
    '#5AE600',  # 42
    '#FFFF00',  # 44
    '#DCC800',  # 46
    '#DC9600',  # 50
    '#DC5000',  # 60
    '#B40000',  # 70
    '#A03232',  # 80
    '#C86464'   # 90
  )
)

# Not needed, but potentially useful
rgb_df <- col2rgb(curvature_colors$hex_color) |>
  t()

curvature_colors <- curvature_colors |>
  dplyr::bind_cols(rgb_df)

# the absolute power scale that is used to determine the contour levels for
# curvature maps. The scale is not linear
absolute_scale <- c(seq(10, 30, by = 2.5),
                    seq(30, 46, by = 0.5),
                    seq(46, 50, by = 1),
                    seq(50, 90, by = 2.5)) |> unique()

# should be 61
length(absolute_scale)

curvature_colors <- data.frame(hex_color = colorRampPalette(curvature_colors$hex_color)(61),
                               breaks = absolute_scale) |>
  #mutate(breaks = breaks + 5) |>
  dplyr::arrange(desc(breaks))


####################
## Make Available
####################

usethis::use_data(elevation_colors, overwrite = T, internal = F)
usethis::use_data(pachymetry_colors, overwrite = T, internal = F)
usethis::use_data(curvature_colors, overwrite = T, internal = F)
usethis::use_data(absolute_scale, overwrite = T, internal = F)
