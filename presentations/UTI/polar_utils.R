cart2polar <- function(data, origin, x = 'x', y = 'y') {
  # this function works similarly to rticulate::transform_coord()
  # as of version 2.0.1, transform_coord() is not fully correct.
  if (!all(c(x, y) %in% colnames(data))) {
    stop(str_glue("The data.frame 'data' does not have the specified x and y columns\ 
              {x} and {y}"))
  }
  if (!(is.numeric(origin))) {
    stop(tr_glue("'origin' must be a numeric vector of length 2.\
    provided: {origin}"))
  }
  if (length(origin) != 2) {
    stop(tr_glue("'origin' must be a numeric vector of length 2.\
    provided: {origin}"))
  }
  transformed_data <- data %>%
    mutate(
      angle := atan2(.data[[y]] - origin[2], .data[[x]] - origin[1]),
      radius := sqrt((.data[[x]] - origin[1]) ^ 2 + (.data[[y]] - origin[2]) ^ 2)
    ) %>% 
    mutate(angle := case_when(
      angle < 0 ~ 2 * pi + angle,
      TRUE ~ angle
    ))
  return(transformed_data) 
}

# plotting in polar coordinates
radial_plot <- function(max_radius = 80) list(
  ylim(0, max_radius),
  coord_radial(theta = "x",
               start = -0.5 * pi, end = 0.5 * pi, 
               direction = -1, expand = FALSE),
  scale_x_continuous(breaks = pi/4 * (0:4),
                     labels = expression(0, pi/4, pi/2, 3/4*pi, pi),
                     limits = c(0, pi)),
  theme(axis.title=element_blank(),
        panel.border = element_blank())
)

