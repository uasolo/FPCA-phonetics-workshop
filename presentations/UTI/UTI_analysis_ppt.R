library(tidyverse)
source("polar_utils.R")
library(ggnewscale)

theme_set(theme_light())
tongue <- readRDS("example_tongue.RDS") %>% 
  rename(x = X, y = Y)

pl <- ggplot(tongue) +
  aes(x, y) +
  geom_path() +
  geom_point(color = 'orangered') +
  ylim(18, 82) +
  geom_text(aes(x, y, label = knot), color = 'orangered', nudge_x = -4)

pl
ggsave(filename = "plots/tongue_example.png", plot = pl, width = 3.5, height = 3)  

pl <- ggplot(tongue) +
  aes(x, y) +
  geom_path() +
  ylim(18, 82) +
  geom_vline(xintercept = 40, linetype = "dashed", color = 'blue')

pl
ggsave(filename = "plots/tongue_example_not_f.png", plot = pl, width = 3.5, height = 3)  

origin <- tibble(x = 70, y = 20)
pl <- ggplot(tongue) +
  aes(x, y) +
  geom_path() +
  geom_point(color = 'orangered') +
  geom_point(data = origin, color = 'orangered') +
  geom_text(data = origin, aes(x, y), label = "Origin", color = 'orangered', nudge_x = 8) +
  geom_segment(aes(x, y, xend = origin$x, yend = origin$y),
               color = 'blue', linetype = "dashed") +
  ylim(18, 82) 

pl
ggsave(filename = "plots/tongue_example_fan.png", plot = pl, width = 3.5, height = 3)  

tongue <- tongue %>% 
  cart2polar(origin %>% as.numeric(), 'x', 'y')

pl <- ggplot(tongue) +
  aes(angle, radius) +
  geom_path() +
  geom_point(color = 'orangered') +
  geom_text(aes(label = knot), color = 'orangered',nudge_x = 0.1, nudge_y = 2) +
  scale_x_reverse()

pl
ggsave(filename = "plots/tongue_example_rad_flat.png", plot = pl, width = 3.5, height = 3)  

pl <- ggplot(tongue) +
  aes(angle, radius) +
  geom_path() +
  geom_point(color = 'orangered') +
  radial_plot(70) 

pl
ggsave(filename = "plots/tongue_example_rad.png", plot = pl, width = 3.5, height = 3)  

tongues <- readRDS("examples_tongue.RDS" )%>% 
  rename(x = X, y = Y)
tongues <- tongues %>% 
  cart2polar(origin %>% as.numeric(), 'x', 'y')

angle_range <- tongues  %>% 
  filter(knot %in% range(knot)) %>% # first and last knot
  group_by(knot) %>% 
  summarise(median_angle = median(angle)) %>% 
  pull(median_angle) %>% 
  sort()


tongue_examples_rad_flat <- ggplot(tongues) +
  aes(angle, radius, group = frame_id, color = frame_id) +
  geom_path() +
  geom_vline(xintercept = angle_range, color = 'blue', linetype = "dashed") +
  scale_x_reverse() +
  theme(legend.position = "none")

tongue_examples_rad_flat
ggsave(filename = "plots/tongue_examples_rad_flat.png", plot = tongue_examples_rad_flat, width = 3.5, height = 3)  

tongue_examples_rad <- ggplot(tongues) +
  aes(angle, radius, group = frame_id, color = frame_id) +
  geom_path() +
  geom_vline(xintercept = angle_range, color = 'blue', linetype = "dashed") +
  theme(legend.position = "none") +
  radial_plot(70) 

tongue_examples_rad
ggsave(filename = "plots/tongue_examples_rad.png", plot = tongue_examples_rad, width = 3.5, height = 3)  

## Linear angle norm

angle_grid <- seq(angle_range[1], angle_range[2], length.out=11)

curves  <- tongues  %>% 
  group_by(frame_id) %>% 
  # linear angle normalization on the angle_range interval
  mutate(angle = (angle - min(angle)) * diff(angle_range) / diff(range(angle)) + angle_range[1]) %>% 
  # linear interpolation on angle_grid (optional for GAM) 
  reframe(bind_cols(
    knot = seq_along(angle_grid), # fake knots, needed for AR1 term in GAM
    approx(angle, radius, angle_grid) %>% as_tibble()
  )) %>% 
  rename(angle = x, radius = y)

tongues_lin_rad_flat <- ggplot(curves) +
  aes(angle, radius, group = frame_id, color = frame_id) +
  geom_path() +
  scale_x_reverse() +
  theme(legend.position = "none")

tongues_lin_rad_flat
ggsave(filename = "plots/tongues_lin_rad_flat.png", plot = tongues_lin_rad_flat, width = 3.5, height = 3)  

tongues_lin_rad <- ggplot(curves) +
  aes(angle, radius, group = frame_id, color = frame_id) +
  geom_path() +
  theme(legend.position = "none") +
  radial_plot(70) 

tongues_lin_rad
ggsave(filename = "plots/tongues_lin_rad.png", plot = tongues_lin_rad, width = 3.5, height = 3)  

## Procrustean
curves  <- tongues  %>% 
  group_by(frame_id) %>% 
  # interpolation
  reframe(bind_cols(
    knot = seq_along(angle_grid), # fake knots, needed for AR1 term in GAM
    approx(angle, radius, angle_grid, rule = 2) %>% as_tibble()
  )) %>% 
  rename(angle = x,radius = y)

tongues_proc_rad_flat <- ggplot(curves) +
  aes(angle, radius, group = frame_id, color = frame_id) +
  geom_path() +
  scale_x_reverse() +
  theme(legend.position = "none")

tongues_proc_rad_flat
ggsave(filename = "plots/tongues_proc_rad_flat.png", plot = tongues_proc_rad_flat, width = 3.5, height = 3)  

tongues_proc_rad <- ggplot(curves) +
  aes(angle, radius, group = frame_id, color = frame_id) +
  geom_path() +
  theme(legend.position = "none") +
  radial_plot(70) 

tongues_proc_rad
ggsave(filename = "plots/tongues_proc_rad.png", plot = tongues_proc_rad, width = 3.5, height = 3)  

#### knots axis

pl <- tongue %>% 
  select(knot, x, y) %>% 
  pivot_longer(cols = c(x, y), names_to = "axis", values_to = "value") %>% 
  ggplot(aes(knot, value)) +
  geom_line() +
  facet_grid(axis ~ .) +
  ylab("")

pl
ggsave(filename = "plots/tongue_knots.png", plot = pl, width = 3.5, height = 4)  

#### Dynamic 

tongue_time <- readRDS("example_tongue_time.RDS") %>% 
  rename(x = X, y = Y)
tongue_time <- tongue_time %>% 
  cart2polar(origin %>% as.numeric(), 'x', 'y')

angle_range <- tongue_time  %>% 
  filter(knot %in% range(knot)) %>% # first and last knot
  group_by(knot) %>% 
  summarise(median_angle = median(angle)) %>% 
  pull(median_angle) %>% 
  sort()

angle_grid <- seq(angle_range[1], angle_range[2], length.out=11)
curves  <- tongue_time  %>% 
  group_by(frame_id, normTime) %>% 
  # interpolation
  reframe(bind_cols(
    knot = seq_along(angle_grid), 
    approx(angle, radius, angle_grid, rule = 2) %>% as_tibble()
  )) %>% 
  rename(angle = x,radius = y)

pl <- ggplot(tongue_time) +
  aes(x, y, group = normTime, color = normTime) +
  geom_path() +
  theme(legend.position = "bottom")

pl
ggsave(filename = "plots/tongue_time.png", plot = pl, width = 3.5, height = 3)  

pl <- ggplot(tongue_time %>% mutate(knot = factor(knot, levels = 1:11))) +
  aes(x, y, group = normTime) +
  geom_path(aes(color = normTime)) +
  new_scale_color() +
  geom_point(aes(x, y), color = 'orangered', inherit.aes = FALSE, show.legend=FALSE) +
  theme(legend.position = "bottom")

pl
ggsave(filename = "plots/tongue_time_knots.png", plot = pl, width = 3.5, height = 3)  


pl <- ggplot(curves) +
  aes(angle, radius, group = normTime, color = normTime) +
  geom_path() +
  radial_plot(70) +
  theme(legend.position = "bottom")
pl
ggsave(filename = "plots/tongue_time_rad.png", plot = pl, width = 3.5, height = 3) 

pl <- ggplot(curves) +
  aes(x = angle, y= normTime) +
  geom_raster(aes(fill = radius)) +
  scale_x_reverse() +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  theme(legend.position = "bottom")
pl
ggsave(filename = "plots/ongue_time_rad_2D.png", plot = pl, width = 3.5, height = 3) 

# knots 2D

pl <- tongue_time %>% 
  select(knot, normTime,  x, y) %>% 
  pivot_longer(cols = c(x, y), names_to = "axis", values_to = "value") %>% 
  ggplot(aes(knot, normTime)) +
  geom_raster(aes(fill = value)) +
  facet_grid(axis ~ .) +
  scale_fill_gradientn(colours = terrain.colors(10), name = element_blank()) +
  theme(legend.position = "bottom")

pl
ggsave(filename = "plots/tongue_time_knots_2D.png", plot = pl, width = 3.5, height = 4) 
