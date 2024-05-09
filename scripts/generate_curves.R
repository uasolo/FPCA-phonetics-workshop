library(fda)
library(funData)
library(tidyverse)
library(landmarkregUtils)
library(signal)

select <- dplyr::select # conflict with MASS
filter <- dplyr::filter # conflict with signal

mytheme <- theme_light() +
  theme(text = element_text(size = 16))

Category.colors <- c("slategray4", "orangered")


# All curves sampled on t0, one curve per category, 

RefT <- 0.5 # reference total duration
ex1D_curves <- list()
ex1D_curves[[1]] <- bind_rows(
  tibble(Category = "A",
         y = c(2.5, 3, 2.5, 0.5, -0.3, 0.5, 0, 0, 0, -1.5, -2.5)),
  tibble(Category = "B",
         y = c(2.5, 3, 2.5, 0.5, -0.3, 0.5, 2, 0, -.5, -1.5, -2.5))
)
ex1D_curves[[2]] <- ex1D_curves[[1]]
ex1D_curves[[3]] <- bind_rows(
  tibble(Category = "A",
         y = c(2.5, 3, 2.5, 0.5, -0.3, 0.5, 0, 0, 0, -1.5, -2.5)),
  tibble(Category = "B",
         y = c(0.5, 1, 0.5, 0.5, -0.3, 0.5, 2, 0, -.5, -1.5, -2.5))
)
ex1D_curves[[4]] <- bind_rows(
  tibble(Category = "A",
         y = c(2.5, 3, 2.5, 0.5, -0.3, 0.5, 2, 0, -.5, -1.5, -2.5)),
  tibble(Category = "B",
         y = c(2.5, 3, 2.5, 0.5, -0.3, 0.5, 2, 0, -.5, -1.5, -2.5))
)

shift_y <- list()
shift_y[[1]] <- list('A' = \(y) return (y),
                     'B' = \(y) return (y)
                     )
shift_y[[2]] <- list('A' = \(y) {
  y[1:3] <- y[1:3] + rnorm(1, 0, 1)
  return (y)
})
shift_y[[2]][['B']] <- shift_y[[2]][['A']]
shift_y[[3]] <- list('A' = \(y) {
  y[1:3] <- y[1:3] + rnorm(1, 0, 1)
  y[7] <- y[7] + rnorm(1, 0, 1)
  return (y)
})
shift_y[[3]][['B']] <- shift_y[[3]][['A']]
shift_y[[4]] <- shift_y[[1]]

# skip first landmark that is always zero
baseLand <- c(.5, .7, 1) * RefT
ex1D_land <- list() 
ex1D_land[[1]] <- list(A = list(input = baseLand,
                                target = baseLand),
                       B = list(input = baseLand,
                                target = baseLand)
)
ex1D_land[[2]] <- ex1D_land[[1]]
ex1D_land[[3]] <- ex1D_land[[1]]
ex1D_land[[4]] <- list(A = ex1D_land[[1]][['A']],
  B = list(input = c(.2, .4, .6, .8 ,1) * RefT,
                                target = c(.2, .4, .75, .9, 1.15) * RefT
))

ex <- 1
modelCurves <- ex1D_curves[[ex]] %>% 
  group_by(Category) %>% 
  mutate(t0 = seq(0, RefT, length.out = n())) %>% 
  ungroup() %>%
  relocate(t0, .after = Category)
  


ggplot(modelCurves) +
  aes(t0, y, group = Category, color = Category) +
  geom_point() + 
  geom_line() +
  scale_color_manual(values=Category.colors) +
  mytheme

# copy NcurvesPerCategory
# apply shift_y
# add vertical jitter
# store fd smooth object
NcurvesPerCategory <- 50
jitter_y <- 0.4
fdCurves <- modelCurves %>%
  group_by(Category) %>% 
  expand_grid(id_ = 1:NcurvesPerCategory) %>% 
  group_by(Category, id_) %>% 
  mutate(curveId = cur_group_id()) %>%
  ungroup() %>% 
  arrange(curveId, t0) %>% 
  group_by(Category, curveId) %>%
  mutate(y = shift_y[[ex]][[cur_group()$Category]](y)) %>% 
  mutate(y = jitter(y, amount = jitter_y)) %>% 
  reframe(fdObj = list(Data2fd(t0, matrix(y)))) %>% # matrix() because of a bug in Data2fd
  ungroup() %>% 
  mutate(across(c(Category, curveId), ~ factor(.x))) 
  
# plot a smooth curve from fdObj           
fdCurves %>%
  filter(curveId == 77) %>% 
  pull(fdObj) %>% 
  `[[`(1) %>% 
  plot()
  
# finer sampling from fd smooth curves
t1 <- seq(0, RefT, by = .01)
curves <- fdCurves %>% 
  group_by(Category, curveId) %>% 
  reframe(t1 = t1,
          y1 = eval.fd(t1, fdObj[[1]]) %>% as.numeric()
          ) %>% 
  ungroup()

subset_curveId <- curves %>%
  ungroup() %>% 
  distinct(curveId) %>%
  slice_sample(n = 20)

# ggplot(curves) +
ggplot(curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(t2, y2, color = Category, group = curveId) +
  geom_line() +
  # geom_point() +
  scale_color_manual(values=Category.colors) +
  # facet_wrap(~ curveId, nrow = 4) +
  mytheme  +
  theme(legend.position = "bottom")




# add AR1 noise and white noise
n0_sd <- 0.05
n1_sd <- 0.05
rho1 <- 0.8
curves <- curves %>% 
  group_by(Category, curveId) %>%
  mutate(n0 = rnorm(n(), 0, n0_sd),
         n1 = rnorm(n(), 0, n1_sd),
         y2 = y1 + signal::filter(filt = 1, a = c(1, -rho1), x = n1) + n0
  ) %>% 
  select(!c(n0, n1))

# remove random segments
P_gap <- 0.5
curves <- curves %>% 
  group_by(Category, curveId) %>%
  mutate(gap = runif(1) < P_gap) %>% 
  filter(t1 < RefT * 0.1 |
           # (
             t1 > RefT * runif(1, 0.3, 0.4)
             # & t1 < RefT * 0.6
             # )
         |
           # t1 > RefT * runif(1, 0.8, 0.9) |
           !gap) %>% 
  select(!gap)

# produce category/curve specific landmarks, move time around.
# modelLand <- ex1D_land[[ex]] %>% 
#   mutate(across(starts_with("l"), ~ RefT * .x))

jitter_land <- 0.04 * RefT
curves <- curves %>% 
  group_by(Category, curveId) %>%
  mutate(t2 = landmarkreg_timeSamples(t1,
                                     c(0, ex1D_land[[ex]][[cur_group()$Category]][['input']]),
                                     c(0, ex1D_land[[ex]][[cur_group()$Category]][['target']] %>% 
                                       jitter(amount = jitter_land)
                                     ))
  )
           
           
           
    #        {
    # inputMarks <- cur_group() %>%
    #   inner_join(modelLand, by = "Category") %>% 
    #   select(starts_with("l")) %>% 
    #   as.numeric()
    # targetMarks <- c(0, jitter(inputMarks[-1], amount = jitter_land))
    # landmarkreg_timeSamples(t1, inputMarks, targetMarks)
    # }) 
  # mutate(t2 = case_when(
  #   t2 == min(t2) ~ 0,
  #   TRUE ~ t2
  # ))

curves <- curves %>%
  ungroup() %>% 
  rename(time = t2, y = y2) %>%
  select(Category, curveId, time, y)



data_dir <- "data/"
saveRDS(curves, file.path(data_dir, str_c("ex1D", ex, "rds", sep = '.')))

#### interp

resample_curves <- function(curves, sp) {
  res <- curves %>% 
      group_by(curveId) %>% 
      reframe(t2 = seq(0, max(time), by = sp),
              y2 = approx(time, y, t2)$y
      )
  return(res %>% 
           ungroup() %>% 
           select(curveId, t2, y2) %>% 
           rename(time = t2, y = y2)
           )
}

rc <- resample_curves(curves, 0.01)

rc <- rc %>% inner_join(curves %>% distinct(curveId, Category))
