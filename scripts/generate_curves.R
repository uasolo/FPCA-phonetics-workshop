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

MeanT <- 0.5
t0 <- seq(0, MeanT, length.out = 11)
modelCurves <- bind_rows(
  tibble(Category = "A",
         t0 = t0,
         y = c(2.5, 3, 2.5, 0.5, -0.3, 0.5, 0, 0, 0, -1.5, -2.5)),
  tibble(Category = "B",
         t0 = t0,
         y = c(2.5, 3, 2.5, 0.5, -0.3, 0.5, 2, 0, -.5, -1.5, -2.5))
)

# copy NcurvesPerCategory
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
  mutate(y = jitter(y, amount = jitter_y)) %>% 
  reframe(fdObj = list(Data2fd(t0, matrix(y)))) %>% 
  mutate(across(c(Category, curveId), ~ factor(.x))) 
  
# plot a smooth curve from fdObj           
fdCurves %>%
  filter(curveId == 7) %>% 
  pull(fdObj) %>% 
  `[[`(1) %>% 
  plot()
  
# finer sampling from fd smooth curves
t1 <- seq(0, MeanT, by = .01)
curves <- fdCurves %>% 
  group_by(Category, curveId) %>% 
  reframe(t1 = t1,
          y1 = eval.fd(t1, fdObj[[1]]) %>% as.numeric()
          )

ggplot(curves) +
  aes(t2, y2, color = Category, group = curveId) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
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

# produce category/curve specific landmarks, move time around.
modelLand <- tribble(
  ~Category, ~l1, ~l2, ~l3, ~l4,
  "A", 0, .5, .7, 1,
  "B", 0, .5, .7, 1,
) %>% 
  mutate(across(starts_with("l"), ~ MeanT * .x))

jitter_land <- 0.04 * MeanT
curves <- curves %>% 
  group_by(Category, curveId) %>%
  mutate(t2 = {
    inputMarks <- cur_group() %>%
      inner_join(modelLand, by = "Category") %>% 
      select(starts_with("l")) %>% 
      as.numeric()
    targetMarks <- c(0, jitter(inputMarks[-1], amount = jitter_land))
    landmarkreg_timeSamples(t1, inputMarks, targetMarks)
    }) %>% 
  mutate(t2 = case_when(
    t2 == min(t2) ~ 0,
    TRUE ~ t2
  ))

curves <- curves %>%
  ungroup() %>% 
  rename(time = t2, y = y2) %>% 
  select(Category, curveId, time, y)

# remove random segments
curves <- curves %>% 
  group_by(Category, curveId) %>%
  {if (TRUE) filter(., time < 0.05 | time > runif(1, 0.15, 0.2)) else .} %>% 
  {if (TRUE) filter(., time < 0.38 | time > 0.44) else .}


data_dir <- "data/"
ex <- 1
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
