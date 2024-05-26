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
ex1D_curves[[5]] <- ex1D_curves[[4]]

shift_y <- list()
shift_y[[1]] <- function(y, Category, u) {return (y)}
shift_y[[2]] <- function(y, Category, u) {
  y[1:3] <- y[1:3] + u
  return (y)
}
shift_y[[3]] <- function(y, Category, u) {
  y[1:3] <- y[1:3] + u
  y[7] <- y[7] + rnorm(1, 0, 1)
  return (y)
}
shift_y[[4]] <- shift_y[[1]]
shift_y[[5]] <- shift_y[[1]]


ex1D_land <- list() 
ex1D_land[[1]] <- function(Category, u, role) {
  inputLand <- c(0, .5, .7, 1) * RefT
  if (role == 'input') return(inputLand)
  nLand <- length(inputLand)
  targetLand <- inputLand
    
  jitter_land <- 0.04 * RefT 
  # hope they do not cross
  targetLand[2:nLand] <- targetLand[2:nLand] %>%
    jitter(amount = jitter_land)
  return(targetLand)
}

ex1D_land[[2]] <- ex1D_land[[1]]
ex1D_land[[3]] <- ex1D_land[[1]]
ex1D_land[[4]] <- function(Category, u, role) {
  inputLand <- c(0, .2, .4, .6, .8 ,1) * RefT
  if (role == 'input') return(inputLand)
  nLand <- length(inputLand)
  if (Category == 'A') {
    targetLand <- c(0, .2, .4, .9, 1.1, 1.3) * RefT
  } else {
    targetLand <- inputLand
  }
  stopifnot("`inputLand` and `targetLand` have different landmark counts" =
              length(inputLand) == length(targetLand))
  shift <- 0.1 * u * RefT
  # truncate shift such that landmarks do not cross
  shift <- shift %>% max(-.18 * RefT) %>% min(.18 * RefT)
  targetLand[4:6] <- targetLand[4:6] + shift
  jitter_land <- 0.02 * RefT 
  # hope they do not cross
  targetLand[2:nLand] <- targetLand[2:nLand] %>%
    jitter(amount = jitter_land)
  return(targetLand)
}
ex1D_land[[5]] <- ex1D_land[[1]]

ex <- 4
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
jitter_y <- 0.2
fdCurves <- modelCurves %>%
  group_by(Category) %>% 
  expand_grid(id_ = 1:NcurvesPerCategory) %>% 
  group_by(Category, id_) %>% 
  mutate(curveId = cur_group_id()) %>%
  ungroup() %>% 
  arrange(curveId, t0) %>% 
  mutate(across(c(Category, curveId), ~ factor(.x))) %>% 
  group_by(Category, curveId) %>%
  mutate(u = rnorm(1, 0, 1)) %>% 
  mutate(y = shift_y[[ex]](y,
                           cur_group()$Category %>% as.character(),
                           u[1])) %>% 
  mutate(y = jitter(y, amount = jitter_y))

ggplot(fdCurves) +
  aes(t0, y, color = Category, group = curveId) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  theme(legend.position = "bottom")

fdCurves <- fdCurves %>% 
  group_by(Category, curveId) %>%
  reframe(fdObj = list(Data2fd(t0, matrix(y))),
          u = u[1]) %>% # matrix() because of a bug in Data2fd
  ungroup() 

  
# finer sampling from fd smooth curves
# apply landmark reg
curves <- fdCurves %>% #slice_head(n = 6) %>% 
  group_by(Category, curveId) %>% 
  mutate(inputMarks = list(ex1D_land[[ex]](Category %>% as.character(),
                                           u,
                                           'input')),
         targetMarks = list(ex1D_land[[ex]](Category %>% as.character(),
                                            u,
                                            'target'))
         ) %>% 
  reframe(x1 = eval.fd(one_landmarkreg_nocurves(
    inputMarks = inputMarks[[1]],
    targetMarks = targetMarks[[1]]
  ), fdObj[[1]]) %>% as.numeric(),
  t1 = seq(0, last(targetMarks[[1]]), length.out = length(x1))
  )

subset_curveId <- curves %>%
  ungroup() %>% 
  distinct(curveId) %>%
  slice_sample(n = 20)

ggplot(curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(t1, x1, color = Category, group = curveId) +
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
         x2 = x1 + signal::filter(filt = 1, a = c(1, -rho1), x = n1) + n0
  ) %>% 
  select(!c(n0, n1))


do_gaps <- FALSE
if (do_gaps) {
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
}   

curves <- curves %>%
  ungroup() %>% 
  rename(time = t1, y = x2) %>%
  select(Category, curveId, time, y)



data_dir <- "data/"
saveRDS(curves, file.path(data_dir, str_c("ex1D", ex, "rds", sep = '.')))

#### landmark reg

generate_land <- function(curves, ex) {
  if (ex == 4) {
    land <- curves %>% 
      group_by(curveId, Category) %>% 
      reframe(l1 = 0,
             l2 = jitter(0.2, amount = 0.02),
             l3 = max(time))
  } else if (ex == 5) {
    boundary <- c(A = 0.15, B = 0.25)
    land <- curves %>% 
      group_by(curveId, Category) %>% 
      reframe(l1 = 0,
             l2 = jitter(boundary[cur_group()$Category], amount = 0.02),
             l3 = max(time))  }
  return(land)
}

land <- generate_land(curves, ex)

saveRDS(land, file.path(data_dir, str_c("land1D", ex, "rds", sep = '.')))






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
