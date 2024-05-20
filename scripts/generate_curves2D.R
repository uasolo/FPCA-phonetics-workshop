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
ex2D_curves <- list()
ex2D_curves[[1]] <- bind_rows(
  tibble(Category = "A",
         y1 = c(0.7, 0.7, 0.6, 0.5, 0.5, 0.7, 1.3, 1.8, 1.2, 1, 0.6),
         y2 = c(1.5, 2, 1.5, 1.5, 1.3, 1.1, 1.4, 1.8, 2, 2.1, 2.2),
         ),
  tibble(Category = "B",
         y1 = c(0.7, 0.7, 0.6, 0.5, 0.5, 0.7, 1.3, 1.8, 1.2, 1, 0.6),
         y2 = c(2, 2.5, 2, 1.5, 1.3, 1.1, 1.4, 1.8, 2, 2.1, 2.2),
         )
)

ex2D_curves[[2]] <- ex2D_curves[[1]]

ex2D_curves[[3]] <- bind_rows(
  tibble(Category = "A",
         y1 = c(3, 3.2, 3.5, 4.5, 5.5, 6, 6.3, 6.5),
         y2 = c(8.5,  8.7,  9.0,  9.6, 10.2, 10.5, 10.8, 11.0),
  ),
  tibble(Category = "B",
         y1 = c(3, 3.2, 3.5, 4.1, 4.7, 5, 5.3, 5.5),
         y2 = c(8.5, 8.7,  9.0, 10.0, 11.0, 11.5, 11.8, 12.0),
  )
)

shift_y <- list()
shift_y[[1]] <- function(y, Category, Dim, u) {
  if (Dim == 'y2') {
    y[1:3] <- y[1:3] + 0.1 * u
  }
  return (y)
}

shift_y[[2]] <- shift_y[[1]]

shift_y[[3]] <- function(y, Category, Dim, u) {
  if (Dim == 'y1') {
    y[4:8] <- y[4:8] + 0.3 *u * c(.4, .8, 1, 1, 1)
  } else {
    y[4:8] <- y[4:8] - 0.3 * (u + 0.5 * rnorm(1, 0, 1)) * c(.4, .8, 1, 1, 1)
  }
  return (y)
}

ex2D_land <- list()
ex2D_land[[1]] <- function(Category, Dim, u, role) {
  baseLand <- c(0, .2, .4, .6, .7, .8, 1) * RefT
  nLand <- length(baseLand)
  land <- baseLand
  if (role == 'target' & Dim == 'y1') {
    if (Category == 'A') {
      shift <- 0.02 * u * RefT
    } else {
      shift <- (0.02 * u + 0.05) * RefT
    }
    # truncate shift such that landmarks do not cross
    shift <- shift %>% max(-.18 * RefT) %>% min(.18 * RefT)
    land[4:6] <- land[4:6] + shift
  }
  jitter_land <- 0.02 * RefT 
  # hope they do not cross
  land[2:(nLand-1)] <- land[2:(nLand-1)] %>%
    jitter(amount = jitter_land)
  return(land)
}

ex2D_land[[2]] <- function(Category, Dim, u, role) {
  baseLand <- c(0, .2, .4, .6, .7, .8, 1) * RefT
  nLand <- length(baseLand)
  land <- baseLand
  if (role == 'target' & Dim == 'y1') {
    shift <- 0.03 * rnorm(1, 0, 1) * RefT # indep of u and Category
    # truncate shift such that landmarks do not cross
    shift <- shift %>% max(-.18 * RefT) %>% min(.18 * RefT)
    land[4:6] <- land[4:6] + shift
  }
  jitter_land <- 0.02 * RefT 
  # hope they do not cross
  land[2:(nLand-1)] <- land[2:(nLand-1)] %>%
    jitter(amount = jitter_land)
  return(land)
}

ex2D_land[[3]] <- function(Category, Dim, u, role) {
  baseLand <- c(0, .2, .4, .6, .8, 1) * RefT
  nLand <- length(baseLand)
  land <- baseLand
  jitter_land <- 0.02 * RefT 
  # hope they do not cross
  land[2:(nLand-1)] <- land[2:(nLand-1)] %>%
    jitter(amount = jitter_land)
  return(land)
}


ex <- 3
modelCurves <- ex2D_curves[[ex]] %>% 
  group_by(Category) %>% 
  mutate(t0 = seq(0, RefT, length.out = n())) %>% 
  ungroup() %>%
  relocate(t0, .after = Category) 


modelCurves %>% 
  pivot_longer(c(y1, y2), names_to = "Dim", values_to = "y") %>% 
  mutate(Dim = factor(Dim, levels = c("y2", "y1"))) %>% 
  ggplot() +
  aes(t0, y, group = interaction(Category, Dim), color = Category) +
  geom_point() + 
  geom_line() +
  facet_grid(Dim ~ .) +
  scale_color_manual(values=Category.colors) +
  mytheme

# copy NcurvesPerCategory
# apply shift_y
# add vertical jitter
# store fd smooth object
NcurvesPerCategory <- 50
jitter_y <- 0.1
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
  pivot_longer(c(y1, y2), names_to = "Dim", values_to = "y") %>% 
  mutate(Dim = factor(Dim, levels = c("y2", "y1"))) %>%
  group_by(Category, curveId, Dim) %>%
  mutate(y = shift_y[[ex]](y,
                           cur_group()$Category %>% as.character(),
                           cur_group()$Dim %>% as.character(),
                           u[1])) %>% 
  mutate(y = jitter(y, amount = jitter_y)) 

ggplot(fdCurves) +
  aes(t0, y, color = Category, group = curveId) +
  geom_line() +
  facet_grid(Dim ~ .) +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  theme(legend.position = "bottom")

fdCurves <- fdCurves %>% 
  group_by(Category, curveId, Dim) %>%
  reframe(fdObj = list(Data2fd(t0, matrix(y))),
          u = u[1]) %>% # matrix() because of a bug in Data2fd
  ungroup() 
  

# finer sampling from fd smooth curves
# apply landmark reg
curves <- fdCurves %>% #slice_head(n = 6) %>% 
  group_by(Category, curveId, Dim) %>% 
  reframe(x1 = eval.fd(one_landmarkreg_nocurves(
      inputMarks = ex2D_land[[ex]](Category %>% as.character(),
                                 Dim %>% as.character(),
                                 u,
                                 'input'),
    
      targetMarks = ex2D_land[[ex]](Category %>% as.character(),
                                  Dim %>% as.character(),
                                  u,
                                  'target')
    ), fdObj[[1]]) %>% as.numeric()
  ) %>%
  group_by(Category, curveId, Dim) %>% 
  mutate(t1 = seq(0, RefT, length.out = n())) %>% 
  ungroup()

subset_curveId <- curves %>%
  ungroup() %>% 
  distinct(curveId) %>%
  slice_sample(n = 20)

ggplot(curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(t1, x1, color = Category, group = curveId) +
  geom_line() +
  facet_grid(Dim ~ .) +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  theme(legend.position = "bottom")


# add AR1 noise and white noise
n0_sd <- 0.02
n1_sd <- 0.02
rho1 <- 0.8
curves <- curves %>% 
  group_by(Category, curveId, Dim) %>%
  mutate(n0 = rnorm(n(), 0, n0_sd),
         n1 = rnorm(n(), 0, n1_sd),
         x2 = x1 + signal::filter(filt = 1, a = c(1, -rho1), x = n1) + n0
  ) %>% 
  select(!c(n0, n1))

curves <- curves %>%
  ungroup() %>% 
  rename(time = t1, y = x2) %>%
  select(Category, curveId, Dim, time, y)



data_dir <- "data/"
saveRDS(curves, file.path(data_dir, str_c("ex2D", ex, "rds", sep = '.')))
