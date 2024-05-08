library(fda)
library(funData)
library(MFPCA)
library(tidyverse)
library(emmeans)
library(mgcv)
library(itsadug)

# install.packages("devtools")
# devtools::install_github("uasolo/landmarkregUtils")
library(landmarkregUtils)

mytheme <- theme_light() +
  theme(text = element_text(size = 16))

Category.colors <- c(A = "slategray4", B = "orangered")


plots_dir <- "presentations/plots/"
data_dir <- "data/"

ex <- 1 # change according to ex number
raw_curves <- readRDS(file.path(data_dir, str_c("ex1D", ex, "rds", sep = '.'))) %>% ungroup() %>% 
# curves <- read_csv(file.path(data_dir, paste("ex1D", ex, "csv", sep = '.'))) %>% 
  mutate(across(c(curveId, Category), ~ factor(.x)))
# nCurves <- curves %>% distinct(curveId) %>% nrow()

sp <- 0.01 # unified sampling period
maxT <- max(raw_curves$time)
grid <- seq(0, maxT, by = sp) # unified sampling grid 
curves <- raw_curves %>% 
  group_by(curveId, Category) %>% # all the factors at the level of curveId or higher (e.g. speaker)
  reframe(approx(time, y, grid) %>% as_tibble()) %>% # linear interpolation on grid 
  ungroup() %>% 
  rename(time = x, y = y)

maxFixGap <- 4 * sp
bigGaps <- raw_curves %>% 
  group_by(curveId, Category) %>%
  mutate(time_to = lead(time),
         bigGap = time_to - time > maxFixGap) %>%
  filter(bigGap) %>% 
  select(!c(y, bigGap)) %>% 
  rename(time_from = time) %>%
  reframe(cond = cur_data() %>%
            pmap(\(time_from, time_to) str_c("(time < ", time_from, " | time > ", time_to, ")")) %>% 
            str_c(collapse = " & ")
  ) %>% 
  ungroup()


  
library(rlang)

curves <- curves %>%
  # filter(curveId %in% 1:10) %>% 
  nest(curve = c(time, y)) %>% 
  left_join(bigGaps, by = c("Category", "curveId")) %>% 
  group_by(curveId, Category) %>%
  mutate(curve = case_when(
    !is.na(cond) ~ str_c("curve[[1]] %>% filter(", cond[[1]], ") %>% list()") %>% 
      parse_expr() %>% eval(),
    TRUE ~ curve
  )) %>% 
  select(!cond) %>% 
  unnest(curve) %>% 
  ungroup()

# plot a few curves
subset_curveId <- raw_curves %>%
  ungroup() %>% 
  distinct(curveId) %>%
  slice_sample(n = 20)

ggplot(curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = y, group = curveId, color = Category) +
  # geom_line() +
  geom_point() +
  facet_wrap(~ curveId) +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  theme(legend.position = "bottom")

# build a funData object
curvesFun <- long2irregFunData(curves, id = "curveId", time = "time", value = "y") %>% 
    as.funData()

# Compute FPCA
fpca <- PACE(curvesFun)

# Prop of explained var
fpca$values  / sum( fpca$values)

# scores st. dev.
sdFun <- fpca$scores %>% apply(2, sd) 
# PC curves to be plotted
PCcurves <- expand_grid(PC = 1:2,
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, fractionOfStDev) %>%
  reframe(
    funData2long1(fpca$mu + fractionOfStDev * sdFun[PC] * fpca$functions[PC],
                  time = "time", value = "y")
  )
# Plot
ggplot(PCcurves) +
  aes(x = time, y = y, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  facet_wrap(~ PC, nrow = 1,
             labeller = labeller(PC = ~ str_glue("PC{.x}"))) +
  # labs(color = expression(frac(s[k], sigma[k]))) +
  labs(color = "Norm. scores") +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0), color = 'black', linewidth = 1.2) +
  mytheme +
  theme(legend.position = "bottom")

# collect PC scores
PCscores <- fpca$scores %>%
  `colnames<-`( paste0("s", 1:2)) %>%
  as_tibble() %>%
  bind_cols(curves %>% distinct(curveId, Category), .)

# scatterplot PC scores s1 and s2 by Category
ggplot(PCscores) +
  aes(x = s1, y = s2, color = Category) +
  geom_point() +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "right")


# model s1 or s2
mod <- lm(s1 ~ Category, data = PCscores)
emmeans(mod, pairwise ~ Category)

# reconstruct predicted curves
emm <- emmeans(mod, pairwise ~ Category)$emmeans %>%
  as_tibble() %>% 
  select(Category, emmean) %>% 
  rename(s1= emmean)

predCurves <- emm %>% 
  group_by(Category) %>% 
  reframe(funData2long1(fpca$mu + s1 * fpca$functions[1])) %>%
  rename(y = value)

ggplot(predCurves) +
  aes(time, y, color = Category) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "bottom")
