library(fda)
library(funData)
library(MFPCA)
library(tidyverse)
library(emmeans)

# install.packages("devtools")
# devtools::install_github("uasolo/landmarkregUtils")
library(landmarkregUtils)

mytheme <- theme_light() +
  theme(text = element_text(size = 16))

Category.colors <- c("darkslategray", "orangered")


plots_dir <- "presentations/plots/"
data_dir <- "data/"

ex <- 3 # change according to ex number
curves <- readRDS(file.path(data_dir, str_c("ex2D", ex, "rds", sep = '.'))) %>% ungroup() %>%
  mutate(across(c(curveId, Category), ~ factor(.x)))
nCurves <- curves %>% distinct(curveId) %>% nrow()
# For simplicity, skipping the construction of common sampling period,
# as these curves are all sampled on a common grid

# plot a few curves
subset_curveId <- curves %>%
  ungroup() %>% 
  distinct(curveId) %>%
  slice_sample(n = 20)

ggplot(curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = y, group = curveId, color = Category) +
  facet_grid(Dim ~ .) +
  # aes(x = time, y = y, group = interaction(curveId, Dim), color = Category) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  theme(legend.position = "bottom",
        axis.title.y = element_blank())



# build a multiFunData object
curvesFun2D <- lapply(c("y1", "y2"), function(y)
  long2irregFunData(curves %>% filter(Dim == {{y}}),
                    id = "curveId",
                    time = "time",
                    value = "y") %>% 
    as.funData()
) %>% 
  multiFunData()

# Compute FPCA
nPC <- 2
mfpca <- MFPCA(curvesFun2D,
               M = nPC,
               uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA"))
)

# Prop of explained var
mfpca$values  / sum( mfpca$values)

# scores st. dev.
sdFun <- mfpca$values %>% sqrt()
# PC curves to be plotted
PCcurves <- expand_grid(PC = 1:nPC,
                        Dim = 1:2, 
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, Dim, fractionOfStDev) %>%
  reframe(
    funData2long1(
      mfpca$meanFunction[[Dim]] +
        fractionOfStDev * sdFun[PC] * mfpca$functions[[Dim]][PC],
      time= "time", value = "y")
  ) %>% 
  mutate(Dim = factor(Dim, levels = c(2,1), labels = c('y2', 'y1')))
# Plot
ggplot(PCcurves) +
  aes(x = time, y = y, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered",
                        breaks = c(-1, 0 , 1)) +
  facet_grid(Dim ~ PC,
             scales = "free_y",
             labeller = labeller(PC = ~ str_glue("PC{.x}"))) + #,
                                 # Dim = Dimlabels)) +
  labs(color = expression(frac(s[k], sigma[k]))) +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0), color = 'black', linewidth = 1.2) +
  mytheme +
  theme(legend.position = "right",
        axis.title.y = element_blank())


# collect PC scores
PCscores <- mfpca$scores %>%
  `colnames<-`( paste0("s", 1:nPC)) %>%
  as_tibble() %>%
  bind_cols(curves %>% distinct(curveId, Category), .)

# scatterplot PC scores s1 and s2 by Category
ggplot(PCscores) +
  aes(x = s1, y = s2, color = Category) +
  geom_point() +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "right")

# boxplots PC scores by Category
PCscores %>% 
  pivot_longer(cols = s1:all_of(str_glue("s{nPC}")), 
               names_to = "PC", values_to = "score") %>% 
  ggplot() +
  aes(x = Category, y = score, color = Category) +
  geom_boxplot() +
  facet_wrap(~ PC) +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "bottom")

s <- 1 # score index
model_eq <- str_glue("s{s} ~ Category") %>% as.formula()
mod <- lm(model_eq, data = PCscores)
mod %>% summary()
emmeans(mod, pairwise ~ Category)

# model predictions with error bars
emmeans(mod, pairwise ~ Category)$emmeans %>%
  as_tibble() %>% 
  ggplot() +
  aes(x = Category, y = emmean, color = Category) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width = .2) +
  scale_color_manual(values=Category.colors) +
  ylab(str_glue("s{s}")) +
  ggtitle(str_glue("Predicted scores according to regr model: s{s} ~ Category")) +
  mytheme 


# reconstruct predicted curves
predCurves <- emmeans(mod, pairwise ~ Category)$emmeans %>%
  as_tibble() %>%
  expand_grid(Dim = 1:2) %>% 
  group_by(Category, Dim) %>% 
  reframe(bind_cols(
    funData2long1(mfpca$meanFunction[[Dim]] +
                    emmean * mfpca$functions[[Dim]][s], value = "y"),
    funData2long1(mfpca$meanFunction[[Dim]] +
                    lower.CL * mfpca$functions[[Dim]][s], value = "yl") %>% 
      select(yl),
    funData2long1(mfpca$meanFunction[[Dim]] +
                    upper.CL * mfpca$functions[[Dim]][s], value = "yu") %>% 
      select(yu)
  ))

predCurves %>% 
  mutate(Dim = factor(Dim, levels = c(2, 1), labels = c('y2', 'y1'))) %>% 
  ggplot() +
  aes(time, y, color = Category) +
  geom_line() +
  geom_ribbon(aes(x = time, ymin = yl, ymax = yu, fill = Category),
              alpha = 0.3, inherit.aes = FALSE) +
  facet_grid(Dim ~ .) +
  scale_color_manual(values=Category.colors) +
  scale_fill_manual(values=Category.colors) +
  ggtitle(str_glue("Reconstructed curves according to regr model: s{s} ~ Category")) +
  mytheme +
  theme(legend.position = "bottom",
        axis.title.y = element_blank())


# Trajectory representation (mock formants for ex 3)
ggplot(predCurves %>% select(-c(yu, yl)) %>% 
               pivot_wider(names_from = "Dim", values_from = "y", names_prefix = "y")) +
  aes(x = y2, y = y1, group = Category, color = Category) +
  geom_path(linewidth = 1, arrow = arrow(length = unit(0.5, "cm"), type = "closed")) +
  scale_color_manual(values=Category.colors) +
  scale_x_reverse() + scale_y_reverse() +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "bottom")
