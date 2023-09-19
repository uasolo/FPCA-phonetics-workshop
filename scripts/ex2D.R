library(fda)
library(funData)
library(MFPCA)
library(tidyverse)
library(landmarkregUtils)
library(emmeans)

mytheme <- theme_light() +
  theme(text = element_text(size = 16))

Category.colors <- c("slategray4", "forestgreen")


plots_dir <- "presentations/plots/"
data_dir <- "data/"

ex <- 1 # change according to ex number
curves <- read_csv(file.path(data_dir, paste("ex2D", ex, "csv", sep = '.'))) %>% 
  mutate(across(c(curveId, Category), ~ factor(.x)))
nCurves <- curves %>% distinct(curveId) %>% nrow()




# plot a few curves
pl <- ggplot(curves %>% filter(curveId %in% sample(nCurves, 20)) %>%
               pivot_longer(starts_with("y"), names_to = "DIM", values_to = "y")) +
  aes(x = time, y = y, group = curveId, color = Category) +
  facet_wrap(~ DIM, ncol = 1) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("ex2D", ex, "curves", 'png', sep = '.')), pl,
       width = 1800, height = 1600, units = "px"
)

# build a multiFunData object
curvesFun2D <- lapply(c("y1", "y2"), function(y)
  long2irregFunData(curves, id = "curveId", time = "time", value = {{y}}) %>% 
    as.funData()
) %>% 
  multiFunData()

# Compute FPCA
mfpca <- MFPCA(curvesFun2D,
               M = 2,
               uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA"))
)

# Prop of explained var
mfpca$values  / sum( mfpca$values)

# scores st. dev.
sdFun <- mfpca$scores %>% apply(2, sd) 
# PC curves to be plotted
PCcurves <- expand_grid(PC = 1:2,
                        DIM = 1:2, 
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, DIM, fractionOfStDev) %>%
  reframe(
    funData2long1(
      mfpca$meanFunction[[DIM]] +
        fractionOfStDev * sdFun[PC] * mfpca$functions[[DIM]][PC],
      time= "time", value = "y")
  )
# Plot
DIMlabels <- c(`1` = "y1", `2` = "y2")
pl <- ggplot(PCcurves) +
  aes(x = time, y = y, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  facet_grid(DIM ~ PC,
             scales = "free_y",
             labeller = labeller(PC = ~ str_glue("PC{.x}"),
                                 DIM = DIMlabels)) +
  labs(color = expression(frac(s[k], sigma[k]))) +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0), color = 'black', linewidth = 1.2) +
  xlab("registered time") +
  mytheme +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("ex2D", ex, "FPCA_curves", 'png', sep = '.')), pl,
       width = 2000, height = 2000, units = "px"
)

# collect PC scores
PCscores <- mfpca$scores %>%
  `colnames<-`( paste0("s", 1:2)) %>%
  as_tibble() %>%
  bind_cols(curves %>% distinct(curveId, Category), .)

# scatterplot PC scores s1 and s2 by Category
pl <- ggplot(PCscores) +
  aes(x = s1, y = s2, color = Category) +
  geom_point() +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "right")

ggsave(file.path(plots_dir, str_c("ex2D", ex, "PCscores_scatter", 'png', sep = '.')), pl,
       width = 1800, height = 1200, units = "px"
)


# model s1
mod <- lm(s1 ~ Category, data = PCscores)
emmeans(mod, pairwise ~ Category)

# reconstruct predicted curves
emm <- emmeans(mod, pairwise ~ Category)$emmeans %>%
  as_tibble() %>% 
  select(Category, emmean) %>% 
  rename(s1= emmean)

predCurves <- expand_grid(emm, DIM = 1:2) %>% 
  group_by(Category, DIM) %>% 
  reframe(funData2long(mfpca$meanFunction[[DIM]] +
                         s1 * mfpca$functions[[DIM]][1])) %>%
  select(!id) %>% 
  rename(y = value)

pl <- ggplot(predCurves) +
  aes(time, y, color = Category) +
  geom_line() +
  facet_wrap(~ DIM, ncol = 1, labeller = labeller(DIM = DIMlabels)) +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("ex2D", ex, "pred_curves", 'png', sep = '.')), pl,
       width = 1500, height = 1200, units = "px"
)

# Trajectory representation
pl <- ggplot(predCurves %>%
               pivot_wider(names_from = "DIM", values_from = "y", names_prefix = "y")) +
  aes(x = y1, y = y2, group = Category, color = Category) +
  geom_path(linewidth = 1, arrow = arrow(length = unit(0.5, "cm"), type = "closed")) +
  scale_color_manual(values=Category.colors) +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("ex2D", ex, "pred_traj", 'png', sep = '.')), pl,
       width = 1200, height = 1200, units = "px"
)
