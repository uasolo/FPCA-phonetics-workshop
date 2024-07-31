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

Category.colors <- c(A = "darkslategray", B = "orangered")


plots_dir <- "presentations/plots"
data_dir <- "data"

ex <- 1 # change according to ex number
raw_curves <- readRDS(file.path(data_dir, str_c("ex1D", ex, "rds", sep = '.'))) %>% ungroup() %>% 
  mutate(across(c(curveId, Category), ~ factor(.x)))


# Create a common sampling period (sp) 
# important to reduce complexity of FPCA computation
sp <- 0.01 
maxT <- max(raw_curves$time)
grid <- seq(0, maxT, by = sp) # unified sampling grid 
curves <- raw_curves %>% 
  group_by(curveId, Category) %>% # all the factors at the level of curveId or higher (e.g. speaker)
  reframe(approx(time, y, grid) %>% as_tibble()) %>% # linear interpolation on grid 
  ungroup() %>% 
  rename(time = x, y = y)



# plot a few curves
subset_curveId <- raw_curves %>%
  ungroup() %>% 
  distinct(curveId) %>%
  slice_sample(n = 20)

# ylim <- c(-3.8, 4)
ggplot(curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = y, group = curveId, color = Category) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  # ylim(ylim) +
  mytheme  +
  theme(legend.position = "bottom")

# build a funData object
curvesFun <- long2irregFunData(curves, id = "curveId", time = "time", value = "y") %>% 
    as.funData()

# Compute FPCA
fpca <- PACE(curvesFun)

# how many PCs?
fpca$npc

# var of PC scores 
fpca$values # eigenvalues of covariance operator
fpca$scores %>% apply(2, var) # same 

# Prop of explained var
round(fpca$values  / sum( fpca$values) , digits = 3)

# scores st. dev.
sdFun <- fpca$values %>% sqrt()
# PC curves to be plotted
nPC <- 2
PCcurves <- expand_grid(PC = 1:nPC,
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
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered",
                        breaks = c(-1, 0 , 1)) +
  facet_wrap(~ PC, nrow = 1,
             labeller = labeller(PC = ~ str_glue("PC{.x}"))) +
  labs(color = expression(frac(s[k], sigma[k]))) +
  xlab("registered time") +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0), color = 'black', linewidth = 1.2) +
  mytheme +
  theme(legend.position = "bottom")

# collect PC scores
PCscores <- fpca$scores %>%
  `colnames<-`( paste0("s", 1:fpca$npc)) %>%
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
  pivot_longer(cols = s1:all_of(str_glue("s{fpca$npc}")), 
               names_to = "PC", values_to = "score") %>% 
  filter(PC %in% str_c("s", 1:nPC)) %>% 
  ggplot() +
  aes(x = Category, y = score, color = Category) +
  geom_boxplot() +
  facet_wrap(~ PC) +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "bottom")
  

    
# model scores with linear regression

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
  group_by(Category) %>% 
  reframe(bind_cols(
    funData2long1(fpca$mu + emmean * fpca$functions[s], value = "y"),
    funData2long1(fpca$mu + lower.CL * fpca$functions[s], value = "yl") %>% 
      select(yl),
    funData2long1(fpca$mu + upper.CL * fpca$functions[s], value = "yu") %>% 
      select(yu)
    ))
    
 
ggplot(predCurves) +
  aes(time, y, color = Category) +
  geom_line() +
  geom_ribbon(aes(x = time, ymin = yl, ymax = yu, fill = Category),
              alpha = 0.3, inherit.aes = FALSE) +
  scale_color_manual(values=Category.colors) +
  scale_fill_manual(values=Category.colors) +
  # ggtitle(str_glue("Reconstructed curves according to regr model: s{s} ~ Category")) +
  mytheme +
  # xlab("Registered time") +
  theme(legend.position = "bottom")

# diff curve B - A 
diffCurve <- emmeans(mod, pairwise ~ Category)$contrasts %>% 
  as_tibble() %>% 
  reframe(bind_cols(
# emmeans computes A - B, so we need to invert the sign of all curves
    funData2long1(-estimate * fpca$functions[s], value = "y"),
    funData2long1(-(estimate - 1.96 * SE) * fpca$functions[s], value = "yl") %>% 
      select(yl),
    funData2long1(-(estimate + 1.96 * SE) * fpca$functions[s], value = "yu") %>% 
      select(yu)
  ))

ggplot(diffCurve) +
  aes(time, y) +
  geom_line() +
  geom_ribbon(aes(x = time, ymin = yl, ymax = yu),
              alpha = 0.3, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'red') +
  ggtitle(str_glue("Difference B - A according to regr model: s{s} ~ Category")) +
  mytheme 

# GAM
# basic form (no AR1, no curve-specific random smooths)
GAM <- bam(y ~ Category + s(time, by = Category), data = curves)
summary(GAM)
plot_smooth(GAM, view = "time", plot_all = "Category",
            print.summary = FALSE, col = Category.colors)
plot_diff(GAM, view = "time", comp = list(Category = c("B", "A")))
