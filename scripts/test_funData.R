library(fda)
library(funData)
library(MFPCA)
library(tidyverse)
library(emmeans)
library(mgcv)
library(itsadug)

ex <- 1 # change according to ex number
curves <- read_csv(file.path("../data/", paste("ex1D", ex, "csv", sep = '.'))) %>% 
  mutate(across(c(curveId, Category), ~ factor(.x)))

nCurves <- curves %>% select(curveId) %>% n_distinct()

ggplot(curves %>% filter(curveId %in% sample(nCurves, 10))) +
  aes(time, y, color = Category) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  theme_light() +
  theme(text = element_text(size = 15),
        legend.position = "bottom")

curvesFun <- long2irregFunData(curves, "curveId", "time", "y")
curvesFun
plot(curvesFun)
curvesFun <- as.funData(curvesFun)
curvesFun@argvals %>% str()
curvesFun@X %>% dim()

curvesFun %>% nObs()
curvesFun %>% nObsPoints()

curvesFun[1] %>% approxNA() %>% funData2fd() %>% eval.fd(0.025, .)

pcaFun <- PACE(curvesFun, npc = 2)
pcaFun %>% summary()
pcaFun$mu %>% plot()
pcaFun$functions[1] %>% plot()

# Prop of explained var
pcaFun$values / sum( pcaFun$values)
# all.equal(integrate(pcaFun$estVar), sum(pcaFun$values))

sdFun <- pcaFun$scores %>% apply(2, sd) 

# plots
PCcurves <- expand_grid(PC = 1:2,
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, fractionOfStDev) %>%
  reframe(time = pcaFun$mu@argvals[[1]],
            y = (pcaFun$mu + fractionOfStDev * sdFun[PC] * pcaFun$functions[PC])@X %>% as.numeric()
  )

ggplot(PCcurves) +
  aes(x = time, y = y, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  facet_grid(PC ~ .,
             scales = "free_y",
             labeller = labeller(PC = ~ str_glue("PC{.x}"))) +
  labs(color = expression(frac(s[k], sigma[k]))) +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0), color = 'black', linewidth = 1.5) +
  ggtitle("FPCA") +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "bottom")


# single PC(t)s
PC_labeller <- as_labeller(function(x) paste0('PC', x))
PCcurves <- tibble(PC = 1:2) %>% 
  group_by(PC) %>% 
  reframe(time = pcaFun$mu@argvals[[1]],
          value = pcaFun$functions[PC]@X %>% as.numeric()
)
ggplot(PCcurves) +
  aes(time, value) +
  geom_line() +
  facet_grid(PC ~ .,
             scales = "free_y",
             labeller = labeller(PC = PC_labeller)) +
  geom_hline(yintercept=0, linetype="dashed", color="red")



PCscores <- pcaFun$scores %>% `colnames<-`( paste0("s", 1:2)) %>% as_tibble %>%
  bind_cols(curves %>% distinct(curveId, Category))

# plots
Category.colors <- c("slategray4", "orangered")
# scatterplot PC scores s1 and s2 by Category
ggplot(PCscores) +
  aes(x = s1, y = s2, color = Category) +
  geom_point() +
  scale_color_manual(values=Category.colors) +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "bottom")

# boxplot PC scores s1 and s2 by Category
ggplot(PCscores %>% pivot_longer(s1:s2, names_to = "score", values_to = "value")) +
  aes(Category, value) +
  geom_boxplot() +
  facet_grid(score ~ ., scales="free_y") +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "bottom")





mod <- lm(s1 ~ Category, data = PCscores)
emmeans(mod, pairwise ~ Category)

emm <- emmeans(mod, pairwise ~ Category)$emmeans %>%
  as_tibble() %>% 
  select(Category, emmean) %>% 
  rename(s1= emmean)

predCurves <- emm %>% 
  group_by(Category) %>% 
  reframe(funData2long(pcaFun$mu + s1 * pcaFun$functions[1])) %>%
  select(-ID) %>% 
  rename(time = argvals, y = X)


ggplot(predCurves) +
  aes(time, y, color = Category) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  theme_light() +
  theme(text = element_text(size = 15),
        legend.position = "bottom")




  
# 2D

ex <- 1
curves <- read_csv(file.path("../data/", paste("ex2D", ex, "csv", sep = '.')))

ylist <- lapply(c('y1', 'y2'), function(y) {
  funData(
  argvals = curves %>% filter(curveId == 1) %>% pull(time),
  X = curves %>%
    complete(curveId,  time) %>% 
    pivot_wider(id_cols = curveId, names_from = time, values_from = {{y}}) %>% 
    select(!curveId) %>% 
    as.matrix()
)
})

y <- multiFunData(ylist)
pca2D <- MFPCA(y,
               M = 2,
               uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA"))
)
pca2D %>% plot()
pca2D %>% scoreplot()
