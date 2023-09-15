# testing landmark reg
library(fda)
library(funData)
library(MFPCA)
library(landmarkregUtils)
library(tidyverse)
library(emmeans)
library(parallel)
library(doParallel)

mytheme <- theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "bottom")

Category.colors <- c("slategray4", "orangered")

inputMarks <- c(0, 1,  2) * 100
targetMarks <-  c(0, 1.9,  2) * 10
h <- one_landmarkreg_nocurves(inputMarks, targetMarks)

plot(seq(targetMarks[1], targetMarks[length(targetMarks)], length.out = length(h)),
     h,
     type = 'l')
abline(v = targetMarks, col ='red')
abline(h = inputMarks, col ='red')

t <- seq(targetMarks[1], targetMarks[length(targetMarks)], length.out = length(h))
h_fd <- funData(t, matrix(h, nrow = 1))
autoplot(h_fd)

t_reg <- landmarkreg_timeSamples(10*(0:20), inputMarks, targetMarks)
points(t_reg, 10*(0:20), col = 'blue')

inputMarks <- matrix(c(0, 0.5, 1, 2, 0, 0.7, 1.2, 1.9, 0, 0.4, 1.1, 2.2), nrow=3, byrow = T)

for (i in 1:5) inputMarks <- rbind(inputMarks, inputMarks)

system.time({
reg <- landmarkreg_nocurves(inputMarks, njobs = 1)
})

ex <- 1 # change according to ex number
curves <- read_csv(file.path("../data/", paste("ex1D", ex, "csv", sep = '.')))
curves3 <- curves %>% 
  filter(curveId %in% 1:3)

curves3 <- curves3 %>% filter(
  curveId == 1 & time < 1.5 | curveId == 2 & time < 1.7 | curveId == 3
)

curves3 <- curves3 %>% mutate(curveId = as.factor(curveId))


land3 <- tribble(
~curveId,  ~l1, ~l2, ~l3,
  1, 0, 1, 1.5,
  2, 0, 1.7/2, 1.7,
  3, 0, 4/3, 2
) %>% 
  mutate(curveId = as.factor(curveId))

ggplot(curves3) +
  aes(time, y, color = curveId) +
  facet_grid(curveId ~ .) +
  geom_line() +
  geom_vline(aes(xintercept = value),
             data = land3 %>% pivot_longer(cols = starts_with("l"))) +
  mytheme



reg <- landmarkreg_nocurves(land3 %>% select(!curveId), c(0,1,2), compute_hinv = TRUE)
reg$h %>% plot()
reg$h %>% fd2funData(seq(0, 2, by=0.1)) %>% autoplot()

curves3reg <- applyReg(curves3, reg, grid = seq(0, 2, by=0.1),
                       id = "curveId", time = "time" , value = "y")

#### exLand

reg <- landmarkreg_nocurves(land %>% select(starts_with("l")), compute_hinv = TRUE)
t <- seq(0, 2, by=0.1)
curvesReg <- applyReg(curves %>% mutate(curveId=factor(curveId)), reg, grid = t,
                      id = "curveId", time = "time" , value = "y")

curvesReg <- curvesReg %>% inner_join(curves %>% distinct(curveId, Category), by = "curveId")

ggplot(curvesReg) +
  aes(time, y, group = curveId, color = Category) +
  geom_line()

curvesRegFd <- long2irregFunData(curvesReg) %>% as.funData()
curvesRegFd

lograteFd <- reg$logvelfd %>% fd2funData(t)

curvesRegFd2D <- multiFunData(curvesRegFd, lograteFd)

curvesFPCA <- MFPCA(curvesRegFd2D,
                    M = 2,
                    uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA"))
)

PCscores <- curvesFPCA$scores %>%
  `colnames<-`(str_c('s', 1:2)) %>% 
  as_tibble() %>% 
  bind_cols(curves %>% distinct(curveId, Category), .)


ggplot(PCscores %>% pivot_longer(cols = s1:s2, names_to = "score", values_to = "value")) +
  aes(Category, value) +
  facet_wrap(~ score, ncol = 1, scales = "free_y") +
  geom_boxplot() +
  mytheme

mod <- lm(s1 ~ Category, data = PCscores)
emm <- emmeans(mod, pairwise ~ Category)$emmean %>% 
  as_tibble() %>% 
  select(Category, emmean) %>% 
  rename(s1 = emmean)

predCurves <- emm %>% 
  group_by(Category) %>% 
  reframe(funData2long(curvesFPCA$meanFunction[[1]] + s1 * curvesFPCA$functions[[1]][1])) %>%
  select(-id) %>% 
  rename(y = value)

ggplot(predCurves) +
  aes(time, y, color = Category) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  mytheme

# compute durations
# example on a lograte curve:
lograteFd[1]
# landmark on reg time is at 0.5
# duration first segment:

lograteFd[1] %>% `*`(-1) %>% exp() %>% funData2fd() %>% defint.fd(c(0, 0.5))
# write a function in landmarkregUtils
# and copy defint.fd 

lr <- curvesFPCA$meanFunction[[2]] - 0.596 * curvesFPCA$functions[[2]][1]
lr %>% `*`(-1) %>% exp() %>% funData2fd() %>% defint.fd(c(0, 0.5))
