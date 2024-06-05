library(fda)
library(funData)
library(MFPCA)
library(tidyverse)
library(landmarkregUtils)
library(emmeans)
library(mgcv)
library(itsadug)

zeroFun <- function(argvals) {
  argvals <- as.numeric(argvals)
  return(funData(argvals, matrix(0, ncol = length(argvals))))
}

reconstruction <- function(scores, basis, mu=NULL) {
  if (is.null(mu)) {
    res <- zeroFun(basis@argvals[[1]])
  } else {
    res <- mu
  }
  for (i in seq_along(scores)) {
    res <- res + (scores[i] * basis[i])
  }
  return(res)
}

mytheme <- theme_light() +
  theme(text = element_text(size = 16))

Category.colors <- c( "darkslategray",  "orangered")   #"firebrick1"  "slategray4" "cadetblue",


plots_dir <- "presentations/plots/"
data_dir <- "data/"

ex <- 1 # change according to ex number
raw_curves <- readRDS(file.path(data_dir, str_c("ex1D", ex, "rds", sep = '.'))) %>% ungroup() %>% 
  mutate(across(c(curveId, Category), ~ factor(.x)))


sp <- 0.01 # unified sampling period
maxT <- max(raw_curves$time)
grid <- seq(0, maxT, by = sp) # unified sampling grid 
curves <- raw_curves %>% 
  group_by(curveId, Category) %>% # all the factors at the level of curveId or higher (e.g. speaker)
  reframe(approx(time, y, grid) %>% as_tibble()) %>% # linear interpolation on grid 
  ungroup() %>% 
  rename(time = x, y = y)

pl <- raw_curves %>% 
  filter(curveId==51) %>% 
  ggplot() +
  aes(time, y) +
  geom_point() +
  ylim(ylim) +
  xlab("time (s)") + ylab("y") +
  mytheme +
  theme(legend.position = "none")

pl

curvesFun <- long2irregFunData(curves, id = "curveId", time = "time", value = "y") %>% 
  as.funData()

arith.colors <- c("black",
                  # "red",
                  "blue") # op1 , op2, result

tx <- curvesFun %>% argvals() %>% `[[`(1)
f2 <- funData(tx, matrix(10 * (tx - 0.25), nrow = 1))
plot(f2)

id <- 51
operands <- list( # change operands here
  curvesFun[id],
  fpca$mu 
  # + fpca$scores[id, 1] * fpca$functions[1]
  # + fpca$scores[id, 2] * fpca$functions[2]
  # + fpca$scores[id, 3] * fpca$functions[3]
  # 0.5 * (curvesFun[51] + curvesFun[1])
  # f2,
) %>% 
  lapply(function(f) {
    funData2long(f) %>% select(!id)
  }) %>% 
  bind_rows(.id = "id") %>% 
  mutate(id = factor(id))

ylim <- c(-3.8, 4)
ggplot(operands) +
  aes(time , value, group = id, color = id) +
  geom_line(linewidth = 1) +
  scale_color_manual(values=arith.colors) +
  ylim(ylim) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("time (s)") + ylab("y") +
  mytheme +
  theme(legend.position = "none")



ggsave(file.path(plots_dir, str_c("ex1D.1_GAM_smooth", '.png')), # pl,
       width = 1500, height = 1200, units = "px"
)

# plot a few curves
subset_curveId <- raw_curves %>%
  ungroup() %>% 
  distinct(curveId) %>%
  slice_sample(n = 20)

pl <- ggplot(
  # fpca$mu %>%
  # fpca$functions[3] %>% 
    # funData2long(time = "time", value = "y", id = "curveId")
  curves %>% inner_join(subset_curveId, by = "curveId")  %>%
               group_by(curveId) %>%
               mutate(y = y - fpca$mu %>% funData2long1() %>% pull(value))
             ) +
  aes(x = time, y = y, group = curveId) + #, color = Category) +
  ylim(ylim) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  theme(legend.position = "bottom")

pl

fpca <- PACE(curvesFun)

id <- 1
fpca$scores[id, 1:3] %>% round(2)

pl <- fpca$functions[1:3] %>% 
  funData2long(value = 'y', id = 'PC') %>% 
  mutate(PC = factor(PC, labels = str_c('PC', 1:3))) %>% 
  ggplot() +
  aes(x = time, y = y) + 
  ylim(ylim) +
  geom_line(linewidth = 1) +
  facet_grid(~ PC) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  mytheme  +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("ex1D.1_PCcolor", '.png')), #pl,
       width = 2500, height = 1200, units = "px"
)

ggplot(curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = y, group = curveId) +
  geom_line() +
  ylim(ylim) +
  mytheme


