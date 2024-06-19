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



ggsave(file.path(plots_dir, str_c("ex1D.4_land2D_pred_curve_land", '.png')), # pl,
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

ggsave(file.path(plots_dir, str_c("ex1D.4_land_curves_land", '.png')), #pl,
       width = 2500, height = 1200, units = "px"
)

ggplot(curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = y, group = curveId) +
  geom_line() +
  ylim(ylim) +
  mytheme


# lin time norm
raw_curves <- raw_curves %>%
  group_by(curveId) %>%
  mutate(time = time/max(time))

#### land reg ex1D 4

# lin norm 

land_y <- land %>%
  pivot_longer(l1:last_col(), names_to = "landmark", values_to = "time") %>% 
  inner_join(subset_curveId, by = "curveId") %>% 
  group_by(curveId) %>% 
  mutate(time = time/max(time)) %>% # lin norm
  mutate(y = {
    y <- curves %>%
      inner_join(cur_group(), by = 'curveId') %>% 
      pull(y)
    approx(grid, y, time, rule = 2)$y
  }) 


#### land reg ex 4 
# plot aligned landmarks in color dots

land_y <- tibble(landmark = reg$landmarks %>% names, time = reg$landmarks) %>% 
  expand_grid(subset_curveId) %>% 
  group_by(curveId) %>% 
  mutate(y = {
    y <- curvesReg %>%
      inner_join(cur_group(), by = 'curveId') %>% 
      pull(y)
    approx(seq(0, last(reg$landmarks), length.out = 100), y, time, rule = 2)$y
  }) 
  
ggplot(curvesReg %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = y, group = curveId) +
  geom_line(linewidth = 0.6, color = 'slategray4') +
  geom_point(data = land_y,
             mapping = aes(time, y, color = landmark, group = curveId),
             inherit.aes = FALSE,
             size = 2) +
  # facet_wrap(~ curveId) +
  scale_color_brewer(palette = "Dark2") + 
  xlab("Registered time") +
  # scale_color_manual(values=c('blue', 'green', 'magenta', 'black')) +
  mytheme  +
  theme(legend.position = "bottom")


ggplot(predCurves) +
  aes(time, y, color = Category) +
  geom_line() +
  geom_ribbon(aes(x = time, ymin = yl, ymax = yu, fill = Category),
              alpha = 0.3, inherit.aes = FALSE) +
  scale_color_manual(values=Category.colors) +
  scale_fill_manual(values=Category.colors) +
  geom_vline(xintercept = reg$landmarks) + 
  xlab("registered time") +
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = reg$landmarks %>% names())) +
  
  mytheme +
  theme(legend.position = "bottom")


### h, ex1D.4

# pick a short curve
i_short_long <- c(8, 98)

# with landmark position
land_y <- land %>%
  pivot_longer(l1:last_col(), names_to = "landmark", values_to = "time") %>% 
  filter(curveId %in% i_short_long) %>% 
  group_by(curveId) %>% 
  mutate(y = {
    y <- curves %>%
      inner_join(cur_group(), by = 'curveId') %>% 
      pull(y)
    approx(grid, y, time, rule = 2)$y
  }) 

ggplot(curves %>% filter(curveId %in% i_short_long)) +
  aes(x = time, y = y, group = curveId, color = Category) +
  geom_line(linewidth = 0.8) +
  geom_point(data = land_y,
             inherit.aes = TRUE,
             size = 3) +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  theme(legend.position = "none")

ggplot(curvesReg %>% filter(curveId %in% i_short_long)) +
  aes(time, y, group = curveId, color = Category) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  geom_vline(xintercept = reg$landmarks) + 
  xlab("registered time") +
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = reg$landmarks %>% names())) +
  mytheme +
  theme(legend.position = "none")


land %>% 
  select(curveId, Category) %>% 
  filter(curveId %in% i_short_long) %>% 
  group_by(curveId, Category) %>% 
  reframe(time = seq(0, last(reg$landmarks), length.out = 100),
          h = eval.fd(time, reg$h[cur_group()$curveId]) %>% as.numeric()) %>% 
  ggplot() +
  aes(time, h, group = curveId, color = Category) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  geom_vline(xintercept = reg$landmarks) + 
  xlab("registered time") +
  ylab("time") +
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = reg$landmarks %>% names())) +
  mytheme +
  theme(legend.position = "none")

  
land %>% 
  select(curveId, Category) %>% 
  filter(curveId %in% i_short_long) %>% 
  group_by(curveId, Category) %>% 
  reframe(time = seq(0, last(reg$landmarks), length.out = 100),
          lograte = eval.fd(time, reg$lograte[cur_group()$curveId]) %>% as.numeric()) %>% 
  ggplot() +
  aes(time, lograte, group = curveId, color = Category) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  geom_vline(xintercept = reg$landmarks) + 
  xlab("registered time") +
  ylab("Log rate") +
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = reg$landmarks %>% names())) +
  mytheme +
  theme(legend.position = "none")


#### TODO: move to function
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

########

