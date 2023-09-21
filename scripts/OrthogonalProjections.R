library(fda)
library(funData)
library(MFPCA)
library(tidyverse)
library(landmarkregUtils)
library(emmeans)
library(mgcv)
library(itsadug)


plotIntegral <- function(fd, plusColor = 'darkblue', minusColor = 'lightskyblue') {
  return(
    fd %>% 
      funData2long1() %>% 
      ggplot() +
      aes(time, value) +
      geom_line() +
      geom_ribbon(aes(ymax = ifelse(value > 0, value, 0), ymin = 0), fill = plusColor) +
      geom_ribbon(aes(ymin = ifelse(value < 0, value, 0), ymax = 0), fill = minusColor) +
      xlab("") +
      ylab("") 
  )
}

gplots::col2hex("lightskyblue") 
# "#00008B"  "#87CEFA"

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

Category.colors <- c("slategray4", "orangered")


plots_dir <- "presentations/plots/"
data_dir <- "data/"

pl <- ggplot() +
  xlim(-3, 3) +
  ylim(-1, 10) +
  xlab("x") + ylab("y") +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_function(fun = ~ .x **2 + 3, linewidth=1.5) +
  annotate(geom = "segment", x=1.7, xend=1.7, y=0, yend=5.89, linetype = 'dashed',
           color = 'blue') +
  annotate(geom = "segment", x=0, xend=1.7, y=5.89, yend=5.89, linetype = 'dashed',
           color = 'blue') +
  annotate(geom = "text", x=1.7, y = -0.5, label = "1.7", color = 'blue', size = 6) +
  annotate(geom = "text", x=-.5, y = 5.89, label = "5.89", color = 'blue', size = 6) +
  annotate(geom = "text", x=1, y=7.5, label="y = f(x)", size = 6) +
  mytheme

ggsave(file.path(plots_dir, str_c("f", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

ex <- 1 # change according to ex number
curves <- read_csv(file.path(data_dir, paste("ex1D", ex, "csv", sep = '.'))) %>% 
  mutate(across(c(curveId, Category), ~ factor(.x)))

# curves <- curves %>% 
#   mutate(Category = fct_recode(Category, ONE_PEAK = "NO_PEAK", TWO_PEAKS = "PEAK")) 
  
# write_csv(curves, file.path(data_dir, paste("ex1D", ex, "csv", sep = '.')))


# GAM
mod <- bam(y ~ Category + s(time, by = Category),
           data = curves  %>% 
             mutate(Category = factor(Category)))
plot_smooth(mod, view = "time", plot_all = "Category", col = Category.colors)
plot_diff(mod, view = "time", comp = list(Category = c("ONE_PEAK", "TWO_PEAKS")))
# save manually, 750x500

curvesFun <- long2irregFunData(curves, id = "curveId", time = "time", value = "y") %>% 
  as.funData()


curve <- curvesFun[1] %>% approxNA()
curve2 <- curvesFun[51] %>% approxNA()

# Orth basis
argvals<-seq(0,2,0.01)
ob <- eFun(argvals,M=20,type="Poly")


arith.colors <- c("black", "red", "blue")

# Compose func operation elements



operands <- list( # change operands here
  curve # op1 
  , ob[2] # op2  ,0.2 * ob[2]  curve2
  # ,curve * ob[2] # , 0.5 * (curve + curve2)  # result
  ) %>% 
  lapply(function(f) {
    funData2long(f) %>% select(!id)
  }) %>% 
  bind_rows(.id = "id") %>% 
  mutate(id = factor(id))


ylim <- c(-0.3, 0.5)
ylim <- c(-1.2, 1.2)
pl <- ggplot(operands) +
  aes(time , value, group = id, color = id) +
  geom_line(linewidth = 1) +
  scale_color_manual(values=arith.colors) +
  ylim(ylim) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("time (s)") + ylab("y") +
  mytheme +
  theme(legend.position = "none")

ggsave(file.path(plots_dir, str_c("scalarprod1", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

pl <- plotIntegral(curve * ob[2]) +
  xlab("time (s)") + ylab("y") +
  mytheme +
  theme(legend.position = "none")




pl <- ggplot(curve2 %>% funData2long()) +
  aes(time , value) +
  geom_line(col = 'black', linewidth = 1) +
  ylim(ylim) +
  # geom_point(col = 'blue') +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("curve2", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

pl <- ggplot((0.5*(curve + curve2)) %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'black', linewidth = 1) +
  ylim(ylim) +
  # geom_point(col = 'blue') +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("curveMean2", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)



pl <- ggplot((0.1 * ob[3]) %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'black', linewidth = 1) +
  ylim(ylim) +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("add3", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

pl <- ggplot((0.1 * ob[3] + curve) %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'black', linewidth = 1) +
  ylim(ylim) +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("curve_plus_add3", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

pl <- ggplot((-0.5 * curve) %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'black', linewidth = 1) +
  ylim(ylim) +
  # geom_point(col = 'blue') +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("curvem05", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)






pl <- ggplot(ob[2] %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'black', linewidth = 1) +
  ylim(-1.3, 1.3) +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("mul2", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

pl <- ggplot((ob[2] * curve) %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'black', linewidth = 1) +
  # ylim(-0.5, 0.5) +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("curve_times_mul2", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)



pl <- ggplot(ob %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'red', linewidth = 1.5) +
  facet_grid(~ ID) +
  xlab("time") + ylab("") +
  mytheme
  
ggsave(file.path(plots_dir, str_c("Poly4", '.png')), pl,
       width = 2500, height = 1000, units = "px"
       )

for (i in 1:4) {
  pl <- ggplot((curve * ob[i]) %>% funData2long1()) +
    aes(time, value) +
    geom_line(col = 'black', linewidth = 1.2) +
    xlab("time") + ylab("") +
    mytheme
  ggsave(file.path(plots_dir, str_c("curvePoly", i, '.png')), pl,
         width = 1500, height = 1200, units = "px"
  )
}


for (i in 1:4) {
  pl <- plotIntegral(curve * ob[i]) +
    xlab("time") + ylab("") +
    mytheme
  ggsave(file.path(plots_dir, str_c("curvePoly", i, 'Int', '.png')), pl,
         width = 1500, height = 1200, units = "px"
  )
  
}

for (i in 1:4) {
  pl <- ggplot(ob[i] %>% funData2long1()) +
    aes(time, value) +
    geom_line(col = 'red', linewidth = 1.2) +
    xlab("time") + ylab("") +
    mytheme
  ggsave(file.path(plots_dir, str_c("Poly", i, '.png')), pl,
         width = 1500, height = 1200, units = "px"
  )
}

s <- scalarProduct(curve, ob)
s %>% round(2)

for (i in c(1:4, 8, 12, 16, 20)) {
  pl <- reconstruction(s[1:i], extractObs(ob, 1:i)) %>%
    funData2long1() %>% 
    ggplot() +
    aes(time, value) +
    geom_line(data = curve %>% funData2long1(), col = 'black', linewidth = 1) +
    geom_line(col = 'blue', linewidth = 1.2) +
    xlab("time") + ylab("") +
    mytheme
  ggsave(file.path(plots_dir, str_c("curveRecPoly", i, '.png')), pl,
         width = 1500, height = 1200, units = "px"
  )
}

pl <- reconstruction(s[1:4], extractObs(ob, 1:4)) %>%
  funData2long1() %>% 
  ggplot() +
  aes(time, value) +
  geom_line(col = 'blue', linewidth = 1) +
  xlab("time (s)") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("curveRecPoly.png")), pl,
       width = 1500, height = 1200, units = "px"
)

# shape descr

s <- scalarProduct(curve, ob)
s %>% round(2)

ylim <- c(-0.2, 0.8)

pl <- ggplot((curve) %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'black', linewidth = 1) +
  ylim(ylim) +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("curve_shape", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)



s <- scalarProduct(curve + 0.3, ob)
s[1:4] %>% round(2)

pl <- ggplot((curve + 0.3) %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'black', linewidth = 1) +
  ylim(ylim) +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("curve_plus03_shape", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)


curveRev <- curve
nSamp <- length(curve@argvals[[1]])
curveRev@X[1,] <- curve@X[1, nSamp:1]
s <- scalarProduct(curveRev, ob)
s[1:4] %>% round(2)

pl <- ggplot((curveRev) %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'black', linewidth = 1) +
  ylim(ylim) +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("curveRev_shape", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

pl <- ggplot(ob %>% extractObs(obs = 1:4) %>% funData2long()) +
  aes(time, value) +
  geom_line(col = 'red', linewidth = 1) +
  facet_wrap(~ ID, ncol = 1, labeller = labeller(ID = ~ str_glue("B{.x}(t)"))) +
  xlab("") + ylab("") +
  theme_light() +
  theme( axis.text=element_blank(),axis.ticks=element_blank())

ggsave(file.path(plots_dir, str_c("PolyAll4", '.png')), pl,
       width = 400, height = 1000, units = "px"
)

pl <- ob %>%
  extractObs(obs = 1:4) %>%
  funData2long() %>% 
  rename(B = ID, time = argvals, y = X) %>% 
  ggplot(aes(time, y)) +
  geom_line(col = 'red', linewidth = 1) +
  facet_wrap(~ B, ncol = 2, labeller = labeller(B = ~ str_glue("B{.x}(t)"))) +
  mytheme

ggsave(file.path(plots_dir, str_c("Poly2by2", '.png')), pl,
       width = 1600, height = 1200, units = "px"
)





# stats with orth proj

nCurves <- curves %>% select(curveId) %>% n_distinct()
ylim <- c(-0.2, 0.5)
# plot a few curves
curveSampleId <- sample(nCurves, 12)
pl <- ggplot(curves %>% filter(curveId %in% curveSampleId)) +
  aes(x = time, y = y, group = curveId, color = Category) +
  scale_color_manual(values=Category.colors) +
  geom_line(linewidth = 0.5) +
  ylim(ylim) +
  theme_light() +
  theme(text = element_text(size = 15),
        legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("curves", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

meanCurve <- meanFunction(curvesFun, na.rm = TRUE)
pl <- meanCurve %>% 
  funData2long() %>% 
  ggplot(aes(time, value)) +
  geom_line(linewidth = 1) +
  xlab("time") + ylab("y") +
  ylim(ylim) +
  mytheme

ggsave(file.path(plots_dir, str_c("meanCurve", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

pl <- funData2long((curvesFun %>% extractObs(curveSampleId)) - meanCurve) %>% 
  ggplot(aes(argvals, X, group = ID)) +
  geom_line(linewidth = 0.5) +
  xlab("time") + ylab("y") +
  ylim(ylim) +
  mytheme

ggsave(file.path(plots_dir, str_c("curvesCentred", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)


polyScores <- sapply(1:4, function(i) {
  scalarProduct(curvesFun - meanCurve, ob[i])
}) %>%
  `colnames<-`(str_c('s', 1:4)) %>% 
  as_tibble() %>% 
  bind_cols(curves %>% distinct(curveId, Category), .) %>% 
  filter(!if_any(starts_with("s"), ~ is.na(.x)))

write_csv(polyScores, "../presentations/data/polyScores.csv")

modList <- lapply(1:4, function(s) {
  lm(as.formula(str_glue("s{s} ~  Category")), data = polyScores)
})

emmList <- lapply(1:4, function(s) {
  emmeans(modList[[s]], pairwise ~ Category)
})

emm <- lapply(seq_along(emmList),
              function(k) {
                emmList[[k]]$emmeans %>%
                  as_tibble() %>% 
                  select(Category, emmean) %>% 
                  rename(!!str_c("s", k) := emmean)
              }) %>% 
  reduce(left_join, by = "Category")


predCurves <- emm %>% 
  group_by(Category) %>% 
  # reframe(funData2long(meanCurve + s1 * ob[1])) %>%
  reframe(funData2long(meanCurve + s1 * ob[1] + s2 * ob[2] + s3 * ob[3] + s4 * ob[4])) %>%
  # reframe(funData2long(meanCurve + mean(polyScores$s1, na.rm = TRUE) * ob[1] +
  #                        mean(polyScores$s2, na.rm = TRUE) * ob[2] +
  #                        mean(polyScores$s3, na.rm = TRUE) * ob[3] +
  #                        s4 * ob[4])) %>% 
  select(-ID) %>% 
  rename(time = argvals, y = X)



ylim <- c(-0.2, 0.5)
pl <- ggplot(predCurves) +
  aes(time, y, color = Category) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  ylim(ylim) +
  theme_light() +
  theme(text = element_text(size = 15),
        legend.position = "bottom")
  
ggsave(file.path(plots_dir, str_c("emmCurves_s1234", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

# FPCA

pcaFun <- PACE(curvesFun, npc = 4)
pl <- pcaFun$functions %>% 
  funData2long() %>% 
  rename(PC = ID, time = argvals, y = X) %>% 
  ggplot(aes(time, y)) +
  geom_line(col = 'red', linewidth = 1) +
  facet_wrap(~ PC, ncol = 2, labeller = labeller(PC = ~ str_glue("PC{.x}(t)"))) +
  mytheme

ggsave(file.path(plots_dir, str_c("PC", '.png')), pl,
       width = 1600, height = 1200, units = "px"
)

(pcaFun$values / sum( pcaFun$values)) %>% `*`(100) %>% round(2)

pcScores <- pcaFun$scores %>%
  `colnames<-`(str_c('s', 1:4)) %>% 
  as_tibble() %>% 
  bind_cols(curves %>% distinct(curveId, Category), .) %>% 
  filter(!if_any(starts_with("s"), ~ is.na(.x)))

write_csv(pcScores, "../presentations/data/pcScores.csv")

modPC <- lm(s1 ~ Category, data = pcScores)
emmPC <- emmeans(modPC, pairwise ~ Category)$emmeans %>%
  as_tibble() %>% 
  select(Category, emmean) %>% 
  rename(s1= emmean)

predCurves <- emmPC %>% 
  group_by(Category) %>% 
  reframe(funData2long(pcaFun$mu + s1 * pcaFun$functions[1])) %>%
  select(-ID) %>% 
  rename(time = argvals, y = X)

pl <- ggplot(predCurves) +
  aes(time, y, color = Category) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  ylim(ylim) +
  theme_light() +
  theme(text = element_text(size = 15),
        legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("emmCurves_PCs1", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

# 2D FPCA



#### create ex2D.3

ex <- 1 # change according to ex number
curves <- read_csv(file.path(data_dir, paste("ex2D", ex, "csv", sep = '.'))) %>% 
  mutate(across(c(curveId, Category), ~ factor(.x)))
# scramble curves in dim y2, leave them as is in y1
curves <- curves %>% 
  distinct(curveId) %>% 
  mutate(id2 = sample(curveId, replace = FALSE)) %>% 
  left_join(curves %>% select(!c(y1, Category)), by = "curveId", multiple = "all") %>%
  select(!curveId) %>% 
  rename(curveId = id2) %>% 
  inner_join(curves %>% select(!y2), by = c("curveId", "time")) %>% 
  arrange(curveId, time) %>% 
  relocate(y1, .before = y2) %>% 
  relocate(Category, .before = time)

##### splines

# rng    <- c(0,2) 
# nbasis <- 20
# norder <- 4
# basis2 <- create.bspline.basis(rng, nbasis, norder)
# fdParObj <- fdPar(fdobj = basis2, Lfdobj = 2, lambda = 0)
# curvesFun1 <- smooth.basis(argvals = curvesFun[1]@argvals[[1]],
#                            y = curvesFun[1]@X %>% as.numeric(),
#                            fdParobj = fdParObj)$fd
# 
# sb <- fd(diag(nrow = nbasis), basis2)
# sb <- fd2funData(sb, argvals)
# 
# 
# curvesFun <- approxNA(curvesFun)
