library(fda)
library(funData)
library(MFPCA)
library(tidyverse)

ex <- 6 # change according to ex number
curves <- read_csv(file.path("../data/", paste("ex1D", ex, "csv", sep = '.')))


curvesFun <- funData(#argvals = seq(0, 2, by = 0.01),
                     argvals = curves %>% filter(curveId == 1) %>% pull(time),
                     X = curves %>%
                       complete(curveId,  time) %>% 
                       pivot_wider(id_cols = curveId, names_from = time, values_from = y) %>% 
                       select(!curveId) %>% 
                       as.matrix()
                     )

curvesFunIrr <- irregFunData(
  argvals = curves %>%
    select(curveId, time) %>%
    pivot_wider(names_from = curveId, values_from = time, values_fn = list) %>%
    as.list() %>%
    unlist(recursive = FALSE, use.names = FALSE),
  X = curves %>%
    select(curveId, y) %>%
    pivot_wider(names_from = curveId, values_from = y, values_fn = list) %>%
    as.list() %>%
    unlist(recursive = FALSE, use.names = FALSE)
  )

nCurves <- curvesFun %>% nObs() 
curvesFun %>% nObsPoints()
curvesFun %>% plot()

curvesFun[1] %>% approxNA() %>% funData2fd() %>% eval.fd(0.025, .)

pcaFun <- PACE(curvesFun, npc = 2)
pcaFun %>% summary()
pcaFun$values / sum( pcaFun$values)
all.equal(integrate(pcaFun$estVar), sum(pcaFun$values))

sdFun <- pcaFun$scores %>% apply(2, sd) 

# red and blue PC polts
PCcurves <- expand_grid(PC = 1:2,
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, fractionOfStDev) %>%
  reframe(time = pcaFun$mu@argvals[[1]],
            value = (pcaFun$mu + fractionOfStDev * sdFun[PC] * pcaFun$functions[PC])@X %>% as.numeric()
  )
# see ex1D.R

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
