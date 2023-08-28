library(fda)
library(funData)
library(MFPCA)
library(tidyverse)

ex <- 1 # change according to ex number
curves <- read_csv(file.path("../data/", paste("ex1D", ex, "csv", sep = '.')))


curvesFun <- funData(argvals = seq(0, 2, by = 0.01),
                     X = curves %>%
                       complete(curveId,  time) %>% 
                       pivot_wider(id_cols = curveId, names_from = time, values_from = y) %>% 
                       select(!curveId) %>% 
                       as.matrix()
                     )

nCurves <- curvesFun %>% nObs() 
curvesFun %>% nObsPoints()
curvesFun %>% plot()

curvesFun[1] %>% approxNA() %>% funData2fd() %>% eval.fd(0.025, .)

pcaFun <- PACE(curvesFun, npc = 2)

sdFun <- pcaFun$scores %>% apply(2, sd) 

PCcurves <- expand_grid(PC = 1:2,
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, fractionOfStDev) %>%
  reframe(time = pcaFun$mu@argvals[[1]],
            value = (pcaFun$mu + fractionOfStDev * sdFun[PC] * pcaFun$functions[PC])@X %>% as.numeric()
  )

scores <- pcaFun$scores %>% `colnames<-`( paste0("s", 1:2)) %>% as_tibble %>%
  bind_cols(curves %>% distinct(curveId, Category))


