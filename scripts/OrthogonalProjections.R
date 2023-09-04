library(fda)
library(funData)
library(tidyverse)

funData2long1 <- function(fd) {
  return(tibble(argvals = fd@argvals[[1]],
                X = fd@X %>% as.numeric()))
}

funData2long <- function(fd) {
  return(
    tibble(ID = seq_len(fd@X %>% nrow())) %>% 
      group_by(ID) %>% 
      reframe(funData2long1(fd[ID])) %>% 
      ungroup() %>% 
      mutate(ID = factor(ID))
  )
}

plotIntegral <- function(fd, plusColor = 'darkgreen', minusColor = 'lawngreen') {
  return(
    fd %>% 
      funData2long1() %>% 
      ggplot() +
      aes(argvals, X) +
      geom_line() +
      geom_ribbon(aes(ymax = ifelse(X > 0, X, 0), ymin = 0), fill = plusColor) +
      geom_ribbon(aes(ymin = ifelse(X < 0, X, 0), ymax = 0), fill = minusColor) +
      xlab("") +
      ylab("") 
  )
}

gplots::col2hex("lawngreen")

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

plots_dir <- "../presentations/plots/"


segments <- tribble(
  ~x, ~xend, ~y, ~yend,
  1.7, 1.7, 0, 1.7**2 + 3,
  0, 1.7, 1.7**2 + 3, 1.7**2 + 3
)

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
curves <- read_csv(file.path("../data/", paste("ex1D", ex, "csv", sep = '.')))

curvesFun <- funData(argvals = seq(0, 2, by = 0.01),
                     X = curves %>%
                       complete(curveId,  time) %>% 
                       pivot_wider(id_cols = curveId, names_from = time, values_from = y) %>% 
                       select(!curveId) %>% 
                       as.matrix()
) %>% 
  approxNA()

curve <- curvesFun[1]

pl <- ggplot(curve %>% funData2long()) +
  aes(argvals, X) +
  geom_line(col = 'black', linewidth = 1) +
  xlab("time") + ylab("") +
  mytheme

ggsave(file.path(plots_dir, str_c("curve", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

argvals<-seq(0,2,0.01)
ob <- eFun(argvals,M=4,type="Poly")

pl <- ggplot(ob %>% funData2long()) +
  aes(argvals, X) +
  geom_line(col = 'red', linewidth = 1.5) +
  facet_grid(~ ID) +
  xlab("time") + ylab("") +
  mytheme
  
ggsave(file.path(plots_dir, str_c("Poly4", '.png')), pl,
       width = 2500, height = 1000, units = "px"
       )

for (i in 1:4) {
  pl <- ggplot((curve * ob[i]) %>% funData2long1()) +
    aes(argvals, X) +
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
    aes(argvals, X) +
    geom_line(col = 'red', linewidth = 1.2) +
    xlab("time") + ylab("") +
    mytheme
  ggsave(file.path(plots_dir, str_c("Poly", i, '.png')), pl,
         width = 1500, height = 1200, units = "px"
  )
}

s <- scalarProduct(curve, ob)
for (i in 1:4) {
  pl <- reconstruction(s[1:i], extractObs(ob, 1:i)) %>%
    funData2long1() %>% 
    ggplot() +
    aes(argvals, X) +
    geom_line(data = curve %>% funData2long1(), col = 'black', linewidth = 1) +
    geom_line(col = 'red', linewidth = 1.2) +
    xlab("time") + ylab("") +
    mytheme
  ggsave(file.path(plots_dir, str_c("curveRecPoly", i, '.png')), pl,
         width = 1500, height = 1200, units = "px"
  )
}






##### splines

rng    <- c(0,2) 
nbasis <- 20
norder <- 4
basis2 <- create.bspline.basis(rng, nbasis, norder)
fdParObj <- fdPar(fdobj = basis2, Lfdobj = 2, lambda = 0)
curvesFun1 <- smooth.basis(argvals = curvesFun[1]@argvals[[1]],
                           y = curvesFun[1]@X %>% as.numeric(),
                           fdParobj = fdParObj)$fd

sb <- fd(diag(nrow = nbasis), basis2)
sb <- fd2funData(sb, argvals)


curvesFun <- approxNA(curvesFun)
