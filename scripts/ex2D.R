# Exercise on Functional PCA 
# Author: Michele Gubian
# Last revision: September 2020


# adjust path or setwd
source('header.R')
library(abind)
ex = 2 # change according to ex number

# load data

curves <- read_csv(file.path(data_dir,  paste("ex2D", ex, "csv", sep = '.')))
# all curves start at time == 0 an end at time == 2
curveRange <- c(0, 2)
nCurves <- curves$curveId %>% n_distinct

# plot a few curves
ggplot(curves %>% filter(curveId %in% sample(nCurves, 20)) %>%
         pivot_longer(starts_with("y"), names_to = "DIM", values_to = "y")) +
  aes(x = time, y = y, group = curveId, color = Category) +
  facet_grid(DIM ~ .) +
  geom_line() +
  theme_light() +
  theme(text = element_text(size = 15),
        legend.position = "bottom")

# smoothing

# create common basis
nKnots <- 12 # try many
lambda <- 1e-6 # try many

Lfdobj <- 2 # 2 + order of derivative expected to be used. E.g. order = 1 if you need 1st deriv 
nOrder <- 2 + Lfdobj  # a fixed relation about B-splines
nBasis <- nKnots + nOrder - 2 # a fixed relation about B-splines
basis <- create.bspline.basis(curveRange, nBasis, nOrder)

fdParObj <- fdPar(fdobj = basis, Lfdobj = Lfdobj, lambda = lambda)

# pick nKnots and lambda by visual inspection of single curves
# here we just construct y1 or y2 curves separately, for the purpose of nKnots and lambda selection
i <- 10 # pick one curveId
curve <- curves %>% filter(curveId == i)
curve_fd <- smooth.basis(argvals = curve$time, y = curve$y1, fdParobj = fdParObj)$fd # smoothed curve

# familiarise with fd objects
curve_fd %>% names
curve_fd$coefs

# here using base plot as it's quicker for just one curve (fda has a plot.fd function)
plot(curve_fd, ylim = c(-0.6, 0.1)) # adjust ylim
title(main = paste("nKnots =", nKnots, "; lambda =", lambda)) 
points(curve$time, curve$y1, col = "red") # original samples

# parameters set to:
# nKnots = 12
# lambda = 1e-6

# smooth all curves
# smooth.basis() does not accept different time samples for different curves.
# Thus we create smooth curves one by one on the same basis,
# store the spline coefficients and compose an fd object at the end.

# in multi-dim case, fd() requires coef to be a  nBasis-by-nCurves-by-nDimensions array 
# (in 2D exercises, nDimensions == 2 of course)

coef <- curves %>%
  pivot_longer(starts_with("y"), names_to = "DIM", values_to = "y") %>%
  group_by(curveId, DIM) %>%
  summarise(coef = c(smooth.basis(argvals = time, y = y, fdParObj)$fd$coefs),
            coefId = seq_len(nBasis)) %>% # coefId to make pivot_wider work (formerly 'spread') 
  pivot_wider(names_from = coefId, values_from = coef) %>%
  group_by(DIM) %>%
  group_map(~ select(.x, -curveId)) %>%
  abind(along = 3) %>%
  aperm(c(2,1,3))

# check class and dim
coef %>% class
coef %>% dim
 

y_fd = fd(coef=coef, basisobj=basis) # all curves in one fd object (same syntax as for 1D case)

op <- par(mfrow=c(2,1))
plot(y_fd) # just checking
par(op)

# FPCA
# same syntax as for 1D case
lambda_pca <- lambda 
pcafdPar  <- fdPar(fdobj = basis, Lfdobj = 2, lambda = lambda_pca) # here Lfdobj = 2 since we don't need derivatives of PC curves.
y_pcafd <- pca.fd(y_fd, nharm=2, pcafdPar) # compute first 2 PCs

# familiarise with pca.fd object
y_pcafd %>% names 
# "harmonics" stands for PCs
# "scores" can be extracted directly, but for higher dim curves use getPCscores (because of a bug in fda)
y_pcafd$varprop # percentage of explained variance by each PC

# plot mean and PC curves
par(mfrow=c(2,1))
plot(y_pcafd$meanfd) ; title(main = "mean(t)")
par(mfrow=c(1,1))
plot(y_pcafd$harmonics[1]) ; title(main = "PC1(t)")
plot(y_pcafd$harmonics[2]) ; title(main = "PC2(t)")

# plot PC scores effect using default plotting from fda 
op <- par(mfrow=c(2,2))
plot.pca.fd(y_pcafd, nx=40)
par(op)

# plot PC scores effect using ggplot (nicer but harder)
tx <- seq(0, 2, length.out = 35) # re-sampling smooth curves at regular intervals
# compute st dev of PC scores, will plot variation -/+ 1 st dev
sdScores <- y_pcafd %>% getPCscores %>% apply(2, sd) 
# construct example curves by applying reconstruction formula
PCcurves <- expand_grid(PC = 1:2,
                        DIM = c(1, 2), # new wrt 1D case
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, DIM, fractionOfStDev) %>%
  summarise(time = tx,
            # linear combination of spline coefs of mean + score * PC curve
            # note the extra dimension in coefs due to DIM
            value = (y_pcafd$meanfd$coefs[,1,DIM] + # mean 
                       fractionOfStDev * sdScores[PC] * # PC score
                       y_pcafd$harmonics$coefs[,PC, DIM]) %>% # PC curve
              fd(coef = ., basisobj = y_pcafd$meanfd$basis) %>% # make it a fd object
              eval.fd(tx, .) %>% # sample it at time = tx
              as.numeric) # otherwise you get a matrix as column (dunno why)

# actual plot
PC_labeller <- as_labeller(function(x) paste0('PC', x))
DIM_labeller <- as_labeller(function(x) paste0('y', x))
ggplot(PCcurves) +
  aes(x = time, y = value, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  facet_grid(DIM ~ PC,
             scales = "free_y",
             labeller = labeller(PC = PC_labeller, DIM = DIM_labeller)) +
  labs(color = expression(frac(s[k], sigma[k]))) +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0), color = 'black', size = 1.5) +
  xlab("Normalised time") +
  ylab("") +
  ggtitle("FPCA") +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "bottom")




# gather PC scores 
PCscores <- y_pcafd %>% getPCscores %>% `colnames<-`( paste0("s", 1:2)) %>% as_tibble %>%
  # curveId and Category corresponding to PC scores
  # (make sure there are no reorderings along the way!)
  bind_cols(., (curves %>% distinct_at(vars(curveId, Category)))) 

# scatterplot PC scores s1 and s2 by Category
ggplot(PCscores) +
  aes(x = s1, y = s2, color = Category) +
  geom_point() +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "bottom")



# FPCA-based curve reconstruction

		
# linear model 

lm1 <- lm(s1 ~ Category, data = PCscores)
lm1 %>% summary

# get estimated marginal means for each level of Category
emm <- emmeans::emmeans(lm1, pairwise ~ Category)


# aux function: reconstruct spline coef given fpcaObj, PC, dimension and scores.
# multi-D version (not general for 1D and multi-D)
# scores is a vector whose indices are PC indices
# Used internally in reconstructCurve() and reconstructDurations()
reconstructCoef <- function(fpcaObj, scores, dimension) {
  fpcaObj$meanfd$coefs[, 1, dimension] + 
     sapply(seq_along(scores), function(PCidx) {
       scores[PCidx] * fpcaObj$harmonics$coefs[, PCidx, dimension]
     }) %>% apply(1, sum)
}

# aux function: reconstruct a curve given fpcaObj, PC, dimension, scores and time samples
# multi-D version (not general for 1D and multi-D)
# scores is a vector whose indices are PC indices
reconstructCurve <- function(fpcaObj, scores, dimension, tx) {
  reconstructCoef(fpcaObj, scores, dimension) %>% 
    fd(coef = ., basisobj = fpcaObj$meanfd$basis) %>%
    eval.fd(tx, .) %>% 
    as.numeric
}


EMMcurves <- emm$emmeans %>%
  as_tibble %>%
  group_by(Category) %>%
  summarise(time = tx,
            # linear combination of spline coefs of mean + score * PC curve
            y1 = reconstructCurve(y_pcafd, emmean, 1, tx),
            y2 = reconstructCurve(y_pcafd, emmean, 2, tx)
  ) 

ggplot(EMMcurves %>% pivot_longer(cols = starts_with("y"), names_to = "DIM", values_to = "y")) +
  aes(time, y, color = Category, group = Category) +
  geom_line() +
  facet_grid(DIM ~ .) +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "bottom"
        )


# Trajectory representation
ggplot(EMMcurves) +
  aes(x = y1, y = y2, group = Category, color = Category) +
  geom_path(size = 1, arrow = arrow(length = unit(0.5, "cm"), type = "closed")) +
  theme_light() +
  theme(text = element_text(size = 16),
        legend.position = "bottom")



