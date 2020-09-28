# Exercise on Functional PCA 
# Author: Michele Gubian
# Last revision: September 2020


# adjust path or setwd
source('header.R')
ex = 1 # change according to ex number

# load data
curves <- read_csv(file.path(data_dir,  paste("exLand", ex, "csv", sep = '.')))
land <- read_csv(file.path(data_dir,  paste("exLand", ex, "land", "csv", sep = '.')))


# plot a few curves (one by one + landmark position)
curveSample <- land %>% pull(curveId) %>% sample(20)
ggplot(curves %>% filter(curveId %in% curveSample)) +
  aes(x = time, y = y, color = Category) +
  geom_line() +
  geom_vline(data = land, mapping = aes(xintercept = l2, color = Category)) 

# landmark reg (based only on land, not on curves!)
reg <- landmarkreg.nocurve(land %>% select(starts_with("l")) %>% as.matrix,
                           nhknots = 8, hlambda=1e-8, wlambda =1e-8)

# plot time warps and log rates (just checking)
reg$warpfd %>% plot
reg$logvelfd %>% plot

# create common basis
nKnots <- 12 # try many
lambda <- 1e-6 # try many
curveRange <- reg$land %>% range

Lfdobj <- 2 # 2 + order of derivative expected to be used. E.g. order = 1 if you need 1st deriv 
nOrder <- 2 + Lfdobj  # a fixed relation about B-splines
nBasis <- nKnots + nOrder - 2 # a fixed relation about B-splines
basis <- create.bspline.basis(curveRange, nBasis, nOrder)

fdParObj <- fdPar(fdobj = basis, Lfdobj = Lfdobj, lambda = lambda)

coef <- curves %>% group_by(curveId) %>%
  summarise(coef = {
    range_i <- range(c( reg$hfunmat[,curveId], time) ) # prevent rounding errors
    basis_i <- create.bspline.basis(range_i,nBasis,nOrder)
    fdPar_i <- fdPar(basis_i,Lfdobj,lambda)
    y_nonreg_fd <- smooth.basis(time,y,fdPar_i)$fd
    # the next line is illustrated in notes/landmark_reg.pdf
    c(smooth.basis(reg$x, eval.fd(reg$hfunmat[,curveId[1]], y_nonreg_fd),fdParObj)$fd$coefs)
  },
  coefId = seq_len(nBasis)) %>% # coefId to make pivot_wider work (formerly 'spread') 
  # fd() requires coef to be a nBasis-by-nCurves matrix
  pivot_wider(names_from = curveId, values_from = coef) %>%
  select(-coefId) %>%
  as.matrix

y_fd = fd(coef=coef, basisobj=basis) # all curves in one fd object

plot(y_fd) # just checking
abline(v = reg$land[2])

# jump to ex1D.R
# compare plot(y_fd) you obtain ignoring registration
# (where possible, or take care of linear time normalisation first)
# carry on in ex1D.R to analyse registered curve in 1D style

# Joint FPCA on registered curves and their relative speech rate
# i.e. 2D style analysis, with rime warp logvel as second dimension

# Build a 2D fd object

# convert logvelfd to basis
tx <- reg$x # or any other sampling comb
logvelfd <- smooth.basis(tx, eval.fd(tx, reg$logvelfd),fdParObj)$fd
# assemble (same basis) spline coefficients
# reminder: in multi-dim case, fd() requires coef to be a  nBasis-by-nCurves-by-nDimensions array 
Y_coefs <- array(dim = c(dim(y_fd$coefs),2))
Y_coefs[,,1] <- y_fd$coefs
Y_coefs[,,2] <- logvelfd$coefs 
# make the two different units comparable
w2 <-  mean(eval.fd(tx,sd.fd(y_fd))) / mean( eval.fd(tx,sd.fd(logvelfd)))
Y_coefs[,,2] <- Y_coefs[,,2] * w2
# bundle all into a 2D fd object
y_fd <- fd(coef=Y_coefs, basisobj=basis)

# jump to ex2D.R for a 'classic' 2D analysis

# representing logvel as segments (i.e. inter-landmark intervals) durations
Nsegments <- (land %>% select(starts_with("l")) %>% ncol) -1
# PC curves -> durations
tx <- seq(0, 2, length.out = 35) # re-sampling smooth curves at regular intervals
# compute st dev of PC scores, will plot variation -/+ 1 st dev
sdScores <- y_pcafd %>% getPCscores %>% apply(2, sd) 
# construct example curves by applying reconstruction formula
PCdurations <- expand_grid(PC = 1:2,
                           fractionOfStDev = seq(-1, 1, by=.25),
                           Interval = seq_len(Nsegments)) %>% # n_segments = n_landmarks -1
  group_by(PC, fractionOfStDev, Interval) %>%
  summarise(# linear combination of spline coefs of mean + score * PC curve
            value = (y_pcafd$meanfd$coefs[,1,2] + # mean of dimension 2, i.e. logvel
                       fractionOfStDev * sdScores[PC] * # PC score
                       y_pcafd$harmonics$coefs[,PC, 2] * # PC curve  of dimension 2, i.e. logvel
                       (1/w2) * # undo the dimension weighting
                       (-1)) %>% # see reverse logvel formula
              fd(coef = ., basisobj = y_pcafd$meanfd$basis) %>% # make it a fd object
              eval.fd(tx, .) %>% # sample it at time = tx
              as.numeric %>%
              exp %>%  # see reverse logvel formula
              smooth.basis(tx, ., fdParObj) %>%
              .$fd %>%
              defint.fd(c(reg$land[Interval],reg$land[Interval+1])))  # see reverse logvel formula

PC_labeller <- as_labeller(function(x) paste0('PC', x))
ggplot(PCdurations) +
  aes(x = fractionOfStDev %>% factor(labels = ""), y = value, color = fractionOfStDev) + 
  geom_bar(stat="identity", fill = 'white') +
  geom_text(aes(label = value %>% round(digits = 2)), size = 5, position = position_stack(vjust = 0.5)) +
  facet_grid(~ PC, labeller = PC_labeller) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  labs(color = expression(frac(s[k], sigma[s[k]]))) +
  ylab("Duration") +
  xlab("") +
  coord_flip() + 
  theme_light() +
  theme(text = element_text(size = 18),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "bottom") 


