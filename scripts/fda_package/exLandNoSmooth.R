# Exercise on Functional PCA 
# Author: Michele Gubian
# Date: 15 March 2023

# adjust path or setwd
source('header.R')
ex = 1 # change according to ex number

# load data
curves <- read_csv(file.path(data_dir,  paste("exLand", ex, "csv", sep = '.')))
land <- read_csv(file.path(data_dir,  paste("exLand", ex, "land", "csv", sep = '.')))


# plot a few curves (one by one + landmark position)
curveSample <- land %>% pull(curveId) %>% sample(5)
ggplot(curves %>% filter(curveId %in% curveSample)) +
  aes(x = time, y = y, color = Category, group = curveId) +
  geom_line() +
  geom_vline(data = land %>% filter(curveId %in% curveSample),
             mapping = aes(xintercept = l2, color = Category)) 

# landmark reg (based only on land, not on curves!)
reg <- landmarkreg.nocurve(land %>% select(starts_with("l")) %>% as.matrix,
                           nhknots = 8, hlambda=1e-8, wlambda =1e-8)

# Apply registration directly on the samples without performing any smoothing.
# In the general case, each curve has a different duration and a different 
# number and location of time samples.

# The operation consists in adding column 'regTime' to curves which reports 
# the mapping from each original time point on 'time' column to the registered time axis.

# Recover number of basis from reg as a reasonable guess, i.e.
# the number of basis used to model h() should be ok for inverse h() too.

nBasis <- reg$warpfd$basis$nbasis
Lfdobj <- 2 # 2 + order of derivative expected to be used. E.g. order = 1 if you need 1st deriv 
nOrder <- 2 + Lfdobj  # a fixed relation about B-splines
lambda <- 1e-6 # a guess, will be empirically checked later

curves %<>%
  group_by(curveId) %>% # curveId should be ordered in the same way as in land
  mutate(regTime = {
    basis_h_inv <- create.bspline.basis(range(time),nBasis,nOrder)
    fdPar_h_inv <- fdPar(basis_h_inv,Lfdobj,lambda)
    h_inv <- smooth.basis(reg$hfunmat[,cur_group_id()], reg$x, fdPar_h_inv)$fd
    eval.fd(time, h_inv) %>% drop()
  }) %>% 
  mutate(regTime = case_when( # correct rounding errors at boundaries
    time == min(time) ~ reg$warpfd$basis$rangeval[1],
    time == max(time) ~ reg$warpfd$basis$rangeval[2],
    TRUE ~ regTime
  )) %>%
  ungroup()
  
# Tune lambda by eye inspection of plots of time vs regTime and their fd object in reg
id <- 10 # pick a curve
plot(x = curves %>% filter(curveId == id) %>% pull(regTime), 
     y = curves %>% filter(curveId == id) %>% pull(time),
     col = "green",
     type = 'p')
lines(reg$warpfd[id], lwd = 2)

# plot the curve before/after landmark reg
curves %>% 
  filter(curveId == id) %>%
  pivot_longer(cols = c(time, regTime), names_to = "Time axis", values_to = "time") %>% 
  ggplot(aes(x = time, y = y, group = `Time axis`, color = `Time axis`)) +
  geom_line()



