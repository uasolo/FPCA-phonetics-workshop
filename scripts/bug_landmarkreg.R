library(fda)
packageVersion('fda')
# [1] ‘6.1.4’
# create a fd object with only one curve
basis <- create.bspline.basis(rangeval = c(0,1), nbasis = 4)
x <- seq(0, 1, by = 0.1)
y <- x ** 2
y_fd <- smooth.basis(x, y, basis)$fd

# landmark register with only one interior landmark,
# let the function determine WfdPar
reg <- landmarkreg(y_fd, 0.5, 0.7)
# Error in create.bspline.basis(rangeval, Wnbasis) : 
#   nbasis must be at least norder;  nbasis = 3;  norder = 4
# more than one interior landmark ok
reg <- landmarkreg(unregfd = y_fd,
                   ximarks = c(0.5, 0.6),
                   x0marks = c(0.6, 0.7)
)
                   
                   
reg <- landmarkreg(unregfd = y_fd,
                   ximarks = c(0.5, 0.6),
                   x0marks = 0.1 + c(0.5, 0.6),
                   x0lim = c(0, 1.5))




