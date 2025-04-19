parab_tongue <- function(a_left, a_right, scale = c(1,1), centre = c(0, 0), npoints = 11,theta = 0) {
  # constructs a fake tongue contour by connecting two half parabolas at the vertex.
  # a_<side> are the coefficients of the quadratic term of the left/right halves,
  # which should be positive;
  # they are then multiplied by -1 so that the contour is concave (upside-down-U).
  # theta is a rotation angle pivoting at the vertex.
  
  # Returns a tibble with knot, x and y columns and as many rows as npoints,
  # where knot is a fake knot index imitating the output of AAA.
  x <- seq(-2, 2, length.out = npoints)
  y <- - c(a_left * (x[x<0] **2) , a_right * (x[x>=0] **2))
  tongue <- tibble(x=x, y=y) %>% 
    rowid_to_column("knot")
  # rotation
  if (theta != 0) {
    rot <- matrix(c( cos(theta), sin(theta), - sin(theta), cos(theta)), 
                  nrow = 2,
                  dimnames = list(c('x', 'y'),c('x', 'y')) )
    tongue <- bind_cols(tongue %>% select(knot),
                        tongue %>% select(c(x, y)) %>% as.matrix() %>% `%*%`(t(rot)) %>% as_tibble()
    )
  }
  # scale and shift
  tongue <- tongue %>% 
    mutate(x = x * scale[1] + centre[1],
           y = y * scale[2] + centre[2])
  return(tongue)
}