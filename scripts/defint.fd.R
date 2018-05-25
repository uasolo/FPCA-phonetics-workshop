defint.fd <- function(fdobj, rng=NULL)
{
# Carry out a definite integral for each of the uni-dimensional functions in fdobj from rng[1] to rng[2].
# Defaults to the entire interval.
# Returns a vector of length = number of functions in fdobj.

# TODO: input format check
if (is.null(rng)) { rng = fdobj$basis$rangeval }
return (apply(inprod(fdobj,fdobj$basis,rng=rng),1,sum))

}

