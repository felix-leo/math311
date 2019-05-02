##########################################################################
#
# Gaussian Chebyshev Quadrature
#   
##########################################################################

gausschebyshev = function(f,n) {
  ## The abscissa and weights are REALLY easy here
  ## The weights are pi/n and the abscissa are cos((2i-1)/(2n)*pi)
  xs = cos((2*(1:n)-1)/(2*n)*pi)
  return(pi/n*sum(f(xs)))
}


##########################################################################
#
# Gaussian Chebyshev Quadrature Examples
#   
##########################################################################


g=function(x) { tan(x)^2 }
gausschebyshev(g,3)

results=matrix(rep(0,30),ncol=1)
for (i in 1:30)
  results[i,1]=gausschebyshev(g,i)
results  #it doesn't change beyond 20!

gausschebyshev(g,20)

##########################################################################
