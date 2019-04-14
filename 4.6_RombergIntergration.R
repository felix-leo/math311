############################################################################
#
# Romberg Integration - This evaluates the integral of f over the
#   interval [a,b] using n as the number of nodes evaluated
#  written by Scott Hyde
############################################################################

romberg = function(a, b, f, n=5) {
  #
  # Inputs
  #   a,b = endpoints of interval
  #   f   = function to integrate
  #   n   = number of rows of Romberg to generate
  h = b-a
  R = matrix(NA, n, n)
  R[1,1] = h*(f(a) + f(b))/2 #Trapezoid Rule n=1
  for (i in 2:n) {
    k=1:2^(i-2) #k is a vector
    R[i,1]=1/2*(R[i-1,1] + h*sum(f(a+(k-.5)*h)))
    for (j in 2:i) {
      R[i,j] = R[i,j-1] + (R[i,j-1] - R[i-1,j-1])/(4^(j-1)-1)
    }
    h = h/2
  }
  return(list(R=R,intg=R[n,n])) #return the R matrix and best approx
}

# f=function(x) { x^2*log(x)}
# rm=romberg(1,1.5,f,4);rm
# 
# f=function(x) {sqrt(1+cos(x)^2)}
# rm=romberg(0,48,f,10);rm
# 
# rm=romberg(0,pi/2,f,6);rm

rm2=romberg(30*pi/2,48,f,6);rm2

