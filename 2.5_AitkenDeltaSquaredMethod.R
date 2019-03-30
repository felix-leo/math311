########################################################################
#
# Aitken's Delta Squared Method
#   Implements the Aitken's Delta Squared method for accelerating
#   the convergence of a sequence to a fixed point.
#
########################################################################

aitkens = function(p) {
  ## Inputs
  ##   p = sequence of numbers (a vector)
  ## Outputs
  ##   result = sequence of numbers accelerated by Aitken's method
  ##
  p_new = function(p) { p[1]-(p[2]-p[1])^2/(p[3]-2*p[2]+p[1]) }
  result = 0
  for (i in 1:(length(p)-2)) {
    result[i] = p_new( p[i: (i+2)] ) 
  }
  return(result)
}

############################################################################
## Test of fixed point for Aitken's Delta Squared Method.  Comparison
## to normal fixed point
g = function(x) { .5*sqrt(10-x^3) }

#Next few lines implement the fixed point method
#  (but it doesn't check for convergence)
iter = 30 #number of iterations
fixpt = 1 #fixpt=initial guess; values after will be added
for (i in 2:iter) fixpt[i] = g(fixpt[i-1]) #fixed point method

phat = aitkens(fixpt)

#Compare the fixpt with Aitkens (NA is added to end of phat so they
#have the same size
p=1
for (i in 2:300) p = g(p)  # find the fixed pt to analyze the absolute error

compare = cbind( fixpt, abs(fixpt-p), c(phat,NA,NA), abs(c(phat,NA,NA)-p ))
colnames(compare) = c("fixed pt", "abs err", "Aitkens", "abs error")
options(width=100)
noquote(formatC(compare, digits=18, format="f"))
############################################################################
#  g=function(x) { sqrt(exp(x)/3) }
#  fixpt = rep(.75,7) #fixpt[1]=initial guess; values after will be replaced by next line > for (i in 2:length(fixpt)) fixpt[i] = g(fixpt[i-1]) #one line fixed point method!
#  A=cbind(aitkens(fixpt)) #cbind combines vectors as columns of a matrix
#  rownames(A)=1:(length(fixpt)-2) #change the row names to the iteration count
#  colnames(A)=c("phat") #call the column name phat!
#  A
# 
# phat
# 1 0.9078586
# 2 0.9095675
# 3 0.9099169
# 4 0.9099888
# 5 0.9100037

############################################################################
