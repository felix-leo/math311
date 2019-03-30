########################################################################
#
# Steffensen's Method 
#   Implements the Steffensen method for finding a fixed point
#   This method is quadratically convergent without computing
#   derivatives at the expense of two function evaluations per step.
#
########################################################################

steff = function(p0, g, eps=1e-7, n=50, iter=FALSE, complete=FALSE) {
  ##
  ## Inputs
  ##   p = initial guess
  ##   f = function to find the root of
  ##   eps = tolerance (defaults to 1e-9 = .000000001)
  ##   n = iterations allowed (defaults to 30)
  ##   complete = include fixed pt iterations too.
  ## Outputs
  ##   fixed point of g
  aitken_s = function(vector) { vector[1]-(vector[2]- vector[1])^2 /(vector[3] - 2*vector[2] + vector[1]) }
  failureMode=TRUE
  labelName = c("current_p", "p_1", "p_2",  "abs(p - p0)")
 
  p1 = NULL
  p2 = NULL
  
  perIteration = matrix(c(p0, NA, NA, NA), 1, length(labelName))  #save iterations
  
  for (i in 1:n) {
    p1 = g(p0);  p2 = g(p1);     #fixed point twice
    
    p = aitken_s(c(p0, p1, p2))  #Aitken's method once
    
    perIteration = rbind(perIteration, c(p0, p1, p2, abs(p - p0))) #save stuff
    if ( abs(g(p)-p) < eps ) {failureMode=FALSE; break }
    p0 = p
  }
  
  if (failureMode) warning(paste("Failed to converge after", i, "iterations"))
  
  if (iter) { 
    dimnames(perIteration)=list(0:i, labelName);
    return(perIteration) 
  }
  
  return(p)
}

############################################################################
## Testing the code
## Function to find the fixed point of
g = function(x) { cos(x) }

# g = function(x){ sqrt(10/(x+4)) }

steff(1,g)
bb=steff(1,g,iter=TRUE)
bb
noquote(formatC(bb,digits=11,format="f")) #error is formatted differently

steff(1,g,1e-15)
steff(1,g,1e-15,100,iter=TRUE)
steff(1,g,1e-15,100,iter=TRUE,complete=TRUE)
############################################################################
