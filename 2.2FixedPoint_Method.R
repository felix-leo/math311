########################################################################
#
# Fixed Point Method
#   Implements the fixed point method for finding a fixed point of g(x)
#     This gives g(p)=p, which also solves f(p)=0
#
########################################################################

fixedPoint = function(x0, g, eps=1e-7, n=30, iter=FALSE) {
  ##
  ## Inputs
  ##   x0 = initial guess
  ##   g = function to find the fixed point of
  ##   eps = tolerance (defaults to 1e-9 = .000000001)
  ##   n = iterations allowed (defaults to 30)
  ##   finite = number of iterations to perform before backing out.
  
  failureMode = TRUE
  p_Old = x0 #first guess
  labelName = c("p_Old", "abs (P - p_Old)") #make labelName as a vector
  perIteration = matrix(c(x0,0), 1, length(labelName))  #matrix to save iterations
  
  for (i in 1:n) {
    p = g(p_Old)
    perIteration = rbind(perIteration, c(p, abs(p - p_Old))) #save stuff
    if ( g(p)==p || abs(p - p_Old) < eps ) { failureMode=FALSE; break }
    p_Old = p
  }
  
  if (failureMode) warning(paste("Failed to converge after",i,"iterations"))
  if (iter) { 
    dimnames(perIteration) = list(0:i, labelName); 
    p=perIteration # swap the answer to per iteration matrix
  }
  return(p) # bring out the answer
}

############################################################################
# To perform fixed point for only a few iterations:

fxpt = function(x0, g, n) {
  ## Inputs
  ##  x0 = initial guess
  ##  g = function to find the fixed point of
  ##  n = number of iterations
  
  #the x0's will get replaced
  p=rep(x0, n+1) # rep(), repeat elements of a vector
  for (i in 1:n) {
    p[i+1] = g(p[i])
  
    piter=matrix( c(p, c(0, abs(diff(p))) ), ncol=2)
    colnames(piter)=c("p_Old", "abs(Î”p_n)")
    rownames(piter)=0:(dim(piter)[1]-1) #Iteration numb.
  }
  return(piter)
}

############################################################################
# Testing the code
# g=function(x) { (3+x-2*x^2)^(1/4) }

# fxpt(1,g,4)  #apply the fixed point for 4 iterations.

# fixedPoint(1,g,iter=TRUE)
# fixedPoint(1,g,n=1000,iter=TRUE)
############################################################################
