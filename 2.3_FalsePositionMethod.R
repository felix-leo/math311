########################################################################
#
# False Position Method
#   Implements this method for finding a root
#
########################################################################

falseposition = function(p0, p1, f, eps=1e-7, n=50, iter=FALSE) {
  ##
  ## Inputs
  ##   [p0, p1] = real root is between p0 and p1
  ##   f = function to find the root of
  ##   eps = tolerance (defaults to 1e-9 = .000000001)
  ##   n = iterations allowed (defaults to 30)
  
  Q0 = f(p0)
  Q1 = f(p1)
  
  failureMode = TRUE
  labelName = c("new_p", "[p0", "p1]", "err (p-p1)")
  perIteration = matrix(0,0,length(labelName))  #initialize a matrix to save iterations
  for (i in 1:n) {
    p = p1 - Q1*(p1-p0)/(Q1-Q0)
    perIteration = rbind(perIteration,c(p, p0, p1, abs(p-p1))) #save stuff
    if ( f(p)==0 || abs(p-p1) < eps ) { failureMode=FALSE; break }
    Q = f(p)
    if ( Q*Q1 < 0 ) { p0=p1; Q0=Q1 } #switch p0&p1 sometimes
    p1=p; Q1=Q 
  }
  if (failureMode) warning(paste("Failed to converge after",i,"iterations"))
  if (iter) { dimnames(perIteration)=list(1:i,labelName); p=perIteration }
  return(p)
}


############################################################################
## Testing the code
## Function to find the zero of
f = function(x) { x^4+2*x^2-3*x-3 }

## Examples
falseposition(1,2,f,iter=TRUE)

############################################################################
