########################################################################
#
# Newton's Method
#   Implements the Newton-Raphson method for finding a root
#
########################################################################

newtons = function(p0, f, fp, eps=1e-7, n=50, iter=FALSE) {
  ##
  ## Inputs
  ##   f = function to find the root of (global function)
  ##   fp = derivative of f (global function)
  ##   p0 = initial guess
  ##   eps = tolerance (defaults to 1e-9 = .000000001)
  ##   n = iterations allowed (defaults to 30)
  ##
  
  failureMode=TRUE
  labelName = c("P_value", "abs (P - P0)")
  perIteration = matrix(c(p0,NA),1,length(labelName))  #initialize a matrix to save iterations
  for (i in 1:n) {
    diff = f(p0)/fp(p0)
    p = p0 - diff
    perIteration = rbind(perIteration, c(p, abs(p- p0))) #save stuff
    if ( f(p)==0 || abs(p- p0) < eps ) { failureMode=FALSE; break }
    p0 = p
  }
  if (failureMode) warning(paste("Failed to converge after",i,"iterations"))
  if (iter) { dimnames(perIteration)=list(0:i,labelName); p=perIteration }
  return(p)
}


############################################################################
## Testing the code
## Function to find the zero of

#f = function(x) { x^2-3 }
#fp = function(x) { 2*x }

## Examples:
newtons(1)

bb=newtons(1,iter=TRUE)
bb
noquote(formatC(bb,digits=11,format="f")) #error is formatted differently

newtons(1,1e-15)
newtons(1,1e-15,100)
newtons(1,1e-15,100,iter=TRUE)
############################################################################
