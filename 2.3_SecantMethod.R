########################################################################
#
# Secant Method
#   Implements the Secant method for finding a root
#
########################################################################

secant = function(p0, p1, f, eps=1e-7, n=50, iter=FALSE) {
  ##
  ## Inputs
  ##   p0, p1 = initial guesses
  ##   f = function to find the root of
  ##   eps = tolerance (defaults to 1e-9 = .000000001)
  ##   n = iterations allowed (defaults to 30)
  
  Q0 = f(p0)
  Q1 = f(p1)
  
  failureMode=TRUE
  labelName = c("current_p", "p_0", "p_1", "err (p0-p1)")
  perIteration = matrix(0, 0, length(labelName))  #initialize a matrix to save iterations
  for (i in 1:n) {
    p = p1 - Q1*(p1-p0)/(Q1-Q0)  #secant step
    perIteration = rbind(perIteration,c(p, p0, p1, abs(p0-p1))) #save stuff
    if ( f(p)==0 || abs(p-p1) < eps ) { failureMode=FALSE; break }
    p0=p1; Q0=Q1; #move p1 to p0
    p1=p; Q1=f(p1) #save new p1
  }
  if (failureMode) warning(paste("Failed to converge after",i,"iterations"))
  if (iter) { dimnames(perIteration)=list(1:i, labelName); p=perIteration }
  return(p)
}

############################################################################
## Testing the code
## Function to find the zero of
f = function(x) { x^2-3 }

## Examples
secant(1,2,f)

bb=secant(1,2,f,iter=TRUE)
bb
noquote(formatC(bb,digits=11,format="f")) #error is formatted differently

secant(1,2,f,1e-15)
secant(1,2,f,1e-15,100)
secant(1,2,f,1e-15,100,iter=TRUE)
############################################################################
