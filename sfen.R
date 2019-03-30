########################################################################
#
# Steffensen's Method 
#   Implements the Steffensen method for finding a fixed point
#   This method is quadratically convergent without computing
#   derivatives at the expense of two function evaluations per step.
#
########################################################################

steff0 = function(p,g,eps=1e-7,n=50,iter=FALSE,complete=FALSE) {
  ##
  ## Inputs
  ##   p = initial guess
  ##   f = function to find the root of
  ##   eps = tolerance (defaults to 1e-9 = .000000001)
  ##   n = iterations allowed (defaults to 30)
  ##   complete = include fixed pt iterations too.
  ## Outputs
  ##   fixed point of g
  delta = function(p) { p[1]-(p[2]-p[1])^2/(p[3]-2*p[2]+p[1]) }
  fail=TRUE
  save=c("p_n","g(p)","abs(g(p)-p)")
  piter = matrix(c(p,g(p),abs(g(p)-p)),1,length(save))  #save iterations
  for (i in 1:n) {
    p1 = g(p); p2 = g(p1); #fixed point twice
    p = delta(c(p,p1,p2))  #Aitken's method once
    if (complete) piter = rbind(piter,c(p1,g(p1),abs(g(p1)-p1)),
                                c(p2,g(p2),abs(g(p2)-p2)))
    piter = rbind(piter,c(p,g(p),abs(g(p)-p))) #save stuff
    if ( abs(g(p)-p) < eps ) {fail=FALSE; break }
  }
  if (fail) warning(paste("Failed to converge after",i,"iterations"))
  if (iter) { dimnames(piter)=list(0:i,save);
    return(piter) }
  return(p)
}

############################################################################
## Testing the code
## Function to find the fixed point of
g = function(x) { cos(x) }

steff0(1,g)
bb=steff0(1,g,iter=TRUE)
bb
noquote(formatC(bb,digits=11,format="f")) #error is formatted differently

steff0(1,g,1e-15)
steff0(1,g,1e-15,100,iter=TRUE)
steff0(1,g,1e-15,100,iter=TRUE,complete=TRUE)
############################################################################
