########################################################################
#
# Bisection Method
#   Implements the bisection method for finding a root
#
########################################################################
bisection = function(a,b,f,eps=1e-7,n=30,iter=FALSE) {
  ## Inputs
  ##   [a,b] = interval root is bracketed by
  ##   f = function to find the root of
  ##   eps = tolerance (defaults to 1e-9 = .000000001)
  ##   n = iterations allowed (defaults to 30)
  
  if ( sign(f(a)) * sign(f(b)) > 0 )
    stop(paste("root does not exist in [", a, ",", b, "]", sep=""))
  
  failureMode = TRUE
  labelName = c("Midpt (p)","LeftB (a)","RightB (b)","err (b-a)")
  perIteration = matrix(0, 0, length(labelName))  #initialize a matrix to save iterations
 
   for (i in 1:n) {
    p = (a+b)/2
    #  rbind() function combines vector, matrix or data frame by rows.
    perIteration = rbind(perIteration,c(p,a,b,abs(b-a)/2)) #save stuff
    # writeLines("Per Iteration: ")
    # print( paste(perIteration) )
    if ( f(p)==0 || (b-a)/2 < eps ) { failureMode=FALSE;  break  }
    if ( sign(f(a))*sign(f(p)) > 0 ) { a=p } else { b=p }
   }
  
  if (failureMode){
    warning( paste("Failed to converge after",i,"iterations") )
  }
  
  if (iter) { 
    # The dimnames()set both rows and columns names of a matrix.
    dimnames(perIteration)=list(1:i,labelName); 
    p=perIteration # swap the answer to per iteration matrix
  }
  return(p) # bring out the answer
}
############################################################################
## Testing the code
## Function to find the zero of
# f = function(x) { x^2-3 }

## Examples:
# bisection(1, 2, f)

# bb=bisection(1,2,f,iter=TRUE)
# bb
# noquote(formatC(bb, digits=11, format="f")) #error is formatted differently

# bisection(1,2,f,1e-15)
# bisection(1,2,f,1e-15,100)
# bisection(1,2,f,1e-15,100,iter=TRUE)
############################################################################