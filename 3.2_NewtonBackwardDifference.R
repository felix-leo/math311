############################################################################
#
# Newton's Forward Difference Formula & Backward Difference Formula
#   This returns the coefficients of the interpolatory polynomial P
#   where f is stored in table form as A, which contains the 
#     x-values and f(x-values) in the first and second columns 
#  written by Scott Hyde
############################################################################

nfdf = function(A,xs) {
  ##
  ## Inputs
  ##   A = table of nodes with function values
  ##   xs = value to interpolate from the table 
  ##
  n = dim(A)[1]-1  
  h = diff(A[1:2,1])  #find the difference between first two x's
  s = (xs-A[1,1])/h   #find s
  #check if all the differences are equal
  if (!all(diff(A[,1])-h < 1e-10)) stop("Nodes are not equally spaced")
  # Compute Deltaf for each k.
  deltaf = rep(0,n+1)
  for (k in 0:n) deltaf[k+1] = sum((-1)^(k-0:k)*choose(k,0:k)*A[0:k+1,2])
  # note that deltaf/(factorial(0:k)*h^(0:k)) are the values
  # returned in the divided differences table
  return(list(coef=deltaf,iterp=sum(choose(s,0:n)*deltaf)))
}

nbdf = function(A,xs) {
  ## This is easy!  Just apply the forward difference to the reverse
  ## of x and y
  return(nfdf(apply(A,2,rev),xs))
}

############################################################################
## Example
x=c(1,1.3,1.6,1.9,2.2)
y=besselJ(x,0)
A=cbind(x,y)

stuff = nfdf(A,1.5)  
stuff
stuff$table
stuff$coef
stuff$interp
############################################################################