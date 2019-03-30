############################################################################
#
# Hermite Interpolation -
#   This returns the coefficients of the Hermite interpolatory polynomial H
#   where f is stored in table form as A, which contains the
#     x-values in the first column, the f(x-values) in the second column,
#     and the f'(x-values) in the third column
#  written by Scott Hyde
############################################################################

hermite = function(A, xs) {
  ##
  ## Inputs
  ##   A = table of nodes with function values and deriv values
  ##   xs = value to interpolate from the table
  ##
  ## Note that n below is really n+1 from algorithm
  ##  This is because R does not use 0 as an index, so
  ##  everything has to be increased by one (except when
  ##  subtracting i-j, which needs to be increased by 1
  ##  as well
  n = dim(A)[1]
  x = A[, 1]
  z = rep(x, each = 2)
  Q = matrix(NA, 2 * n, 2 * n)
  Q[, 1] = rep(A[, 2], each = 2)
  Q[seq(2, 2 * n, by = 2), 2] = A[, 3]
  Q[seq(3, 2 * n - 1, by = 2), 2] = diff(A[, 2]) / diff(x)
  
  for (i in 3:(2 * n)) {
    for (j in 3:i) {
      Q[i, j] = (Q[i, j - 1] - Q[i - 1, j - 1]) / (z[i] - z[i - j + 1])
    }
  }
  
  #  The below code follows the algorithm in the book explicitly:
  #	n=dim(A)[1]-1
  #	z=rep(0,2*n+1)
  #	Q=matrix(NA,2*n+2,2*n+2)
  #
  #	## note how the code adds one to EACH location value
  #	for (i in 0:n) {
  #		row=i+1
  #		evnrow=2*i+1  #even row + 1 (no zero row or column)
  #	    oddrow=2*i+1+1 #odd row
  #		z[evnrow]=A[row,1]
  #		z[oddrow]=A[row,1]
  #		Q[evnrow,0+1]=A[row,2]
  #		Q[oddrow,0+1]=A[row,2]
  #		Q[oddrow,1+1]=A[row,3]
  #		if ( i != 0 ) Q[evnrow,1+1]=(Q[evnrow,0+1]-Q[oddrow-2,0+1])/(z[evnrow]-z[oddrow-2])
  #	}
  #
  #	for (i in 2:(2*n+1) ) {
  #		for (j in 2:i) {
  #			# note EACH location has +1 added to it Q[i,j] becomes Q[i+1,j+1] because
  #			#   R cannot do matrices with 0 indices.
  #			Q[i+1,j+1]=(Q[i+1,j-1+1]-Q[i-1+1,j-1+1])/(z[i+1]-z[i-j+1])
  #		}
  #	}
  
  # The coefficients of the Newton Interpolary Divided Difference
  # Formula are the diagonal entries of Q
  coef = diag(Q)
  # The next two lines use the coefficients to figure out the
  # interpolation.  First line creates a vector of the product of
  # x-xj, then the second finds the dot product of them.
  zvec = c(1, cumprod(xs - z[-length(z)]))
  interp = sum(coef * zvec)
  
  return(list(
    table = Q,
    coef = coef,
    interp = interp
  ))
}

# for my homework testing
# f = function(x){
#   h = 1e-8
#   herm=rep(0,length(x))
#   for (i in 1:length(x)) {
#     herm[i]=(hermite(A,x[i]+h)$interp-hermite(A,x[i])$interp)/h-242/3
#   }
#   return(herm)
# }
# firsttime=bisection(5,6,f,eps=1e-12,n=100)
# firsttime

