############################################################################
#
# Cubic Spline Interpolation -
#   This returns the coefficients, aj,bj,cj,dj for i=0,...,n-1 for the
#     cubic spline interpolatry function
#          S_j = aj+bj(x-xj)+cj(x-xj)^2+dj(x-xj)^3
#   given f is stored in table form as A, which contains the 
#     x-values in the first column, the f(x-values) in the second
#     column,
#   The return is a matrix with the columns [ xj aj bj cj dj ]
#  written by Scott Hyde from algorithm 3.4&3.5
#    additions now plot the function with splineplot below
############################################################################

spline = function(A,type="free",fp0,fpn) {
  ##
  ## Inputs
  ##   A = table of nodes with function values and deriv values
  ##   type = "free" or "clamped".  Free is default.
  ##   fp0 = f'(x0)  only used with clamped
  ##   fpn = f'(xn)  only used with clamped
  ## For the code below, the vector (x0,x1,...,xn) has been
  ## reindexed to start at 1 - (x1,x2,...,x(n+1)).  The index of the
  ## for loops have been modified from the original algorithm.
  x=A[,1]
  a=A[,2]
  n=dim(A)[1]-1
  
  ## This will do Step one of the algorithm
  h=diff(x)
  if ( type == "clamped" ) {
    alpha= 3*(a[2]-a[1])/h[1]-3*fp0; #first alpha
    alpha[n+1]= 3*fpn - 3*(a[n+1]-a[n])/h[n]
  } else alpha=rep(0,n-1) # Set alpha[1]=alpha[n+1]=0 too
  for (i in 2:n) alpha[i] = 3/h[i]*(a[i+1]-a[i])-3/h[i-1]*(a[i]-a[i-1]) 
  if ( type == "clamped" ) {
    l=2*h[1]; mu=0.5; z=alpha[1]/l[1]
  } else {
    l=2; mu=0; z=0;  # sets the first element
  }
  for (i in 2:n) {
    l[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*mu[i-1]
    mu[i] = h[i]/l[i]
    z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i]
  }
  c=b=d=rep(0,n+1)
  if ( type == "clamped" ) {
    l[n+1] = h[n]*(2-mu[n]);
    z[n+1] = (alpha[n+1]-h[n]*z[n])/l[n+1];
    c[n+1] = z[n+1]
  } else {
    l[n+1] = 1; z[n+1] = 0;
  }
  for (j in n:1 ) {
    c[j] = z[j] - mu[j]*c[j+1];
    b[j] = (a[j+1]-a[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3;
    d[j] = (c[j+1]-c[j])/(3*h[j])
  }
  mat = cbind(x, a, b, c, d)#[-(n+1),] #drop last row
  mat[n + 1, 3:5] = NA
  return(mat)
}

splinefunc = function(M,x) {
  ## This is the spline function defined by the matrix M (the return
  ## of the spline function.)
  y=c()
  for ( i in 1:length(x) ) {
    whichf=if (all(M[,1]>=x[i])) 1 else max(which(M[,1]<=x[i]))
    y[i] = sum(M[whichf,]*c(0,1,cumprod(rep(x[i]-M[whichf,1],3))))
  }
  return(y)
}

splineplot = function(M,xmin=min(M[,1]),xmax=max(M[,1])) {
  ## Use this to get a plot of the spline over the range defined by
  ## the matrix.  If you want to go beyond the limits, then choose a
  ## new xmin and xmax
  par(pty="s")
  xx=seq(xmin,xmax,length=1000)
  plot(xx,splinefunc(M,xx),type="l",xlab="x",ylab="spline",asp=1)
  
}

