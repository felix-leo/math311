############################################################################
#
# Neville's Method - this function interpolates the value f(xs),
#   where xs, where f is stored in table form as A, which contains the
#     x-values and f(x-values) in the first and second columns
#  written by Scott Hyde
############################################################################

neville = function(A, xs, nodes = 0:(dim(A)[1] - 1)) {
  ## Inputs
  ##   A = table of nodes with function values
  ##   xs = value to interpolate from the table
  ##   nodes = nodes to use for the approximation.  Should be
  ##        created using :, so 2:4, 1:3, 0:2 (Uses zero index, so
  ##        first node is zero.
  if (all(diff(nodes) == 1)) {
    # This checks that the nodes are in a sequence separated by 1.
    n = length(nodes)
    nodes = nodes + 1
    x = A[nodes, 1]
    Q = matrix(, length(nodes), length(nodes))
    Q[, 1] = A[nodes, 2]
    
    for (i in 2:(dim(Q)[1])) {
      for (j in 2:i) {
        Q[i, j] = ((xs - x[i - j + 1]) * Q[i, j - 1] - (xs - x[i]) * Q[i - 1, j -
                                                                         1]) / (x[i] - x[i - j + 1])
      }
    }
    return(list(table = Q, interp = Q[n, n]))
  } else{
    stop(
      "Nodes need to be consecutive.  Use 2:4, 1:3, 0:2, etc\n  (Use zero index, so the first node is zero)"
    )
  }
}

############################################################################
## Testing the code

## Find approximations of exp(x) over 1.0, 1.1, 1.2, 1.3, 1.4, 1.5

x = 1 + (0:5) / 10  #[1] 1.0 1.1 1.2 1.3 1.4 1.5
y = exp(x)
A = cbind(x, y)
#       x     y
# [1,] 1.0 2.718282
# [2,] 1.1 3.004166
# [3,] 1.2 3.320117
# [4,] 1.3 3.669297
# [5,] 1.4 4.055200
# [6,] 1.5 4.481689

## Find approx to exp(1.5) using all nodes or only 1:2 (linear)
neville(A, 1.5)     # Use all nodes
neville(A, 1.5, 1:2) # Use only nodes 1:2


## Find approx to Bessel function (first kind) of order 0
x = c(1, 1.3, 1.6, 1.9, 2.2)
y = besselJ(x, 0)
A = cbind(x, y)
## Find approx to J0(1.5) using all nodes
neville(A, 1.5)     #Use all nodes

## Find approx to function defined below
x = c(1.00, 1.05, 1.10, 1.15)
y = c(.1924, .2414, .2933, .3492)
A = cbind(x, y)
## What is f(1.09)?
neville(A, 1.09)


