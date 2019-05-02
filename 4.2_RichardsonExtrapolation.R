############################################################################
# this is formula works for the even power of h.
# which is O^(2j)
############################################################################
richardson = function (x0, h, f, n) {
  options(digits=12)
  
  # Three Point Central Difference (Eq. 4.5)
  df=function(x,h) { (f(x+h) - f(x-h))/(2*h)}
  
  tabl=matrix (, n, n)
  tabl[1,1] = df(x0,h)
  h=h/2
  
  for (i in 2: n) {
    tabl[i,1] = df(x0, h)
    for (j in 2:i)
      tabl[i,j]=tabl[i,j-1]+(tabl[i,j-1]-tabl[i-1, j-1])/(4^(j-1) -1)
    h=h/2
  }
  return(tabl)
}

# ************************************
f=function(x) {x + exp(x)}
richardson(0, 0.4, f, 3)
# [,1]         [,2]          [,3]
# [1,] 2.02688081451           NA            NA
# [2,] 2.00668001271 1.9999464121            NA
# [3,] 2.00166750020 1.9999966627 2.00000001274

# ************************************
