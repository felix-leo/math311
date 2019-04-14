G = function(x) {
  z = (sin(x) - x + x ^ 3 / 6) / x ^ (1 / 4)
  
  z[is.na(z)] = 0
  z
}
options(digits = 12)

simp = simpsons(0, 1, G, n = 4, plot = true)

simp + 166/315


G2 = function(x){
  z = (exp(x) - (1+x +(x^2)/2 +(x^3)/6 +(x^4)/24) ) / sqrt(x)
  z[x == 0] = 0;
  z
}

simp = simpsons(0, 1, G2, n = 6, plot = true)
simp
exp(-1)*(simp + 11051/3780)

f = function(x) { 1/(9* x^2 +1) }
simpsons(0, 1, f, 4, plot=TRUE)
# ************************************************************

# 4.9 problem 4a
f4 = function(t){
  (1 -(t-1) +(t-1)^2 -(t-1)^3 +(t-1)^4 ) / ((t-1)^(3/4))
}
integrate(f4, 1, 2)
# 3.5720462487 with absolute error < 0.00024

G4 = function(t){
  z = (1/t -(1 -(t-1) +(t-1)^2 -(t-1)^3 +(t-1)^4 )) / ((t-1)^(3/4))
  z[t == 1] = 1;
  z
}

simpsons(1, 2, G4, 4, plot = true)
# [1] -0.0209950220411

(3.5720462487 - 0.0209950220411)/4
# [1] 0.887762806665

