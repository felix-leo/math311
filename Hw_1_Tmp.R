# Example of the Effects of Floating Point Error 
# on Numerical Computations
h=1
n=70
x0=1
expx0 = true = exp(x0)          #store exp(x0) into expx0 and true
hvec = rep(0,n)     #create a n dimentional vector
errvec = rep(0,n)
for (i in 1:n) {
  ddif = (exp(x0+h)-expx0)/h    #Compute divided difference approx.
  abs_error = abs(ddif - true)  #Compute abs. value of approx. error
  hvec[i] = h
  errvec[i] = abs_error
  h = h/2                       #decrease step size h
}

plot(log10(hvec),log10(errvec),type="l",xlab="log(h)",ylab="log(abs error(h))")

title("Divided Difference Approx Error vs. h ")
# **********************************
          # problem_B
# Part_a
par(mar=c(1,1,1,1))
#import(svglite)
library(MASS) # for fractions() function
library(ggplot2)
library(svglite)
xs=c(256, 1/1024, 0)
n=70
minh = c(0,0,0)
for(x0 in xs) {
  hvec = 2^ - ((1:n)-1)
  errvec = abs( (exp(x0+hvec) - exp(x0)) /hvec-exp(x0) )
  # which() function will return the postion of the elements
  # in a logical vector that is true
  name=paste("probset1_", which(xs == x0), ".svg", seq="")
  # svg is a type of a img file, like png, jpeg
  svglite::svglite(name, width=14, height=14)
  png(name, width=14, height=14)
  if (which(xs==x0)==3){
    plot(log10(hvec), errvec, type="1", xlab="log(h)", ylab="|error(h)|")
  }else{
    plot(log10(hvec), log10(errvec), type="1", xlab="log(h)", ylab="|error(h)|")
  }
  stx = toString(fractions(x0))
  # bquote() function quotes its argument 
  # except that terms wrapped in ' .() '
  # 
  title(bquote(paste("Divided Diff. Error of Deriv of ",
                     e^x[0], " vs. h",
                     "with", x[0], "=", 
                     .(stx)
                     )
               )
        )
  
  dev.off()
  #graohics.off()
  minh[which(xs==x0)] = hvec [max(which(errvec==min(errvec)))] #save min h
  }

# **********************************
          # problem_C
 xs = c(256,1/1024,0)
 min(which(2^(1:2000) == Inf)) -1;
 .Machine$double.exponent

# **********************************
          # problem_D
# Part_a
b = -1e9
c = -1
x_pos = (-b + sqrt(b^2- 4*c)) /2
x_pos

x_neg = (-b - sqrt(b^2 - 4*c )) /2
x_neg

# Part_b

