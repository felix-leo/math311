# Problem B

xs=c(256,1/1024,0)
n=70
minh=c(0,0,0)

graphics.off() # to avoid figure margins too large error

for(x0 in xs) {
  hvec = 2^-((1:n)-1) 
  errvec = abs((exp(x0+hvec)-exp(x0))/hvec-exp(x0))
  #name = paste("probset1_",which(xs==x0),".svg",sep="")
  #svg(name, width=14,height=14)
  
  if(which(xs==x0)==3)
    plot(log10(hvec),errvec,type="l",xlab="log(h)",ylab="|error(h)|")
  else
    plot(log10(hvec),log10(errvec),type = "l",xlab="log(h)",ylab="|error(h)|")
  stx = toString(fractions(x0))
  title(bquote(paste("Divided Diff. Error of Derv of ",
                     e^x[0],
                     " vs. h",
                     " with ",
                     x[0], "=",.(stx))))
  minh[which(xs==x0)]=hvec[max(which(errvec==min(errvec)))]
}
