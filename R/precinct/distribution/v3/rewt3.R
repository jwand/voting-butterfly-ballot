source("../rewtanh2.R")


# tmp <-
# cbind(us16.buch.100.1$w,us16.buch.100.1$pop,us16.buch.100.1$dep*us16.buch.100.1$pop,us16.buch.100.1$indep)
# write.table(as.data.frame(tmp),
#            sep=",",row.names=F,col.names=F,file="us16.buch.100.1.dat")

rest3 <- rewtanh(us16.buch.100.1$par,us16.buch.100.1$s0^2,
  us16.buch.100.1$stand.res,
  us16.buch.100.1$pop,
  us16.buch.100.1$dep*us16.buch.100.1$pop,
  us16.buch.100.1$indep)

summary(rest3$glmobj, dispersion=rest3$disp)
print("sandwich SEs")
print(sqrt(diag(rest3$rcov)))

print(rest3)

rest3$sres/sqrt(1-rest3$h)

cbind(us16.buch.100.1$w, rest3$w)

round(cbind(1:length(us16.buch.100.1$w),us16.buch.100.1$w, rest3$w,
            us16.buch.100.1$raw.res,us16.buch.100.1$student.res,
            us16.buch.100.1$dep,us16.buch.100.1$pred,
            us16.buch.100.1$dep-us16.buch.100.1$pred)[ us16.buch.100.1$w != rest3$w ,],4)
