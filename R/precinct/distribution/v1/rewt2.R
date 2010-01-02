source("../rewtanh2.R")


# tmp <-
# cbind(ss35.buch.100.1$w,ss35.buch.100.1$pop,ss35.buch.100.1$dep*ss35.buch.100.1$pop,ss35.buch.100.1$indep)
# write.table(as.data.frame(tmp),
#            sep=",",row.names=F,col.names=F,file="ss35.buch.100.1.dat")

rest2 <- rewtanh(ss35.buch.100.1$par,ss35.buch.100.1$s0^2,
  ss35.buch.100.1$stand.res,
  ss35.buch.100.1$pop,
  ss35.buch.100.1$dep*ss35.buch.100.1$pop,
  ss35.buch.100.1$indep)

summary(rest2$glmobj, dispersion=rest2$disp)
print("sandwich SEs")
print(sqrt(diag(rest2$rcov)))

print(rest2)

rest2$sres/sqrt(1-rest2$h)

cbind(ss35.buch.100.1$w, rest2$w)

round(cbind(1:105,ss35.buch.100.1$w, rest2$w,
            ss35.buch.100.1$raw.res,ss35.buch.100.1$student.res,
            ss35.buch.100.1$dep,ss35.buch.100.1$pred,
            ss35.buch.100.1$dep-ss35.buch.100.1$pred)[ ss35.buch.100.1$w != rest2$w ,],4)
