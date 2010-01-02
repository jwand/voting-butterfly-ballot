source("../../lqd/rewtanh2.R")


# tmp <-
# cbind(pbc.buch.100.1$w,pbc.buch.100.1$pop,pbc.buch.100.1$dep*pbc.buch.100.1$pop,pbc.buch.100.1$indep)
# write.table(as.data.frame(tmp),
#            sep=",",row.names=F,col.names=F,file="pbc.buch.100.1.dat")

rest3 <- rewtanh(pbc.buch.100.1$par,pbc.buch.100.1$s0^2,
  pbc.buch.100.1$stand.res,
  pbc.buch.100.1$pop,
  pbc.buch.100.1$dep*pbc.buch.100.1$pop,
  pbc.buch.100.1$indep)

summary(rest3$glmobj, dispersion=rest3$disp)
print("sandwich SEs")
print(sqrt(diag(rest3$rcov)))

print(rest3)

cbind(pbc.buch.100.1$w, rest3$w)

round(cbind(1:length(pbc.buch.100.1$w),pbc.buch.100.1$w, rest3$w,
            pbc.buch.100.1$raw.res,pbc.buch.100.1$student.res,
            pbc.buch.100.1$dep,pbc.buch.100.1$pred,
            pbc.buch.100.1$dep-pbc.buch.100.1$pred)[ pbc.buch.100.1$w != rest3$w ,],4)
