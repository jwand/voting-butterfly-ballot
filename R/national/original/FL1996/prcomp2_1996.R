############################################################################
# Principal components for each state...let's add one orthogonal component #
############################################################################

#per capita race

library(mva);

#This file needs to be run from a precise location in natplot.tanh2.R

#let's get the orthogonal data.
#The data to use:
#inc: MEDIAN HOUSEHOLD MONEY INCOME 1989 DOL        
#abo: POPULATION BY RACE 1990 (1C) ABS
#asian: POPULATION BY RACE 1990 (1C) ABS
#college: POPULATION BY RACE 1990 (1C) ABS
#black: POPULATION BY RACE 1990 (1C) ABS
#hispanic: HISPANIC ORIGIN POPULATION PERCENT OF TOTAL POPULATION
#white: POPULATION BY RACE 1990 (1C) ABS
#density: 1992 POPULATION PER 1990 SQUARE MILES RTE
#pop: POPULATION 1992 ABS

ndemographics _ 9;
demographics _ as.data.frame(matrix(0,nrow=obs,ncol=ndemographics));
nnames _ c("inc","abo","asian","college","black","hispanic","white","density","pop");
names(demographics) _ nnames;
a1 _ lm(us.tmp$inc ~ prep92 + pperot92);
demographics$inc _ a1$residuals;
a1 _ lm(us.tmp$abo/us.tmp$pop ~ prep92 + pperot92);
demographics$abo _ a1$residuals;
a1 _ lm(us.tmp$asian/us.tmp$pop ~ prep92 + pperot92);
demographics$asian _ a1$residuals;
a1 _ lm(us.tmp$college ~ prep92 + pperot92);
demographics$college _ a1$residuals;
a1 _ lm(us.tmp$black/us.tmp$pop ~ prep92 + pperot92);
demographics$black _ a1$residuals;
a1 _ lm(us.tmp$hispanic ~ prep92 + pperot92);
demographics$hispanic _ a1$residuals;
a1 _ lm(us.tmp$white/us.tmp$pop ~ prep92 + pperot92);
demographics$white _ a1$residuals;
a1 _ lm(us.tmp$density ~ prep92 + pperot92);
demographics$density _ a1$residuals;
a1 _ lm(us.tmp$pop ~ prep92 + pperot92);
demographics$pop _ a1$residuals;

prcomp.results  _ prcomp(demographics, scale=TRUE);
cat("\nprcomp results:\n")
print(prcomp.results)
cat("\n");
print(summary(prcomp.results));
pr1 _ (prcomp.results$x)[,"PC1"];
pr2 _ (prcomp.results$x)[,"PC2"];
pr3 _ (prcomp.results$x)[,"PC3"];
us.tmp$pr1 _ pr1;
us.tmp$pr2 _ pr2;
us.tmp$pr3 _ pr3;
