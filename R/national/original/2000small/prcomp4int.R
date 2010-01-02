############################################################################
# Principal components for each state...let's add one orthogonal component #
############################################################################

#1999 estimates
#per capita race

library(mva);

#This file needs to be run from a precise location in natplot.tanh2.R

#let's get the orthogonal data.



ndemographics _ 9;
demographics _ as.data.frame(matrix(0,nrow=obs,ncol=ndemographics));
nnames _ c("inc","cen00.native","cen00.asian","college","cen00.black","cen00.hispanic","cen00.white","densityUP","cen00.population");
names(demographics) _ nnames;

density _ us.tmp$density/us.tmp$pop*us.tmp$cen00.population;

a1 _ lm(us.tmp$inc ~ int1 + int2 + int3 + prep96 + pperot96);
demographics$inc _ a1$residuals;

a1 _ lm(us.tmp$cen00.native/us.tmp$cen00.population ~ int1 + int2 + int3 + prep96 + pperot96);
demographics$cen00.native _ a1$residuals;

a1 _ lm(us.tmp$cen00.asian/us.tmp$cen00.population ~ int1 + int2 + int3 + prep96 + pperot96);
demographics$cen00.asian _ a1$residuals;

a1 _ lm(us.tmp$college ~ int1 + int2 + int3 + prep96 + pperot96);
demographics$college _ a1$residuals;

a1 _ lm(us.tmp$cen00.black/us.tmp$cen00.population ~ int1 + int2 + int3 + prep96 + pperot96);
demographics$cen00.black _ a1$residuals;

a1 _ lm(us.tmp$cen00.hispanic/us.tmp$cen00.population ~ int1 + int2 + int3 + prep96 + pperot96);
demographics$cen00.hispanic _ a1$residuals;

a1 _ lm(us.tmp$cen00.white/us.tmp$cen00.population ~ int1 + int2 + int3 + prep96 + pperot96);
demographics$cen00.white _ a1$residuals;

a1 _ lm(density ~ int1 + int2 + int3 + prep96 + pperot96);
demographics$densityUP _ a1$residuals;

a1 _ lm(us.tmp$cen00.population ~ int1 + int2 + int3 + prep96 + pperot96);
demographics$cen00.population _ a1$residuals;

prcomp.results  _ prcomp(demographics, scale=TRUE);

rm(density);
rm(demographics);

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
