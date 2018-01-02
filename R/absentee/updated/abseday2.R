# crunch the Florida 2000 county-level data on absentee and election-day votes for
# Buchanan as a test of the equality of the proportions for each ballot type

abseday <- read.csv(file="../data/abseday.FL2000.csv", header=TRUE);

attach(abseday);

#> names(abseday)
# [1] "X."                "County"            "Total.Pres.Votes" 
# [4] "Total.Buchanan"    "Total.Absentee"    "Buchanan.Absentee"
# [7] "Elecday.Buch.."    "Absentee.Buch.."   "Abs.Elec.Ratio"   
#[10] "Difference"       

totpres  <- Total.Pres.Votes;
totbuch  <- Total.Buchanan;
totabs   <- Total.Absentee;
toteday  <- Total.Pres.Votes - Total.Absentee;
absbuch  <- Buchanan.Absentee;
edaybuch <- Total.Buchanan - Buchanan.Absentee;
detach(abseday)

# Elecday.Buch.. == (100 * edaybuch/toteday)
# Absentee.Buch.. == (100 * absbuch/totabs)
# Abs.Elec.Ratio == (100 * (absbuch/totabs)/(edaybuch/toteday))
# Difference == (100 * ((edaybuch/toteday)-(absbuch/totabs)))

sigma <- 3.81;
pbucha <- absbuch / totabs;
pbuche <- edaybuch / toteday;
pvara <- pbucha * (1-pbucha) / totabs;
pvare <- pbuche * (1-pbuche) / toteday;
pvar <- pvara + pvare;
eapdiff  <- pbuche-pbucha;
zeapdiff <- eapdiff/(sigma*sqrt(pvar));

diffmat <- cbind(100*eapdiff, zeapdiff);
dimnames(diffmat) <- list(as.character(abseday$County), c("pctdiff","Zdiff"));
print("Election-day minus Absentee vote for Buchanan under HO:  equal proportions");
print("Simple percentage difference (col 1), z-score of difference (col 2)");
print("rows ordered by county name");
print(diffmat);

print("Election-day minus Absentee vote for Buchanan under HO:  equal proportions");
print("Simple percentage difference (col 1), z-score difference (col 2)");
print("rows ordered by absolute percentage difference");
print(diffmat[rev(order(diffmat[,1])),]);

print("Election-day minus Absentee vote for Buchanan under HO:  equal proportions");
print("Simple percentage difference (col 1), z-score difference (col 2)");
print("rows ordered by absolute z-score difference");
print(diffmat[rev(order(diffmat[,2])),]);
