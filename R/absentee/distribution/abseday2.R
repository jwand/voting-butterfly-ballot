# crunch the Florida 2000 county-level data on absentee and election-day votes for
# Buchanan as a test of the equality of the proportions for each ballot type

abseday <- read.csv(file="../data/abseday.FL2000.csv", header=TRUE);

attach(abseday);

#> names(abseday)
# [1] "X."                "County"            "Total.Pres.Votes" 
# [4] "Total.Buchanan"    "Total.Absentee"    "Buchanan.Absentee"
# [7] "Elecday.Buch.."    "Absentee.Buch.."   "Abs.Elec.Ratio"   
#[10] "Difference"       

totpres  _ Total.Pres.Votes;
totbuch  _ Total.Buchanan;
totabs   _ Total.Absentee;
toteday  _ Total.Pres.Votes - Total.Absentee;
absbuch  _ Buchanan.Absentee;
edaybuch _ Total.Buchanan - Buchanan.Absentee;

# Elecday.Buch.. == (100 * edaybuch/toteday)
# Absentee.Buch.. == (100 * absbuch/totabs)
# Abs.Elec.Ratio == (100 * (absbuch/totabs)/(edaybuch/toteday))
# Difference == (100 * ((edaybuch/toteday)-(absbuch/totabs)))

sigma _ 3.81;
pbucha _ absbuch / totabs;
pbuche _ edaybuch / toteday;
pvara _ pbucha * (1-pbucha) / totabs;
pvare _ pbuche * (1-pbuche) / toteday;
pvar _ pvara + pvare;
eapdiff  _ pbuche-pbucha;
zeapdiff _ eapdiff/(sigma*sqrt(pvar));

diffmat _ cbind(100*eapdiff, zeapdiff);
dimnames(diffmat) _ list(as.character(County), c("pctdiff","Zdiff"));
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
