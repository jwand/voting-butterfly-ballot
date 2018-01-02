# crunch the Florida 2000 county-level data on absentee and election-day votes for
# Buchanan as a test of the equality of the proportions for each ballot type

abseday <- read.csv(file="../data/abseday.FL2000.csv", header=TRUE);
dim(abseday)

names(abseday)
# [1] "X."                "County"            "Total.Pres.Votes" 
# [4] "Total.Buchanan"    "Total.Absentee"    "Buchanan.Absentee"
# [7] "Elecday.Buch.."    "Absentee.Buch.."   "Abs.Elec.Ratio"   
#[10] "Difference"       

attach(abseday);
totpres  <- Total.Pres.Votes;
totbuch  <- Total.Buchanan;
totabs   <- Total.Absentee;
toteday  <- Total.Pres.Votes - Total.Absentee;
absbuch  <- Buchanan.Absentee;
edaybuch <- Total.Buchanan - Buchanan.Absentee;
detach(abseday)


pdiff   <- edaybuch/toteday-absbuch/totabs

pdf(file="abseday2.pdf",
           width=6,height=6,onefile=F)

plot(totpres,pdiff, type="n", cex=.9, axes=F,
  ylab="Election-Day Minus Absentee Proportion for Buchanan",
  xlab="Total Number of Presidential Votes Cast")
axis(1, at=c(0,100000,200000,300000,400000,500000,600000),
     labels=c("0","100,000","200,000","300,000","400,000","500,000","600,000"))
axis(2)
points(totpres,edaybuch/toteday-absbuch/totabs, cex=.6)
lines(seq(1,max(totpres),500),1/sqrt(seq(1,max(totpres),500)), lty=2)
lines(seq(1,max(totpres),500),-1/sqrt(seq(1,max(totpres),500)), lty=2)
text(totpres[50],pdiff[50]+.0030,"Palm Beach", cex=1.1)
arrows(totpres[50],pdiff[50]+.0025-.0005,totpres[50],pdiff[50]+.00075, length=.075, code=2)

dev.off()
