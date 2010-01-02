dat <- as.data.frame(read.csv("table.out",header=F))
colnames(dat) <- c("precinct","pres","senate","count")

dat$count[is.na(dat$count)] <- 0
aidx <- substr(as.character(dat$precinct),1,1) == "a"

dat$isabs _ ifelse(aidx,1,0)

xtabs(count ~ pres + senate, data=dat)
xtabs(count ~ pres + senate + isabs, data=dat)

dat$pspoiled _ dat$pres=="over" | dat$pres=="under" ;
dat$sspoiled _ dat$senate=="over" | dat$senate=="under" ;

dat$vbuchanan _ factor(ifelse(dat$pres=="buchanan", "ybuch", "nbuch"))
dat$vnelson   _ factor(ifelse(dat$senate=="nelson", "ynels", "nnels"))
dat$vdeckard  _ factor(ifelse(dat$senate=="deckard", "ydeck", "ndeck"))
dat$vsenate   _
  factor(ifelse(dat$senate=="nelson", "nelson",
    ifelse(dat$senate=="deckard", "deckard", "other")))

dat$igore _ ifelse(dat$pres=="gore", 1,0)
dat$ibuchanan _ ifelse(dat$pres=="buchanan", 1,0)
dat$inelson   _ ifelse(dat$senate=="nelson", 1,0)
dat$ideckard  _ ifelse(dat$senate=="deckard", 1,0)

# show 3x2 tables

print("all ballots")
xtabs(count ~ igore + inelson, data=dat);

print("election-day ballots")
etab _ xtabs(count ~ igore + inelson, subset= isabs==0, data=dat);
etab;
etab/apply(etab,1,sum);

print("absentee ballots")
atab _ xtabs(count ~ igore + inelson, subset= isabs==1, data=dat);
atab;
atab/apply(atab,1,sum);

print("election-day ballots, excluding under/overvoted presidential ballots")
etabns _ xtabs(count ~ igore + inelson, subset= isabs==0 & !pspoiled, data=dat);
etabns;
etabns/apply(etabns,1,sum);

print("absentee ballots, excluding under/overvoted presidential ballots")
atabns _ xtabs(count ~ igore + inelson, subset= isabs==1 & !pspoiled, data=dat);
atabns;
atabns/apply(atabns,1,sum);

