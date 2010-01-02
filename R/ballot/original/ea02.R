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

dat$ibuchanan _ ifelse(dat$pres=="buchanan", 1,0)
dat$inelson   _ ifelse(dat$senate=="nelson", 1,0)
dat$ideckard  _ ifelse(dat$senate=="deckard", 1,0)

# dat6: reduced data frame including all ballots
dat6 _ as.data.frame.table(
 xtabs(count ~ precinct + ibuchanan + inelson + ideckard,
   data=dat))
dat6$isabs _ ifelse(substr(as.character(dat6$precinct),1,1) == "a",1,0)

xtabs(Freq ~ ibuchanan + inelson + ideckard, data=dat6)
xtabs(Freq ~ ibuchanan + inelson + ideckard, subset= isabs==0, data=dat6)
xtabs(Freq ~ ibuchanan + inelson + ideckard, subset= isabs==1, data=dat6)

# dat6ns: reduced data frame excluding under/overvoted ballots
dat6ns _ as.data.frame.table(
 xtabs(count ~ precinct + ibuchanan + inelson + ideckard,
   subset= !pspoiled, data=dat))
dat6ns$isabs _ ifelse(substr(as.character(dat6ns$precinct),1,1) == "a",1,0)

xtabs(Freq ~ ibuchanan + inelson + ideckard, data=dat6ns)
xtabs(Freq ~ ibuchanan + inelson + ideckard, subset= isabs==0, data=dat6ns)
xtabs(Freq ~ ibuchanan + inelson + ideckard, subset= isabs==1, data=dat6ns)


# election-day:  all ballots
summary(glm(ibuchanan ~ inelson + ideckard, weight = Freq, family=quasibinomial,
 subset= Freq>0 & isabs==0, data=dat6 ))

# absentee:  all ballots
summary(glm(ibuchanan ~ inelson + ideckard, weight = Freq, family=quasibinomial,
 subset= Freq>0 & isabs==1, data=dat6 ))


# election-day:  excluding under/overvoted presidential ballots
summary(glm(ibuchanan ~ inelson + ideckard, weight = Freq, family=quasibinomial,
 subset= Freq>0 & isabs==0, data=dat6ns ))

# absentee:  excluding under/overvoted presidential ballots
summary(glm(ibuchanan ~ inelson + ideckard, weight = Freq, family=quasibinomial,
 subset= Freq>0 & isabs==1, data=dat6ns ))
