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

xtabs(count ~ vnelson + vbuchanan, subset=isabs==0, data=dat)
xtabs(count ~ vnelson + vbuchanan, subset=isabs==1, data=dat)

chisq.test(xtabs(count ~ vnelson + vbuchanan, subset=isabs==0, data=dat))
chisq.test(xtabs(count ~ vnelson + vbuchanan, subset=isabs==1, data=dat))

xtabs(count ~ vnelson + vbuchanan + isabs, data=dat)
ftable(xtabs(count ~ isabs + vnelson + vbuchanan, data=dat))

xtabs(count ~ vnelson + vbuchanan, subset= !pspoiled & isabs==0, data=dat)
xtabs(count ~ vnelson + vbuchanan, subset= !pspoiled & isabs==1, data=dat)

chisq.test(xtabs(count ~ vnelson + vbuchanan, subset= !pspoiled & isabs==0, data=dat))
chisq.test(xtabs(count ~ vnelson + vbuchanan, subset= !pspoiled & isabs==1, data=dat))

xtabs(count ~ vnelson + vbuchanan + isabs, subset= !pspoiled, data=dat)
ftable(xtabs(count ~ isabs + vnelson + vbuchanan, subset= !pspoiled, data=dat))

attach(dat)

ibuchanan _ ifelse(dat$pres=="buchanan", 1,0)
inelson   _ ifelse(dat$senate=="nelson", 1,0)
ideckard  _ ifelse(dat$senate=="deckard", 1,0)

# election-day:  all ballots
summary(glm(ibuchanan ~ inelson + ideckard, weight = count, family=quasibinomial,
 subset= count>0 & isabs==0 ))

# absentee:  all ballots
summary(glm(ibuchanan ~ inelson + ideckard, weight = count, family=quasibinomial,
 subset= count>0 & isabs==1 ))

# election-day:  excluding under/overvoted presidential ballots
summary(glm(ibuchanan ~ inelson + ideckard, weight = count, family=quasibinomial,
 subset= count>0 & !pspoiled & isabs==0 ))

# absentee:  excluding under/overvoted presidential ballots
summary(glm(ibuchanan ~ inelson + ideckard, weight = count, family=quasibinomial,
 subset= count>0 & !pspoiled & isabs==1 ))
