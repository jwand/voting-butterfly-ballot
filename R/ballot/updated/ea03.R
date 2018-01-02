dat <- as.data.frame(read.csv("../data/table.out",header=F))
colnames(dat) <- c("precinct","pres","senate","count")

dat$count[is.na(dat$count)] <- 0
aidx <- substr(as.character(dat$precinct),1,1) == "a"

dat$isabs <- ifelse(aidx,1,0)

xtabs(count ~ pres + senate, data=dat)
xtabs(count ~ pres + senate + isabs, data=dat)

dat$pspoiled <- dat$pres=="over" | dat$pres=="under" ;
dat$sspoiled <- dat$senate=="over" | dat$senate=="under" ;

dat$vbuchanan <- factor(ifelse(dat$pres=="buchanan", "ybuch", "nbuch"))
dat$vnelson   <- factor(ifelse(dat$senate=="nelson", "ynels", "nnels"))
dat$vdeckard  <- factor(ifelse(dat$senate=="deckard", "ydeck", "ndeck"))
dat$vsenate   <-
  factor(ifelse(dat$senate=="nelson", "nelson",
    ifelse(dat$senate=="deckard", "deckard", "other")))

dat$ibuchanan <- ifelse(dat$pres=="buchanan", 1,0)
dat$inelson   <- ifelse(dat$senate=="nelson", 1,0)
dat$ideckard  <- ifelse(dat$senate=="deckard", 1,0)

# show 3x2 tables

print("all ballots")
xtabs(count ~ vsenate + vbuchanan, data=dat);

print("election-day ballots")
etab <- xtabs(count ~ vsenate + vbuchanan, subset= isabs==0, data=dat);
etab;
etab/apply(etab,1,sum);

print("absentee ballots")
atab <- xtabs(count ~ vsenate + vbuchanan, subset= isabs==1, data=dat);
atab;
atab/apply(atab,1,sum);

print("election-day ballots, excluding under/overvoted presidential ballots")
etabns <- xtabs(count ~ vsenate + vbuchanan, subset= isabs==0 & !pspoiled, data=dat);
etabns;
etabns/apply(etabns,1,sum);

print("absentee ballots, excluding under/overvoted presidential ballots")
atabns <- xtabs(count ~ vsenate + vbuchanan, subset= isabs==1 & !pspoiled, data=dat);
atabns;
atabns/apply(atabns,1,sum);

# build individual-level dataframes

print("all ballots");
ndatrows <- dim(dat)[1]
nobs <- sum(dat$count)
wrkdat <- matrix(0,4,nobs)
j <- 0
for (i in 1:ndatrows) {
  if (dat$count[i] > 0) {
    wrkdat[1:4,j+(1:dat$count[i])] <-
      c(dat$ibuchanan[i], dat$inelson[i], dat$ideckard[i], dat$isabs[i]);
    j <- j + dat$count[i]
  }
}
dimnames(wrkdat) <- list(c("ibuchanan","inelson","ideckard","isabs"), NULL)
wrkframe <- data.frame(t(wrkdat))

print("election-day:  all ballots");
summary(glm(ibuchanan ~ inelson + ideckard, family="binomial",
 subset= isabs==0, data=wrkframe ))

print("absentee:  all ballots");
summary(glm(ibuchanan ~ inelson + ideckard, family="binomial",
 subset= isabs==1, data=wrkframe ))


print("excluding uder/overvoted  ballots");
ndatrows <- dim(dat)[1]
nobs <- sum(dat$count[!dat$pspoiled])
wrkdat <- matrix(0,4,nobs)
j <- 0
for (i in 1:ndatrows) {
  if (!dat$pspoiled[i] && dat$count[i] > 0) {
    wrkdat[1:4,j+(1:dat$count[i])] <-
      c(dat$ibuchanan[i], dat$inelson[i], dat$ideckard[i], dat$isabs[i]);
    j <- j + dat$count[i]
  }
}
dimnames(wrkdat) <- list(c("ibuchanan","inelson","ideckard","isabs"), NULL)
wrkframe <- data.frame(t(wrkdat))

print("election-day:  excluding under/overvoted presidential ballots");
summary(glm(ibuchanan ~ inelson + ideckard, family="binomial",
 subset= isabs==0, data=wrkframe ))

print("absentee:  excluding under/overvoted presidential ballots");
summary(glm(ibuchanan ~ inelson + ideckard, family="binomial",
 subset= isabs==1, data=wrkframe ))
