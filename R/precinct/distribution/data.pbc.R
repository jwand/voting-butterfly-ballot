pbc.pres <- read.table (file="../../data/2000pres.parsed.txt",header=T,as.is=T)
pbc.pres <- pbc.pres[,-1]                # kill observation column

  ptotal <-
    pbc.pres$bushm +
      pbc.pres$gorem +
        pbc.pres$brownem +
          pbc.pres$naderm +
            pbc.pres$harrism +
              pbc.pres$hagelinm +
                pbc.pres$buchananm +
                  pbc.pres$mcreynoldsm +
                    pbc.pres$phillipsm +
                      pbc.pres$mooreheadm
  
  pbuch <- pbc.pres$buchananm/ptotal



pbc.pres.total <- pbc.pres$bushm+pbc.pres$gorem+pbc.pres$brownem+pbc.pres$naderm+pbc.pres$harrism+pbc.pres$hagelinm+pbc.pres$buchananm+pbc.pres$mcreynoldsm+pbc.pres$phillipsm+pbc.pres$mooreheadm
pbc.pres[,-1] <- pbc.pres[,-1]/pbc.pres.total

  useme.ss35_rep(0,times=length(pbc.pres.total))
  useme.us16_rep(0,times=length(pbc.pres.total))  


pbc.senate <- read.table (file="../../data/2000senate.parsed.txt",header=T,as.is=T)
  stotal <-
    pbc.senate$mccollum +
      pbc.senate$nelson +
        pbc.senate$simonetta +
          pbc.senate$deckard +
            pbc.senate$logan +
              pbc.senate$martin +
                pbc.senate$mccormick +
                  pbc.senate$writein
  
  pdeckard <- pbc.senate$deckard / stotal
  pnelson <- pbc.senate$nelson / stotal

pbc.senate.total <- pbc.senate$mccollum+pbc.senate$nelson+pbc.senate$simonetta+pbc.senate$deckard+pbc.senate$logan+pbc.senate$martin+pbc.senate$mccormick+pbc.senate$writein

pbc.senate[,-1] <- pbc.senate[,-1]/pbc.senate.total

pbc.ss35 <- read.table (file="../../data/pbc.statesenate35.all.txt",header=T,as.is=T)
  useme.ss35[as.logical(pbc.ss35$vaughan>-99)]_1
  sstotal <- pbc.ss35$vaughan + pbc.ss35$rossin + pbc.ss35$lowe
  plowe <- pbc.ss35$lowe / sstotal


pbc.ss35.total <- pbc.ss35$vaughan+pbc.ss35$rossin+pbc.ss35$lowe

pbc.ss35[,-1] <- pbc.ss35[,-1]/pbc.ss35.total

##pbc.ss35[pbc.ss35 < 0,] <- 0




pbc.us16 <- read.table (file="../../data/pbc.ushouse16.all.txt",header=T,as.is=T)
  useme.us16[as.logical(pbc.us16$foley>-99)]_1
  shtotal <-
    pbc.us16$foley +
      pbc.us16$brown +
        pbc.us16$mcguire +
          pbc.us16$writein

  pmcguire <- pbc.us16$mcguire / shtotal


pbc.us16.total <- pbc.us16$foley+pbc.us16$brown+pbc.us16$mcguire+pbc.us16$writein

pbc.us16[,-1] <- pbc.us16[,-1]/pbc.us16.total


useme.ss35[532:637] <- 0
useme.us16[532:637] <- 0

useme.ss35[is.na(plowe)] <- 0
useme.us16[is.na(pmcguire)] <- 0



print(cbind(pbc.pres$precinct[useme.ss35==1],
      round(cbind(plowe,pbc.ss35$lowe,
                  plowe-pbc.ss35$lowe,
                  sstotal,
                  ptotal,
                  pbc.ss35.total,
                  pbc.pres.total,
                  useme.ss35)[useme.ss35==1,],4)))

print(cbind(pbc.pres$precinct[useme.us16==1],
            round(cbind(pmcguire,
                  pbc.us16$mcguire,
                  pmcguire-pbc.us16$mcguire,
                  shtotal,ptotal,
                  pbc.us16.total,
                  pbc.pres.total,
                  useme.us16)[useme.us16==1,],4)))


pbc <- as.data.frame(list(precinct=pbc.pres$precinct,
                          pbuch=pbuch,
                          ptotal=ptotal,
                          pdeckard=pdeckard,
                          pnelson=pnelson,
                          stotal=stotal,
                          plowe=plowe,
                          sstotal=sstotal,
                          pmcguire=pmcguire,
                          shtotal=shtotal,
                          useme.ss35=useme.ss35,
                          useme.us16=useme.us16))

pbc.complete <- cbind(pbc.pres,pbc.senate[,-1],pbc.ss35[,-1],pbc.us16[,-1],useme.ss35,useme.us16,pbc.pres.total)

#print(cbind(pbc.complete$precinct,pbc.complete$buchananm,pbc.complete$nelson,pbc.complete$deckard)[useme.ss35 == 1, ])

fit1 <- glm(buchananm~nelson+deckard,data=pbc.complete[useme.ss35==1,],family=binomial(link="logit"),weights=pbc.pres.total)
print (summary(fit1))


print(summary(glm(pbuch~pnelson+pdeckard,data=pbc[pbc$useme.ss35==1,],
          family=binomial(link="logit"),weights=ptotal)))


#print(cbind(pbc.complete$precinct,pbc.complete$lowe,pbc.complete$nelson,pbc.complete$deckard)[useme.ss35 == 1, ])

fit2 <- glm(lowe~nelson+deckard,data=pbc.complete[useme.ss35==1,],family=binomial(link="logit"),weights=pbc.pres.total)
print (summary(fit2))

print(summary(glm(plowe~pnelson+pdeckard,data=pbc[pbc$useme.ss35==1,],
          family=binomial(link="logit"),weights=sstotal)))


#print(cbind(pbc.complete$precinct,pbc.complete$buchananm,pbc.complete$nelson,pbc.complete$deckard)[useme.us16 == 1, ])

fit3 <- glm(buchananm~nelson+deckard,data=pbc.complete[useme.us16==1,],family=binomial(link="logit"),weights=pbc.pres.total)
print (summary(fit3))

print(summary(glm(pbuch~pnelson+pdeckard,data=pbc[pbc$useme.us16==1,],
          family=binomial(link="logit"),weights=ptotal)))


#print(cbind(pbc.complete$precinct,pbc.complete$mcguire,pbc.complete$nelson,pbc.complete$deckard)[useme.us16 == 1, ])

fit4 <- glm(mcguire~nelson+deckard,data=pbc.complete[useme.us16==1,],family=binomial(link="logit"),weights=pbc.pres.total)
print (summary(fit4))

print(summary(glm(pmcguire~pnelson+pdeckard,data=pbc[pbc$useme.us16==1,],
          family=binomial(link="logit"),weights=shtotal)))

