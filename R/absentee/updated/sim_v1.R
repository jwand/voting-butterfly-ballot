
## things we need
## nvoters:   number of voters by mail, at polls
## ndistrict: number of voting units
## p: probability of voting for buch (assume equal for all)
##    p.buch
## nvote.buch.mail
set.seed(391)

nvoters.mail <- 100
nvoters.poll <- 1000
ndistricts <- 1000
p.buch <- .1

nvote.buch.mail <- rbinom( ndistricts, nvoters.mail,  p.buch)
nvote.buch.poll <- rbinom( ndistricts, nvoters.poll,  p.buch)
p.mail <- nvote.buch.mail/nvoters.mail
p.poll <- nvote.buch.poll/nvoters.poll
summary(p.mail)
summary(p.poll)

plot( density(p.mail - p.poll))

