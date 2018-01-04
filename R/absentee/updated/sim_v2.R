
## things we need
## nvoters:   number of voters by mail, at polls
## ndistrict: number of voting units
## p: probability of voting for buch (assume equal for all)
##    p.buch
## nvote.buch.mail
set.seed(391)

nvoters.mail.small <- 100
nvoters.poll.small <- 1000
nvoters.mail.big <- 100*50
nvoters.poll.big <- 1000*50

ndistricts <- 1000
p.buch <- .1


f.simulate.vote <- function(ndistricts, nvoters,  p) {
  nvote  <- rbinom( ndistricts, nvoters,  p)
  pvote <- nvote/nvoters
  return( list( nvote=nvote, pvote=pvote) )
}

mail.small <- f.simulate.vote( ndistricts, nvoters.mail.small, p.buch)
poll.small <- f.simulate.vote( ndistricts, nvoters.mail.small, p.buch)

mail.big <- f.simulate.vote( ndistricts, nvoters.mail.big, p.buch)
poll.big <- f.simulate.vote( ndistricts, nvoters.mail.big, p.buch)

summary(poll.small$pvote)
summary(mail.small$pvote)

plot( density(poll.small$pvote - mail.small$pvote),ylim=c(0,70))
lines(density(poll.big$pvote - mail.big$pvote),)

d.small <- sort(poll.small$pvote - mail.small$pvote)
plot( jitter(d.small) , cex=.1)
lines(density(poll.big$pvote - mail.big$pvote),)



