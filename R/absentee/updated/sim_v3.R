
## things we need
## nvoters:   number of voters by mail, at polls
## ndistrict: number of voting units
## p: probability of voting for buch (assume equal for all)
##    p.buch
## nvote.buch.mail
set.seed(391)

nvoters.mail.small <- 100*10
nvoters.poll.small <- 1000*10
nvoters.mail.big <- 100*50
nvoters.poll.big <- 1000*50

ndistricts <- 10000
p.buch <- .1
p.vote <- .8

f.simulate.vote <- function(ndistricts, nvoters,  pvote, pbuch) {
  nvote  <- rbinom( ndistricts, nvoters,  pvote)
  nbuch  <- rbinom( ndistricts, nvote  ,  pbuch)
  prop.buch <- nbuch/nvote
  return( list( nvote=nvote, nbuch=nbuch, prop.buch=prop.buch) )
}

mail.small <- f.simulate.vote( ndistricts, nvoters.mail.small, p.vote, p.buch)
poll.small <- f.simulate.vote( ndistricts, nvoters.mail.small, p.vote, p.buch)

mail.big <- f.simulate.vote( ndistricts, nvoters.mail.big, p.vote, p.buch)
poll.big <- f.simulate.vote( ndistricts, nvoters.mail.big, p.vote, p.buch)

summary(poll.small$prop.buch)
summary(mail.small$prop.buch)


d.small <- sort(poll.small$prop.buch - mail.small$prop.buch)
d.big <- sort(poll.big$prop.buch - mail.big$prop.buch)
plot( d.small , cex=.1,frame=F,main="Comparing distributions of differences")
points( d.big, cex=.1, col="blue")

y <- 1:ndistricts / ndistricts
plot( d.small , y, cex=.1,frame=F,main="Empirical CDF of differences",xlim=c(-.06,.06),ylab="Fraction")
points( d.big , y, cex=.1, col="blue")

hist( d.small, freq=TRUE) 
hist( d.small, freq=FALSE) 
lines( density(poll.small$prop.buch - mail.small$prop.buch),ylim=c(0,70),frame=F,main="Comparing distributions of differences")

d.big <- density(poll.big$prop.buch - mail.big$prop.buch)
d.small <- d <- density(poll.small$prop.buch - mail.small$prop.buch)
idx <- d$x > -0.035 & d$x < 0.035
xx <- c(d$x[idx], rev(d$x[idx]))
yy <- c(rep(0,sum(idx)), rev(d$y[idx]))
polygon( xx,yy, col="green")

plot( d.small,ylim=c(0,70),frame=F,main="Comparing distributions of differences")
polygon( xx,yy, col="green")
lines( d.big,col="blue")

