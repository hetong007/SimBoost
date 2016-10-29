load('../data/davis_data.Rda')
load('../data/metz_data.Rda')
load('../data/kiba_data.Rda')

par(mfrow=c(1,3))

hist(davis_triplet[,3], main="Davis", xlab="pKd",
     cex.lab = 1.5, col='grey', border = 'white')
abline(v=7.0, col="red",lwd=3)

hist(metz_triplet[,3], main="Metz", xlab="pKd",
     cex.lab = 1.5, col='grey', border = 'white')
abline(v=7.6, col="red",lwd=3)

hist(kiba_triplet[,3], main="KIBA", xlab="KIBA score",
     cex.lab = 1.5, col='grey', border = 'white')
abline(v=12.1, col="red",lwd=3)

par(mfrow=c(1,1))