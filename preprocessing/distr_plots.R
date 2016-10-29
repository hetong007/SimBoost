par(mfrow=c(1,3))

load('../data/davis_data.Rda')
hist(davis_triplet[,3], main="Distribution of values in dataset Davis", xlab="binding affinity in pKd", col='grey', border='white', cex.lab = 1.4)
abline(v=7.0, col="red")

load('../data/metz_data.Rda')
hist(metz_triplet[,3], main="Distribution of values in dataset Metz", xlab="binding affinity in pKi", col='grey', border='white', cex.lab = 1.4)
abline(v=7.6, col="red")

load('../data/kiba_data.Rda')
hist(kiba_triplet[,3], main="Distribution of values in dataset KIBA", xlab="binding affinity as KIBA score", col='grey', border='white', cex.lab = 1.4)
abline(v=12.1, col="red")

par(mfrow=c(1,1))

