par(mfrow = c(1,3))

# Davis

load('davis.inspection.rda')
load('../data/davis_data.Rda')

ord = order(Preds.upper-Preds.lower, decreasing = TRUE)
real = davis_triplet[ord,3]
ord.lower = Preds.lower[ord]
ord.upper = Preds.upper[ord]

triplet_ord = davis_triplet[ord,]

ind = which(triplet_ord[,2]==10)
real = real[ind]
ord.lower = ord.lower[ind]
ord.upper = ord.upper[ind]
ord.mid = (ord.lower + ord.upper)/2
triplet_ord = triplet_ord[ind,]
# ord = order(real)
ord = order(ord.mid)
real = real[ord]
ord.lower = ord.lower[ord]
ord.upper = ord.upper[ord]
ord.mid = ord.mid[ord]
triplet_ord = triplet_ord[ord,]

plot(real, type='l', ylim = range(c(real, ord.lower, ord.upper)),
     main = "Target 10 from Davis", xlab = "Drug Index", ylab = "Score")
lines(ord.lower, col=2, lty = 'dashed')
lines(ord.upper, col=3, lty = 'dotdash')
lines(ord.mid, col=4, lty='dotted')


# Metz

load('metz.inspection.rda')
load('../data/metz_data.Rda')

ord = order(Preds.upper-Preds.lower, decreasing = TRUE)
real = metz_triplet[ord,3]
ord.lower = Preds.lower[ord]
ord.upper = Preds.upper[ord]

triplet_ord = metz_triplet[ord,]

ind = which(triplet_ord[,2]==30)
real = real[ind]
ord.lower = ord.lower[ind]
ord.upper = ord.upper[ind]
ord.mid = (ord.lower + ord.upper)/2
triplet_ord = triplet_ord[ind,]
# ord = order(real)
ord = order(ord.mid)
real = real[ord]
ord.lower = ord.lower[ord]
ord.upper = ord.upper[ord]
ord.mid = ord.mid[ord]
triplet_ord = triplet_ord[ord,]

plot(real, type='l', ylim = range(c(real, ord.lower, ord.upper)),
     main = "Target 30 from Metz", xlab = "Index", ylab = "Score")
lines(ord.lower, col=2, lty = 'dashed')
lines(ord.upper, col=3, lty = 'dotdash')

lines(ord.mid, col=4, lty='dotted')

# Kiba

load('kiba.inspection.rda')
load('../data/kiba_data.Rda')

ord = order(Preds.upper-Preds.lower, decreasing = TRUE)
real = kiba_triplet[ord,3]
ord.lower = Preds.lower[ord]
ord.upper = Preds.upper[ord]

triplet_ord = kiba_triplet[ord,]

ind = which(triplet_ord[,2]==20)
real = real[ind]
ord.lower = ord.lower[ind]
ord.upper = ord.upper[ind]
ord.mid = (ord.lower + ord.upper)/2
triplet_ord = triplet_ord[ind,]
# ord = order(real)
ord = order(ord.mid)
real = real[ord]
ord.lower = ord.lower[ord]
ord.upper = ord.upper[ord]
ord.mid = ord.mid[ord]
triplet_ord = triplet_ord[ord,]
plot(real, type='l', ylim = range(c(real, ord.lower, ord.upper)),
     main = "Target 20 from Kiba", xlab = "Index", ylab = "Score")
lines(ord.lower, col=2, lty = 'dashed')
lines(ord.upper, col=3, lty = 'dotdash')

lines(ord.mid, col=4, lty='dotted')

par(mfrow = c(1,1))
