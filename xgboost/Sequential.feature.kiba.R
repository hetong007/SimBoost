source('Sequential.cv.R')

load('../data/kiba_data.Rda')
triplet = kiba_triplet
triplet = triplet[order(triplet[,1]),]
d.sim = kiba_drug_sim
t.sim = kiba_target_sim
# load('../data/5fold_cv_data_kiba.rda')
# triplet = dataset_triplets
# triplet = triplet[order(triplet[,1]),]
# d.sim = drug_adj_mat
# t.sim = target_adj_mat

res = clear.cold.start(triplet, d.sim, t.sim, 5)
triplet = res[[1]]
d.sim = res[[2]]
t.sim = res[[3]]
nfold = 5
k = 4
latent.dim = 10
threshold = 0.3
threshold.identity = (0:5)/10

set.seed(1024)
seq.feature = sequential.cv.feature(triplet, d.sim, t.sim, latent.dim,
                                    nfold, k, threshold, threshold.identity)
cv.data = seq.feature[[1]]
cv.folds = seq.feature[[2]]

save(triplet, d.sim, t.sim, cv.data, cv.folds, file='../data/sequential.cv.feature.mf.kiba.rda')




