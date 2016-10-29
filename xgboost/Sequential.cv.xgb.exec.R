require(xgboost)
require(methods)
source('utils.R')
source('Sequential.cv.xgb.R')

#Looptime
L = 10

parse.mid = function(res) {
  as.numeric(strsplit(gsub('\\s+','',res),',')[[1]][2])
}

rmse.fun = function(preds, label) {
  mean((preds - label)^2)
}

# Davis

# Load original data
source('Sequential.cv.R')
load('../data/davis_data.Rda')
triplet = davis_triplet
triplet = triplet[order(triplet[,1]),]
d.sim = davis_drug_sim
t.sim = davis_target_sim
res = clear.cold.start(triplet, d.sim, t.sim, 5)
triplet = res[[1]]
d.sim = res[[2]]
t.sim = res[[3]]
nfold = 5
k = 4
latent.dim = 10
threshold = 0.3
threshold.identity = (0:5)/10

davis = list()
rmse = rep(0,L)
params = list(nthread = 8, objective = "reg:linear", eta = 0.05, 
              subsample = 0.6, colsample_bytree = 0.6, #sqrt(feature.p)/feature.p,
              num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

for (i in 1:L) {
  cat('Round ',i,'\n')
  set.seed(i)
  seq.feature = sequential.cv.feature(triplet, d.sim, t.sim, latent.dim,
                                      nfold, k, threshold, threshold.identity)
  cv.data = seq.feature[[1]]
  cv.folds = seq.feature[[2]]
  
  save(triplet, d.sim, t.sim, cv.data, cv.folds, file='../data/sequential.cv.feature.mf.davis.rda')
  
  
  davis[[i]] = sequential.mean.cv('../data/sequential.cv.feature.mf.davis.rda',
                                  '../data/sequential.cv.xgb.mf.davis.rda',
                                  cutoff = 7, params = params, nrounds = 1000, seed = i)
  res = davis[[i]]
  res = cbind(res[,1],res[,2])
  auc(res[,1],res[,2],7)
  write.table(res, file=paste0('davis_',i,'_cont.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
  rmse[i] = rmse.fun(res[,1], res[,2])
  res[,2] = as.numeric(res[,2] > 7)
  write.table(res, file=paste0('davis_',i,'_bin.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
}

save(davis, file='../data/xgb.davis.rda')

sapply(1:10, function(i) auc(davis[[i]][,1],davis[[i]][,2], 7))

auc = rep(0,L)
aupr = rep(0,L)
ci = rep(0,L)
rmse = rep(0,L)
for (i in 1:L) {
  cat(i,'\n')
  res = read.table(paste0('davis_',i,'_cont.txt'))
  rmse[i] = rmse.fun(res[,1],res[,2])
  res = system(paste0('python ../evaluation/evaluate_AUC_complete.py davis_',i,'_bin.txt'), intern = TRUE)
  auc[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_AUPR_complete.py davis_',i,'_bin.txt'), intern = TRUE)
  aupr[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_CI_complete.py davis_',i,'_cont.txt'), intern = TRUE)
  ci[i] = parse.mid(res)
}

davis.mean = c(mean(rmse), mean(auc), mean(aupr), mean(ci))
davis.sd = c(sd(rmse), sd(auc), sd(aupr), sd(ci))
save(davis.mean, davis.sd, file='davis.table.rda')

# Metz


load('../data/metz_data.Rda')
triplet = metz_triplet
triplet = triplet[order(triplet[,1]),]
d.sim = metz_drug_sim
t.sim = metz_target_sim
res = clear.cold.start(triplet, d.sim, t.sim, 5)
triplet = res[[1]]
d.sim = res[[2]]
t.sim = res[[3]]
nfold = 5
k = 4
latent.dim = 10
threshold = 0.3
threshold.identity = (0:5)/10

metz = list()
rmse = rep(0,L)
params = list(nthread = 8, objective = "reg:linear", eta = 0.1, 
              subsample = 0.6, colsample_bytree = 0.6, #sqrt(feature.p)/feature.p,
              num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

for (i in 1:L) {
  cat('Round ',i,'\n')
  set.seed(i)
  seq.feature = sequential.cv.feature(triplet, d.sim, t.sim, latent.dim,
                                      nfold, k, threshold, threshold.identity)
  cv.data = seq.feature[[1]]
  cv.folds = seq.feature[[2]]
  
  save(triplet, d.sim, t.sim, cv.data, cv.folds, file='../data/sequential.cv.feature.mf.metz.rda')
  
  metz[[i]] = sequential.mean.cv('../data/sequential.cv.feature.mf.metz.rda',
                                  '../data/sequential.cv.xgb.mf.metz.rda',
                                  cutoff = 7.6, params = params, nrounds = 1000, seed = i)
  res = metz[[i]]
  res = cbind(res[,1],res[,2])
  auc(res[,1],res[,2],7.6)
  write.table(res, file=paste0('metz_',i,'_cont.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
  rmse[i] = rmse.fun(res[,1], res[,2])
  res[,2] = as.numeric(res[,2] > 7.6)
  write.table(res, file=paste0('metz_',i,'_bin.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
}

save(metz, file='../data/xgb.metz.rda')

auc = rep(0,L)
aupr = rep(0,L)
rmse = rep(0,L)
ci = rep(0,L)
for (i in 1:L) {
  cat(i,'\n')
  res = read.table(paste0('metz_',i,'_cont.txt'))
  rmse[i] = rmse.fun(res[,1],res[,2])
  res = system(paste0('python ../evaluation/evaluate_AUC_complete.py metz_',i,'_bin.txt'), intern = TRUE)
  auc[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_AUPR_complete.py metz_',i,'_bin.txt'), intern = TRUE)
  aupr[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_CI_complete.py metz_',i,'_cont.txt'), intern = TRUE)
  ci[i] = parse.mid(res)
}

metz.mean = c(mean(rmse), mean(auc), mean(aupr), mean(ci))
metz.sd = c(sd(rmse), sd(auc), sd(aupr), sd(ci))
save(metz.mean, metz.sd, file='metz.table.rda')

# Kiba

load('../data/kiba_data.Rda')
triplet = kiba_triplet
triplet = triplet[order(triplet[,1]),]
d.sim = kiba_drug_sim
t.sim = kiba_target_sim
res = clear.cold.start(triplet, d.sim, t.sim, 5)
triplet = res[[1]]
d.sim = res[[2]]
t.sim = res[[3]]
nfold = 5
k = 4
latent.dim = 10
threshold = 0.3
threshold.identity = (0:5)/10


kiba = list()
rmse = rep(0,L)
params = list(nthread = 8, objective = "reg:linear", eta = 0.2, 
              subsample = 0.8, colsample_bytree = 0.8,
              num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

for (i in 1:L) {
  cat('Round ',i,'\n')
  set.seed(i)
  seq.feature = sequential.cv.feature(triplet, d.sim, t.sim, latent.dim,
                                      nfold, k, threshold, threshold.identity)
  cv.data = seq.feature[[1]]
  cv.folds = seq.feature[[2]]
  
  save(triplet, d.sim, t.sim, cv.data, cv.folds, file='../data/sequential.cv.feature.mf.kiba.rda')
  
  kiba[[i]] = sequential.mean.cv('../data/sequential.cv.feature.mf.kiba.rda',
                                 '../data/sequential.cv.xgb.mf.kiba.rda',
                                 cutoff = 12.1, params = params, nrounds = 1500, seed = i)
  res = kiba[[i]]
  res = cbind(res[,1],res[,2])
  auc(res[,1],res[,2],12.1)
  write.table(res, file=paste0('kiba_',i,'_cont.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
  rmse[i] = rmse.fun(res[,1], res[,2])
  res[,2] = as.numeric(res[,2] > 12.1)
  write.table(res, file=paste0('kiba_',i,'_bin.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
}

save(kiba, file='../data/xgb.kiba.rda')

auc = rep(0,L)
aupr = rep(0,L)
rmse = rep(0,L)
ci = rep(0,L)
for (i in 1:L) {
  cat(i,'\n')
  res = read.table(paste0('kiba_',i,'_cont.txt'))
  rmse[i] = rmse.fun(res[,1],res[,2])
  res = system(paste0('python ../evaluation/evaluate_AUC_complete.py kiba_',i,'_bin.txt'), intern = TRUE)
  auc[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_AUPR_complete.py kiba_',i,'_bin.txt'), intern = TRUE)
  aupr[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_CI_complete.py kiba_',i,'_cont.txt'), intern = TRUE)
  ci[i] = parse.mid(res)
}

kiba.mean = c(mean(rmse), mean(auc), mean(aupr), mean(ci))
kiba.sd = c(sd(rmse), sd(auc), sd(aupr), sd(ci))
save(kiba.mean, kiba.sd, file='kiba.table.rda')


