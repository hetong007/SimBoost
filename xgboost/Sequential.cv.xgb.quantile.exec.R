require(xgboost)
require(methods)
source('utils.R')
source('Sequential.cv.xgb.quantile.R')

# Loop time
L = 10

parse.mid = function(res) {
  as.numeric(strsplit(gsub('\\s+','',res),',')[[1]][2])
}

rmse.fun = function(preds, label) {
  mean((preds - label)^2)
}

# Davis

davis = list()
rmse = rep(0,L)
params = list(nthread = 8, eta = 0.02,
              subsample = 0.6, colsample_bytree = 0.6,
              num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

for (i in 1:L) {
  davis[[i]] = sequential.quantile.cv('../data/sequential.cv.feature.mf.davis.rda',
                                      '../data/sequential.cv.xgb.quantile.2.mf.davis.rda',
                                      cutoff = 7, 0.95, 0.05, nrounds = 200,
                                      params.upper = params, params.lower = params, 
                                      seed = i)
  res = davis[[i]]
  res = cbind((res[,1] + res[,2])/2, res[,3])
  write.table(res, file=paste0('davis_',i,'_quantile_cont.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
  rmse[i] = rmse.fun(res[,1], res[,2])
  res[,2] = as.numeric(res[,2] > 7)
  write.table(res, file=paste0('davis_',i,'_quantile_bin.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
}
save(davis, file = '../data/xgb.davis.quantile.rda')

auc = rep(0,L)
aupr = rep(0,L)
rmse = rep(0,L)
ci = rep(0,L)
for (i in 1:L) {
  cat(i,'\n')
  res = read.table(paste0('davis_',i,'_quantile_cont.txt'))
  rmse[i] = rmse.fun(res[,1],res[,2])
  res = system(paste0('python ../evaluation/evaluate_AUC_complete.py davis_',i,'_quantile_bin.txt'), intern = TRUE)
  auc[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_AUPR_complete.py davis_',i,'_quantile_bin.txt'), intern = TRUE)
  aupr[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_CI_complete.py davis_',i,'_quantile_cont.txt'), intern = TRUE)
  ci[i] = parse.mid(res)
}

davis.quantile.mean = c(mean(rmse), mean(auc), mean(aupr), mean(ci))
davis.quantile.sd = c(sd(rmse), sd(auc), sd(aupr), sd(ci))
save(davis.quantile.mean, davis.quantile.sd, file='davis.quantile.table.rda')

# Metz

metz = list()
rmse = rep(0,L)
params = list(nthread = 8, eta = 0.01,
              subsample = 0.6, colsample_bytree = 0.6,
              num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

for (i in 1:L) {
  metz[[i]] = sequential.quantile.cv('../data/sequential.cv.feature.mf.metz.rda',
                                     '../data/sequential.cv.xgb.quantile.2.mf.metz.rda',
                                     cutoff = 7.6, 0.95, 0.05, nrounds = 500,
                                     params.upper = params, params.lower = params, 
                                     seed = i)
  res = metz[[i]]
  res = cbind((res[,1] + res[,2])/2, res[,3])
  write.table(res, file=paste0('metz_',i,'_quantile_cont.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
  rmse[i] = rmse.fun(res[,1], res[,2])
  res[,2] = as.numeric(res[,2] > 7.6)
  write.table(res, file=paste0('metz_',i,'_quantile_bin.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
}
save(metz, file = '../data/xgb.metz.quantile.rda')

auc = rep(0,L)
aupr = rep(0,L)
rmse = rep(0,L)
ci = rep(0,L)
for (i in 1:L) {
  cat(i,'\n')
  res = read.table(paste0('metz_',i,'_quantile_cont.txt'))
  rmse[i] = rmse.fun(res[,1],res[,2])
  res = system(paste0('python ../evaluation/evaluate_AUC_complete.py metz_',i,'_quantile_bin.txt'), intern = TRUE)
  auc[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_AUPR_complete.py metz_',i,'_quantile_bin.txt'), intern = TRUE)
  aupr[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_CI_complete.py metz_',i,'_quantile_cont.txt'), intern = TRUE)
  ci[i] = parse.mid(res)
}

metz.quantile.mean = c(mean(rmse), mean(auc), mean(aupr), mean(ci))
metz.quantile.sd = c(sd(rmse), sd(auc), sd(aupr), sd(ci))
save(metz.quantile.mean, metz.quantile.sd, file='metz.quantile.table.rda')


# Kiba
kiba = list()
rmse = rep(0,L)
params = list(nthread = 8, eta = 0.1,
              subsample = 0.6, colsample_bytree = 0.6,
              num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

for (i in 1:L) {
  kiba[[i]] = sequential.quantile.cv('../data/sequential.cv.feature.mf.kiba.rda',
                                     '../data/sequential.cv.xgb.quantile.2.mf.kiba.rda',
                                     cutoff = 12.1, 0.95, 0.05, nrounds = 400,
                                     params.upper = params, params.lower = params, 
                                     seed = i)
  res = kiba[[i]]
  res = cbind((res[,1] + res[,2])/2, res[,3])
  write.table(res, file=paste0('kiba_',i,'_quantile_cont.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
  rmse[i] = rmse.fun(res[,1], res[,2])
  res[,2] = as.numeric(res[,2] > 12.1)
  write.table(res, file=paste0('kiba_',i,'_quantile_bin.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
  save(i, kiba, file = '../data/xgb.kiba.quantile.rda')
}

auc = rep(0,L)
aupr = rep(0,L)
rmse = rep(0,L)
ci = rep(0,L)
for (i in 1:L) {
  cat(i,'\n')
  res = read.table(paste0('kiba_',i,'_quantile_cont.txt'))
  rmse[i] = rmse.fun(res[,1],res[,2])
  res = system(paste0('python ../evaluation/evaluate_AUC_complete.py kiba_',i,'_quantile_bin.txt'), intern = TRUE)
  auc[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_AUPR_complete.py kiba_',i,'_quantile_bin.txt'), intern = TRUE)
  aupr[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_CI_complete.py kiba_',i,'_quantile_cont.txt'), intern = TRUE)
  ci[i] = parse.mid(res)
}

kiba.quantile.mean = c(mean(rmse), mean(auc), mean(aupr), mean(ci))
kiba.quantile.sd = c(sd(rmse), sd(auc), sd(aupr), sd(ci))
save(kiba.quantile.mean, kiba.quantile.sd, file='kiba.quantile.table.rda')
