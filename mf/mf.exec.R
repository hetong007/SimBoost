source('recosystem.R')
source('mf.cv.R')

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
load('../data/davis_data.Rda')
for (i in 1:L) {
  davis[[i]] = mf.cv(davis_triplet, 10, 5, i)
  res = davis[[i]]
  write.table(res, file=paste0('davis_',i,'_cont.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
  rmse[i] = rmse.fun(res[,1], res[,2])
  res[,2] = as.numeric(res[,2] > 7)
  write.table(res, file=paste0('davis_',i,'_bin.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
}
save(davis, file = '../data/mf.davis.rda')

auc = rep(0,L)
aupr = rep(0,L)
ci = rep(0,L)
for (i in 1:L) {
  res = system(paste0('python ../evaluation/evaluate_AUC_complete.py davis_',i,'_bin.txt'), intern = TRUE)
  auc[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_AUPR_complete.py davis_',i,'_bin.txt'), intern = TRUE)
  aupr[i] = parse.mid(res)
  res = system(paste0('python ../evaluation/evaluate_CI_complete.py davis_',i,'_cont.txt'), intern = TRUE)
  ci[i] = parse.mid(res)
}

c(mean(rmse), mean(auc), mean(aupr), mean(ci))
c(sd(rmse), sd(auc), sd(aupr), sd(ci))

# Metz

metz = list()
rmse = rep(0,L)
load('../data/metz_data.Rda')
for (i in 1:L) {
  metz[[i]] = mf.cv(metz_triplet, 10, 5, i)
  res = metz[[i]]
  write.table(res, file=paste0('metz_',i,'_cont.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
  rmse[i] = rmse.fun(res[,1], res[,2])
  res[,2] = as.numeric(res[,2] > 7.6)
  write.table(res, file=paste0('metz_',i,'_bin.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
}
save(metz, file = '../data/mf.metz.rda')

auc = rep(0,L)
aupr = rep(0,L)
ci = rep(0,L)
for (i in 1:L) {
  cat(i,'\n')
  res = system(paste0('python ../evaluation/evaluate_AUC_complete.py metz_',i,'_bin.txt'), intern = TRUE)
  auc[i] = parse.mid(res)
  cat(auc[i],' ')
  res = system(paste0('python ../evaluation/evaluate_AUPR_complete.py metz_',i,'_bin.txt'), intern = TRUE)
  aupr[i] = parse.mid(res)
  cat(aupr[i],' ')
  res = system(paste0('python ../evaluation/evaluate_CI_complete.py metz_',i,'_cont.txt'), intern = TRUE)
  ci[i] = parse.mid(res)
  cat(ci[i],' ')
}

c(mean(rmse), mean(auc), mean(aupr), mean(ci))
c(sd(rmse), sd(auc), sd(aupr), sd(ci))

# Kiba

kiba = list()
rmse = rep(0, L)
load('../data/kiba_data.Rda')
for (i in 1:L) {
  kiba[[i]] = mf.cv(kiba_triplet, 10, 5, i)
  res = kiba[[i]]
  write.table(res, file=paste0('kiba_',i,'_cont.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
  rmse[i] = rmse.fun(res[,1], res[,2])
  res[,2] = as.numeric(res[,2] > 12.1)
  write.table(res, file=paste0('kiba_',i,'_bin.txt'), col.names = FALSE,quote = FALSE,
              row.names = FALSE)
}
save(kiba, file = '../data/mf.kiba.rda')

auc = rep(0,L)
aupr = rep(0,L)
ci = rep(0,L)
for (i in 1:L) {
  cat(i,'\n')
  res = system(paste0('python ../evaluation/evaluate_AUC_complete.py kiba_',i,'_bin.txt'), intern = TRUE)
  auc[i] = parse.mid(res)
  cat(auc[i],' ')
  res = system(paste0('python ../evaluation/evaluate_AUPR_complete.py kiba_',i,'_bin.txt'), intern = TRUE)
  aupr[i] = parse.mid(res)
  cat(aupr[i],' ')
  res = system(paste0('python ../evaluation/evaluate_CI_complete.py kiba_',i,'_cont.txt'), intern = TRUE)
  ci[i] = parse.mid(res)
  cat(ci[i],' ')
}

c(mean(rmse), mean(auc), mean(aupr), mean(ci))
c(sd(rmse), sd(auc), sd(aupr), sd(ci))