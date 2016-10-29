require(xgboost)
source('utils.R')
load('../data/sequential.cv.feature.mf.kiba.rda')

set.seed(1024)

auc_eval = function(preds, dtrain) {
  labels = getinfo(dtrain, "label")
  auc_res = auc(preds, labels, 12.1)
  return(list(metric = "auc", value = auc_res))
}

nrounds = 1500
feature.p = ncol(cv.data[[1]][[1]])
params = list(nthread = 8, objective = "reg:linear", eta = 0.1, 
              subsample = 0.6, colsample_bytree = 0.6, #sqrt(feature.p)/feature.p,
              num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

nfold = length(cv.data)
nr = nrow(triplet)
Preds = rep(0, nr)
models = vector(nfold, mode='list')
for (i in 1:nfold) {
  train.x = cv.data[[i]][[1]] # $train.x
  train.y = cv.data[[i]][[3]] # $train.y
  test.x = cv.data[[i]][[2]] # $test.x
  test.y = cv.data[[i]][[4]]
  # model = xgboost(data = train.x, label = train.y, params = params, 
  #                 nrounds = nrounds, 
  #                 watchlist = list(train = train.x, valid = test.x))
  dtrain = xgb.DMatrix(data = train.x, label = train.y)
  dtest = xgb.DMatrix(data = test.x, label = test.y)
  model = xgb.train(params = params, data = dtrain,
                    nrounds = nrounds,
                    watchlist = list(train = dtrain, valid = dtest),
                    feval = auc_eval)
  preds = predict(model, test.x)
  Preds[cv.folds[[i]][[2]]] = preds # $test.fold] = preds
  models[[i]] = model
}

rmse(Preds, triplet[,3])
# 0.2007384 
auc(Preds, triplet[,3], 12.1)
# 0.9068976

save(Preds, triplet, nrounds, params, file='../data/sequential.cv.xgb.mf.kiba.rda')
