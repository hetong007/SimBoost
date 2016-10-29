require(xgboost)
source('utils.R')

sequential.mean.cv = function(rdafile, outfile, cutoff,
                              params, nrounds, seed=1024) {
  load(rdafile)
  set.seed(seed)
  
  auc_eval = function(preds, dtrain) {
    labels = getinfo(dtrain, "label")
    auc_res = auc(preds, labels, cutoff)
    return(list(metric = "auc", value = auc_res))
  }
  
  feature.p = ncol(cv.data[[1]][[1]])
  
  nfold = length(cv.data)
  nr = nrow(triplet)
  Preds = rep(0, nr)
  models = vector(nfold, mode='list')
  for (i in 1:nfold) {
    cat(i,'\n')
    train.x = cv.data[[i]][[1]] # $train.x
    train.y = cv.data[[i]][[3]] # $train.y
    test.x = cv.data[[i]][[2]] # $test.x
    test.y = cv.data[[i]][[4]]
    dtrain = xgb.DMatrix(data = train.x, label = train.y)
    dtest = xgb.DMatrix(data = test.x, label = test.y)
    model = xgb.train(params = params, data = dtrain,
                      nrounds = nrounds,
                      watchlist = list(train = dtrain, valid = dtest),
                      feval = auc_eval, maximize = TRUE)
    preds = predict(model, test.x)
    Preds[cv.folds[[i]][[2]]] = preds
    models[[i]] = model
  }
  
  save(Preds, triplet, nrounds, params, file=outfile)
  
  return(cbind(Preds, triplet[,3]))
}
