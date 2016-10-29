require(xgboost)
source('utils.R')

sequential.quantile.cv = function(rdafile, outfile, cutoff,
                                  upper.quantile, lower.quantile,
                                  params.upper, params.lower, nrounds, seed = 1024) {

  load(rdafile)
  set.seed(seed)
  
  auc_eval = function(preds, dtrain) {
    labels = getinfo(dtrain, "label")
    auc_res = auc(preds, labels, cutoff)
    return(list(metric = "auc", value = auc_res))
  }
  
  quant_eval = function(preds, dtrain) {
    labels = getinfo(dtrain, "label")
    n = length(labels)
    sm = sum(preds>labels)
    return(list(metric = "coverage", value = sm/n))
  }
  
  quantile.obj.lower = function(preds, dtrain, alpha = lower.quantile) {
    labels = getinfo(dtrain, "label")
    mask = as.numeric(labels > preds)
    grad = (alpha * mask) - ((1.0 - alpha) * (1-mask))
    grad = grad * labels
    hess = rep(1, length(grad))
    return(list(grad = -grad, hess = hess))
  }
  
  quantile.obj.upper = function(preds, dtrain, alpha = upper.quantile) {
    labels = getinfo(dtrain, "label")
    mask = as.numeric(labels > preds)
    grad = (alpha * mask) - ((1.0 - alpha) * (1-mask))
    grad = grad * labels
    hess = rep(1, length(grad))
    return(list(grad = -grad, hess = hess))
  }
  
  feature.p = ncol(cv.data[[1]][[1]])
  
  nfold = length(cv.data)
  nr = nrow(triplet)
  Preds.lower = rep(0, nr)
  Preds.upper = rep(0, nr)
  models.lower = vector(nfold, mode='list')
  models.upper = vector(nfold, mode='list')
  for (i in 1:nfold) {
    train.x = cv.data[[i]][[1]] # $train.x
    train.y = cv.data[[i]][[3]] # $train.y
    test.x = cv.data[[i]][[2]] # $test.x
    test.y = cv.data[[i]][[4]]

    dtrain = xgb.DMatrix(data = train.x, label = train.y)
    dtest = xgb.DMatrix(data = test.x, label = test.y)
    
    bs.lower = quantile(train.y, 0.5)
    bs.upper = quantile(train.y, 0.5)
    
    params.lower$objective = quantile.obj.lower
    params.lower$base_score = bs.lower
    params.upper$objective = quantile.obj.upper
    params.upper$base_score = bs.upper
    
    cat(i,'\n\n\n')
    model = xgb.train(params = params.lower, data = dtrain,
                      nrounds = nrounds,
                      watchlist = list(train = dtrain, valid = dtest),
                      feval = quant_eval, maximize = TRUE)
    preds = predict(model, test.x)
    Preds.lower[cv.folds[[i]][[2]]] = preds
    models.lower[[i]] = model
    
    model = xgb.train(params = params.upper, data = dtrain,
                      nrounds = nrounds,
                      watchlist = list(train = dtrain, valid = dtest),
                      feval = quant_eval, maximize = TRUE)
    preds = predict(model, test.x)
    Preds.upper[cv.folds[[i]][[2]]] = preds
    models.upper[[i]] = model
  }
  
  Preds = (Preds.lower + Preds.upper)/2
  # res = c(rmse(Preds, triplet[,3]), auc(Preds, triplet[,3], cutoff))
  res = cbind(Preds.lower, Preds.upper, triplet[,3])
  
  save(Preds.lower, Preds.upper, triplet, nrounds, params.lower, params.upper, 
       models.upper, models.lower,
       file=outfile)
  
  return(res)
}

sequential.quantile.plot = function(rdafile, plot.main) {
  
  load(rdafile)
  
  auc.list = c()
  for (x in seq(5,9, length.out = 200)) {
    res = auc(Preds, triplet[,3], x)
    auc.list = c(auc.list, res)
    cat(res,'\n')
  }
  
  plot(seq(5,9, length.out = 200), auc.list, type='l')
  
  
  
  ord = order(Preds.upper-Preds.lower, decreasing = TRUE)
  real = triplet[ord,3]
  ord.lower = Preds.lower[ord]
  ord.upper = Preds.upper[ord]
  
  triplet_ord = triplet[ord,]
  
  ind = which(triplet_ord[,2]==10)
  real = real[ind]
  ord.lower = ord.lower[ind]
  ord.upper = ord.upper[ind]
  triplet_ord = triplet_ord[ind,]
  ord = order(real)
  real = real[ord]
  ord.lower = ord.lower[ord]
  ord.upper = ord.upper[ord]
  triplet_ord = triplet_ord[ord,]
  plot(real, type='l', ylim = range(c(real, ord.lower, ord.upper)),
       main = plot.main, xlab = "Index", ylab = "Score")
  lines(ord.lower, col=2)
  lines(ord.upper, col=3)
}
