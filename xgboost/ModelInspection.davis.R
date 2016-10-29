require(xgboost)
source('utils.R')
source('Sequential.cv.R')
source('Sequential.cv.xgb.R')

load('../data/davis_data.Rda')

# Feature extraction

features = feature.construction(davis_triplet, davis_triplet[1:5,],
                                davis_drug_sim, davis_target_sim,
                                k=4, latent.dim=10, 0.3, (0:5)/10)
train.features = features[[1]]
train.y = features[[3]]

# xgboost training

params = list(nthread = 8, objective = "reg:linear", eta = 0.05, 
              subsample = 0.6, colsample_bytree = 0.6, #sqrt(feature.p)/feature.p,
              num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

auc_eval = function(preds, dtrain) {
  labels = getinfo(dtrain, "label")
  auc_res = auc(preds, labels, 7)
  return(list(metric = "auc", value = auc_res))
}

dtrain = xgb.DMatrix(data = train.features, label = train.y)
model = xgb.train(params = params, data = dtrain,
                  nrounds = 1000,
                  watchlist = list(train = dtrain),
                  feval = auc_eval, maximize = TRUE)

# Model inspection

nms = colnames(train.features)

importance_matrix = xgb.importance(colnames(train.features), model=model)
xgb.ggplot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Relative importance")
xgb.plot.importance(importance_matrix[1:15,], rel_to_first = TRUE, 
                    xlab = "Relative importance of features", cex = 1.5, cex.lab = 2)

# merge MF
ind = grep('mf',importance_matrix[,Feature])
importance_matrix_new = importance_matrix[-ind,]
newval = matrix(colSums(as.data.frame(importance_matrix[ind,])[,-1]),nrow = 1)
to.append = data.frame(Feature = 'MF', newval)
names(to.append) = names(importance_matrix)
importance_matrix_new = rbind(to.append, importance_matrix_new)

xgb.plot.importance(importance_matrix_new[1:15,], rel_to_first = TRUE, 
                    xlab = "Relative importance of features", cex = 1.5, cex.lab = 2)

# Learn Quantile

dtrain = xgb.DMatrix(data = train.features, label = train.y)

quant_eval = function(preds, dtrain) {
  labels = getinfo(dtrain, "label")
  n = length(labels)
  sm = sum(preds>labels)
  return(list(metric = "coverage", value = sm/n))
}

quantile.obj.lower = function(preds, dtrain) {
  labels = getinfo(dtrain, "label")
  mask = as.numeric(labels > preds)
  grad = (0.05 * mask) - ((1.0 - 0.05) * (1-mask))
  grad = grad * labels
  hess = rep(1, length(grad))
  return(list(grad = -grad, hess = hess))
}

quantile.obj.upper = function(preds, dtrain) {
  labels = getinfo(dtrain, "label")
  mask = as.numeric(labels > preds)
  grad = (0.95 * mask) - ((1.0 - 0.95) * (1-mask))
  grad = grad * labels
  hess = rep(1, length(grad))
  return(list(grad = -grad, hess = hess))
}

params.lower = list(nthread = 8, eta = 0.02, base_score = quantile(train.y, 0.5),
                    subsample = 0.6, colsample_bytree = 0.6, objective = quantile.obj.lower,
                    num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

params.upper = list(nthread = 8, eta = 0.02, base_score = quantile(train.y, 0.5),
                    subsample = 0.6, colsample_bytree = 0.6, objective = quantile.obj.upper,
                    num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

model.lower = xgb.train(params = params.lower, data = dtrain,
                        nrounds = 200,
                        watchlist = list(train = dtrain),
                        feval = quant_eval, maximize = TRUE)
Preds.lower = predict(model.lower, dtrain)

model.upper = xgb.train(params = params.upper, data = dtrain,
                        nrounds = 200,
                        watchlist = list(train = dtrain),
                        feval = quant_eval, maximize = TRUE)
Preds.upper = predict(model.upper, dtrain)

real = davis_triplet[,3]
save(Preds.lower, Preds.upper, real, file='davis.inspection.rda')

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
     main = "Davis Prediction", xlab = "Index", ylab = "Score")
lines(ord.lower, col=2)
lines(ord.upper, col=3)

lines(ord.mid, col=4)
