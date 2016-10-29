require(xgboost)
source('utils.R')
source('Sequential.cv.R')
source('Sequential.cv.xgb.R')

load('../data/kiba_data.Rda')

# Feature extraction

features = feature.construction(kiba_triplet, kiba_triplet[1:5,],
                                kiba_drug_sim, kiba_target_sim,
                                k=4, latent.dim=10, 0.3, (0:5)/10)
train.features = features[[1]]
train.y = features[[3]]

# xgboost training

params = list(nthread = 8, objective = "reg:linear", eta = 0.2, 
              subsample = 0.8, colsample_bytree = 0.8,
              num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

auc_eval = function(preds, dtrain) {
  labels = getinfo(dtrain, "label")
  auc_res = auc(preds, labels, 12.1)
  return(list(metric = "auc", value = auc_res))
}

dtrain = xgb.DMatrix(data = train.features, label = train.y)
model = xgb.train(params = params, data = dtrain,
                  nrounds = 1500,
                  watchlist = list(train = dtrain),
                  feval = auc_eval, maximize = TRUE)

# Model inspection

nms = colnames(train.features)

importance_matrix = xgb.importance(colnames(train.features), model=model)
xgb.ggplot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Relative importance")
xgb.plot.importance(importance_matrix[1:15,], rel_to_first = TRUE, 
                    xlab = "Relative importance", cex = 1.3)

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

params.lower = list(nthread = 8, eta = 0.1, base_score = quantile(train.y, 0.5),
                    subsample = 0.6, colsample_bytree = 0.6, objective = quantile.obj.lower,
                    num_parallel_tree = 1, max_depth = 6, min_child_weight = 10)

params.upper = list(nthread = 8, eta = 0.1, base_score = quantile(train.y, 0.5),
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

real = kiba_triplet[,3]
save(Preds.lower, Preds.upper, real, file='kiba.inspection.rda')

# par(mfrow = c(1,3))
# par(mfrow = c(1,2))

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
plot(real, type='l', 
     main = "Target 20 in KIBA", xlab = "Index", ylab = "Score", ylim = c(8,17),
     cex.axis = 1.5, cex.lab = 1.5)
lines(ord.lower, col=2)
lines(ord.upper, col=3)

lines(ord.mid, col=4)

ord = order(Preds.upper-Preds.lower, decreasing = TRUE)
real = kiba_triplet[ord,3]
ord.lower = Preds.lower[ord]
ord.upper = Preds.upper[ord]

triplet_ord = kiba_triplet[ord,]

ind = which(triplet_ord[,2]==40)
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
plot(real, type='l',
     main = "Target 40 in KIBA", xlab = "Index", ylab = "Score", ylim = c(8,17),
     cex.axis = 1.5, cex.lab = 1.5)
lines(ord.lower, col=2)
lines(ord.upper, col=3)

lines(ord.mid, col=4)

ord = order(Preds.upper-Preds.lower, decreasing = TRUE)
real = kiba_triplet[ord,3]
ord.lower = Preds.lower[ord]
ord.upper = Preds.upper[ord]

triplet_ord = kiba_triplet[ord,]

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
plot(real, type='l', 
     main = "Target 30 in KIBA", xlab = "Index", ylab = "Score", ylim = c(8,17),
     cex.axis = 1.5, cex.lab = 1.5)
lines(ord.lower, col=2)
lines(ord.upper, col=3)

lines(ord.mid, col=4)

par(mfrow = c(1,1))


# Analysis Width-Observation

triplet = kiba_triplet
triplet = cbind(triplet, Preds.lower, (Preds.lower+Preds.upper)/2, 
                Preds.upper, Preds.upper-Preds.lower)
ave.width = tapply(triplet[,7], triplet[,2], mean)
ave.obs = tapply(triplet[,7], triplet[,2], length)
plot(ave.width~ ave.obs, pch = 20, xlab = "Number of Observations", 
     ylab = "Average Prediction Interval Width", main = "Number of observations v.s. Interval Width",
     cex.axis = 1.5, cex.lab = 1.5)
model = lm(ave.width~ave.obs)
abline(model$coefficients, col = 2)
