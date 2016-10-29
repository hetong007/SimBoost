require(xgboost)
source('Feature.heterous.R')

similarity.knn = function(sim, k) {
  n = nrow(sim)
  matind = matrix(0, n, n)
  for (i in 1:n) {
    ind = order(sim[i,], decreasing = TRUE)[1:k]
    matind[i, ind] = 1
  }
  matind = 1-(1-matind)*(1-t(matind))
  sim = sim*matind
  diag(sim) = 0
  return(sim)
}

clear.cold.start = function(triplet, d.sim, t.sim, min.num) {
  flag = FALSE
  while (!flag) {
    flag = TRUE
    d.dict = sort(unique(triplet[,1]))
    for (dd in d.dict) {
      ind = which(triplet[,1] == dd)
      if (length(ind) < min.num) {
        flag = FALSE
        triplet = triplet[-ind,]
        # d.sim = d.sim[-dd, -dd]
      }
    }
    t.dict = sort(unique(triplet[,2]))
    for (tt in t.dict) {
      ind = which(triplet[,2] == tt)
      if (length(ind) < min.num) {
        flag = FALSE
        triplet = triplet[-ind,]
        # t.sim = t.sim[-tt, -tt]
      }
    }
  }
  d.dict = sort(unique(triplet[,1]))
  t.dict = sort(unique(triplet[,2]))
  d.sim = d.sim[d.dict, d.dict]
  t.sim = t.sim[t.dict, t.dict]
  triplet[,1] = match(triplet[,1], d.dict)
  triplet[,2] = match(triplet[,2], t.dict)
  return(list(triplet, d.sim, t.sim))
}

feature.construction = function(train.triplet, test.triplet, 
                                d.sim, t.sim, k, latent.dim,
                                threshold, threshold.identity) {
  #d.sim = similarity.knn(d.sim, 4)
  #t.sim = similarity.knn(t.sim, 5)
  
  d.identity = get.identity(train.triplet[,1], train.triplet[,3], 
                            d.sim, threshold.identity)
  colnames(d.identity) = paste0('d.', colnames(d.identity))
  t.identity = get.identity(train.triplet[,2], train.triplet[,3], 
                            t.sim, threshold.identity)
  colnames(t.identity) = paste0('t.', colnames(t.identity))
  
  d.homo = get.homogeneous.feature(train.triplet[,1], train.triplet[,3],
                                   d.sim, 4, threshold, threshold.identity)
  colnames(d.homo) = paste0('d.', colnames(d.homo))
  t.homo = get.homogeneous.feature(train.triplet[,2], train.triplet[,3],
                                   t.sim, 4, threshold, threshold.identity)
  colnames(t.homo) = paste0('t.', colnames(t.homo))
  
  f.hetero = get.heterogenous.feature(train.triplet[,1], train.triplet[,2],
                                      train.triplet[,3], 4,
                                      test.triplet[,1], test.triplet[,2],
                                      d.sim, t.sim, 0.6, 0.6, latent.dim)
  d.hetero = f.hetero[[2]]
  t.hetero = f.hetero[[3]]
  d.t.hetero = f.hetero[[1]]
  
  d.feature = cbind(d.identity, d.homo, d.hetero)
  # colnames(d.feature) = paste0('d.', colnames(d.feature))
  t.feature = cbind(t.identity, t.homo, t.hetero)
  # colnames(t.feature) = paste0('t.', colnames(t.feature))
  
  train.feature = cbind(d.feature[train.triplet[,1],],
                        t.feature[train.triplet[,2],],
                        d.t.hetero[1:nrow(train.triplet),])
  test.feature = cbind(d.feature[test.triplet[,1],],
                       t.feature[test.triplet[,2],],
                       d.t.hetero[1:nrow(test.triplet),])
  train.y = train.triplet[,3]
  test.y = test.triplet[,3]
  return(list(train.feature, test.feature, train.y, test.y))
}

get.cv.fold = function(triplet, nfold) {
  n = nrow(triplet)
  train.fold = vector(nfold, mode='list')
  test.fold = vector(nfold, mode='list')
  
  # label = sample(1:nfold, n, replace = TRUE)
  label = rep(sample(1:nfold), n%/%nfold)
  m = n-length(label)
  label = c(label, sample(1:nfold, m))
  t.dict = sort(unique(triplet[,2]))
  for (tt in t.dict) {
    ind = which(triplet[,2] == tt)
    if (length(unique(label[ind])) == 1) {
      cat(tt,'\n')
      cat(label[ind],'\n')
      flag = FALSE
      l = 1
      while (!flag && l <= length(ind)) {
        ind.x = ind[l]
        dd = triplet[ind.x, 1]
        ind.d = which(triplet[,1] == dd)
        new.lb = nfold + 1 - label[ind.x]
        new.d.label = label[ind.d]
        new.d.label[l] = new.lb
        if (length(unique(new.d.label)) > 1) {
          flag = TRUE
          label[ind.x] = new.lb
        } else {
          l = l+1
        }
      }
      if (!flag) {
        stop("Not able to do cv folding!")
      }
    }
  }
  for (i in 1:nfold) {
    test.fold[[i]] = which(label == i)
  }
  # for (i in 1:n) {
  #   test.fold[[(i %% nfold) + 1]]
  # }
  for (i in 1:nfold) {
    train.fold[[i]] = setdiff(1:n, test.fold[[i]])
  }
  cv.folds = vector(nfold, mode='list')
  for (i in 1:nfold) {
    cv.folds[[i]] = list(train.fold[[i]], test.fold[[i]])
  }
  return(cv.folds)
}

sequential.cv.feature = function(triplet, d.sim, t.sim, latent.dim,
                                 nfold, k, threshold, threshold.identity) {
  n = nrow(triplet)
  cv.folds = get.cv.fold(triplet, nfold)
  cv.models = vector(nfold, mode='list')
  cv.preds = rep(0, n)
  cv.data = vector(nfold, mode='list')
  for (i in 1:nfold) {
    cat(i,'\n')
    train.fold = cv.folds[[i]][[1]] # $train.fold
    test.fold = cv.folds[[i]][[2]] # $test.fold
    train.triplet = triplet[train.fold,]
    test.triplet = triplet[test.fold,]
    features = feature.construction(train.triplet, test.triplet, 
                                    d.sim, t.sim, k, latent.dim,
                                    threshold, threshold.identity)
    train.x = features[[1]]
    train.y = features[[3]]
    
    test.x = features[[2]]
    test.y = features[[4]]
    cv.data[[i]] = list(train.x, test.x, train.y, test.y)
  }
  return(list(cv.data, cv.folds))
}

