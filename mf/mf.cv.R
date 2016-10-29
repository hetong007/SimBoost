source('recosystem.R')

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

mf.cv = function(triplet, latent.dim, nfold, seed) {
  set.seed(seed)
  m = max(triplet[,1])
  n = max(triplet[,2])
  cv.folds = get.cv.fold(triplet, nfold)
  triplet[,1] = triplet[,1] - 1
  triplet[,2] = triplet[,2] - 1
  Preds = rep(0, nrow(triplet))
  for (i in 1:nfold) {
    cat(i,'\n')
    train.fold = cv.folds[[i]][[1]] # $train.fold
    test.fold = cv.folds[[i]][[2]] # $test.fold
    train.triplet = triplet[train.fold,]
    test.triplet = triplet[test.fold,]
    model = libmf(train.triplet, m = m, n = n,
                  k = latent.dim, nthread = 8)
    preds = predict.libmf(model, test.triplet)
    Preds[test.fold] = preds
  }
  return(cbind(Preds, triplet[,3]))
}
