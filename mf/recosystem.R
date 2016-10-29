require(recosystem)

model.loader = function(input) {
  tmp = readLines(input, n = 3)
  sizes = sapply(strsplit(tmp," "), function(x) as.numeric(x[2]))
  m = sizes[2]
  n = sizes[3]
  k = sizes[4]
  
  P = read.table(input, skip = 5, nrows = m, sep = " ", header = FALSE)
  Q = read.table(input, skip = 5+m, nrows = n, sep = " ", header = FALSE)
  P = P[,-c(1,ncol(P))]
  Q = Q[,-c(1,ncol(Q))]
  P = data.matrix(P)
  Q = data.matrix(Q)
  return(list(P, Q))
}

libmf = function(triplet, m, n, k = 10, cost = 0.1, lrate = 0.1, niter = 50, 
                 nthread = 1, nmf = FALSE, verbose = TRUE) {
  trainfile.temp = tempfile()
  # start from 0
#   triplet[,1] = triplet[,1] - min(triplet[,1])
#   triplet[,2] = triplet[,2] - min(triplet[,2])
  if (max(triplet[,1])>=m || max(triplet[,2])>=n) {
    stop("the row and col numbers start from 0")
  }
  
  write.table(triplet, file = trainfile.temp, sep = " ", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  r = Reco()
  r$train(data_file(trainfile.temp),
          opts = list(dim = k, lrate = lrate, niter = niter,
                      nthread = nthread, nmf = nmf, verbose = verbose,
                      costp_l1 = cost, costp_l2 = cost,
                      costq_l1 = cost, costq_l2 = cost))
  modelfile = r$model$path
  model = model.loader(modelfile)
  
  return(list(r,model))
}

libmf.tune = function(triplet, m, n, k = 10, nfold = 10, cost = 0.1, lrate = 0.1, niter = 20, 
                      nthread = 1, nmf = FALSE, verbose = TRUE) {
  trainfile.temp = tempfile()
  if (max(triplet[,1])>=m || max(triplet[,2])>=n) {
    stop("the row and col numbers start from 0")
  }
  write.table(triplet, file = trainfile.temp, sep = " ", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  r = Reco()
  res = r$tune(data_file(trainfile.temp), 
               opts = list(dim = k, lrate = lrate, niter = niter,
                           nthread = nthread, nmf = nmf, verbose = verbose,
                           nfold = nfold,
                           costp_l1 = cost, costp_l2 = cost,
                           costq_l1 = cost, costq_l2 = cost))
  return(res)
}

predict.libmf = function(model, triplet) {
  r = model[[1]]
  
  testfile.temp = tempfile()
  write.table(triplet, file = testfile.temp, sep = " ", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  output.temp = tempfile()
  r$predict(data_file(testfile.temp), out_file(output.temp))
  pred = readLines(output.temp)
  pred = as.numeric(pred)
  return(pred)
}

libmf.cv = function(triplet, m, n, k = 10, nfold = 10, cost = 0.1, lrate = 0.1, niter = 20, 
                    nthread = 1, nmf = FALSE, verbose = TRUE) {
  N = nrow(triplet)
  folds = vector(nfold, mode='list')
  shuf = sample(N)
  mm = N %/% nfold
  for (i in 1:(nfold-1)) {
    folds[[i]] = shuf[1:mm]
    shuf = setdiff(shuf, folds[[i]])
  }
  folds[[nfold]] = shuf
  preds = rep(0,N)
  
  for (i in 1:nfold) {
    teind = folds[[i]]
    trind = setdiff(1:N, teind)
    res = libmf(triplet[trind, ], m = m, n = n, k = k, cost = cost, lrate = lrate,
                niter = niter, nthread = nthread, nmf = nmf, verbose = verbose)
    preds[teind] = predict.libmf(res, triplet[teind,])
  }
  
  return(preds)
}
