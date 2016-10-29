require(ROCR)

rmse = function(label, preds) {
  return(mean((label-preds)^2))
}

auc = function(preds, label, cutoff) {
  label.cut = as.numeric(label > cutoff)
  # if (sum(label.cut) == 0 || sum(label.cut) == length(label))
  pred.obj = prediction(preds, label.cut)
  roc.perf = performance(pred.obj, "auc")
  return(roc.perf@y.values[[1]])
}

