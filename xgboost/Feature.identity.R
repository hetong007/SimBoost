# Helper functions
get.num.obs = function(dict, vec) {
  match.ind = match(vec, dict, 0)
  return (tabulate(match.ind, length(dict)))
}

get.ave.sim = function(dict, sim) {
  return(apply(sim[dict,], 1, mean))
}

get.ave.sim.hist = function(dict, sim, threshold) {
  res = apply(sim[dict,], 1, function(x) {
    threshold = c(threshold, Inf)
    cnt = sapply(threshold, function(thr) sum(x < thr))
    cnt = cnt[-1] - cnt[-length(cnt)]
    return(cnt)
  })
  if (length(threshold) > 1) {
    res = t(res)
  }
  return(res)
}

get.ave.rating = function(dict, inds, rating) {
  ind = match(inds, dict, 0)
  ind = which(ind > 0)
  inds = inds[ind]
  res = tapply(rating[ind], inds, mean)
  return(res)
}

# Get Identity Features
get.identity = function(inds, ratings, sim, threshold) {
  dict = sort(unique(inds))
  
  # Num of observations
  n.obs = get.num.obs(dict, inds)
  
  # Average similarity
  ave.sim = get.ave.sim(dict, sim)
  
  # Average similarity histogram
  ave.sim.hist = get.ave.sim.hist(dict, sim, threshold)
  colnames(ave.sim.hist) = paste0('ave.sim.hist', 1:ncol(ave.sim.hist))
  
  # Average rating
  ave.rating = get.ave.rating(dict, inds, ratings)
  
  # Concatenation
  res = cbind(n.obs, ave.sim, ave.sim.hist, ave.rating)
  colnames(res) = c('n.obs', 'ave.sim', paste0('hist.sim',1:ncol(ave.sim.hist)), 'ave.val')
  return (res)
}

# a = get.identity(dataset_triplets[,1], dataset_triplets[,3], drug_adj_mat, 0.6)
