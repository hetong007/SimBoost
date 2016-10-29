# Weighted ratings from d to t's neighbours (if any)
# Weighted ratings from t to d's neighbours (if any)
# Paths from d to t if d-t is not observed
# Shortest distance to all the other same type nodes
# PageRank

require(igraph)
source('Feature.identity.R')
source('Feature.homo.R')
source('recosystem.R')

# Helper Functions

get.w.rating = function(ind1, ind2, ratings, knn2) {
  dict = get.dict(ind1)
  res = rep(0, length(ratings))
  k = ncol(knn2)
  l = 1
  for (dd in dict) {
    ind = which(ind1 == dd)
    ind2.sub = ind2[ind]
    knn2.sub = knn2[ind2.sub,,drop = FALSE]
    rating.sub = ratings[ind]
    m = nrow(knn2.sub)
    tmp = rep(0, m)
    for (i in 1:m) {
      n.ind = intersect(knn2.sub[i,], ind2.sub)
      if (length(n.ind) > 0) {
        t.ind = match(n.ind, ind2.sub, 0)
        tmp[i] = mean(rating.sub[t.ind])
      } else {
        tmp[i] = 0
      }
    }
    res[l:(l+m-1)] = tmp
    l = l+m
  }
  return(res)
}

get.dict = function(inds) {
  return(sort(unique(inds)))
}

get.hetero.net = function(d.inds, t.inds, ratings) {
  d.dict = get.dict(d.inds)
  t.dict = get.dict(t.inds)
  
  nr = length(d.dict)
  nc = length(t.dict)
  net = matrix(0, nr, nc)
  
  m = length(d.inds)
  for (i in 1:m) {
    net[d.inds[i], t.inds[i]] = 1 #ratings[i]
  }
  return(net)
}

get.mf = function(triplet, latent.dim, ...) {
  m = max(triplet[,1])
  n = max(triplet[,2])
  triplet[,1] = triplet[,1] - 1
  triplet[,2] = triplet[,2] - 1
  model = libmf(triplet = triplet, m = m, n = n, k = latent.dim, verbose = FALSE, ...)
  return(model[[2]])
}

get.hetero.graph = function(d.inds, t.inds) {
  m = length(d.inds)
  D = length(unique(d.inds))
  edges = matrix(0, 2, m)
  wt = rep(1, m)
  for (i in 1:m) {
    edges[,i] = c(d.inds[i], t.inds[i] + D)
  }
  edges = as.vector(edges)
  graph = make_graph(edges, directed = FALSE)
  graph = set.edge.attribute(graph, "weight", value=wt)
  return(graph)
}

get.heterogenous.feature = function(d.inds, t.inds, ratings, k,
                                    d.test.inds, t.test.inds,
                                    d.sim, t.sim,
                                    threshold, threshold.identity,
                                    latent.dim, ...) {
  d.dict = get.dict(d.inds)
  t.dict = get.dict(t.inds)
  
  d.net = get.homo.graph(d.sim, threshold)
  d.knn = get.knn(d.net, k)
  # d.graph = get.graph.from.net(d.net)
  
  t.net = get.homo.graph(t.sim, threshold)
  t.knn = get.knn(t.net, k)
  # t.graph = get.graph.from.net(t.net)
  
  m = length(d.test.inds)
  d.t.w.rating = get.w.rating(c(d.inds, d.test.inds),
                              c(t.inds, t.test.inds),
                              c(ratings, rep(0,m)), t.knn)
  t.d.w.rating = get.w.rating(c(t.inds, t.test.inds),
                              c(d.inds, d.test.inds),
                              c(ratings, rep(0, m)), d.knn)
  
  # hetero.net = get.hetero.net(d.inds, t.inds, ratings)
  # hetero.graph = get.hetero.graph(d.inds, t.inds)
  # cl = igraph::closeness(hetero.graph)
  # d.cl = cl[1:length(d.dict)]
  d.cl = rep(1, length(d.dict))
  # t.cl = cl[-(1:length(d.dict))]
  t.cl = rep(1, length(t.dict))
  # pr = igraph::page_rank(hetero.graph, directed = FALSE)
  # d.pr = pr$vector[1:length(d.dict)]
  d.pr = rep(1, length(d.dict))
  # t.pr = pr$vector[-(1:length(d.dict))]
  t.pr = rep(1, length(t.dict))
  
  triplet = cbind(d.inds, t.inds, ratings)
  mf.latentvector = get.mf(triplet, latent.dim, ...)
  d.mf = mf.latentvector[[1]]
  colnames(d.mf) = paste0('d.mf', 1:ncol(d.mf))
  t.mf = mf.latentvector[[2]]
  colnames(t.mf) = paste0('t.mf', 1:ncol(t.mf))
  
  weighted.rating = cbind(d.t.w.rating, t.d.w.rating)
  colnames(weighted.rating) = c('d.t.ave', 't.d.ave')
  d.feature = cbind(d.cl, d.pr, d.mf)
  t.feature = cbind(t.cl, t.pr, t.mf)
  res = list(weighted.rating, d.feature, t.feature)             
  return(res)
}

# debug(get.heterogenous.feature)
# a = get.heterogenous.feature(dataset_triplets[,1], dataset_triplets[,2], 
#                              dataset_triplets[,3], 4,
#                              drug_adj_mat, target_adj_mat,
#                              0.6, 0.6)
