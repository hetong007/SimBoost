source('Feature.identity.R')
require(igraph)
# Helper function

# Homogeneous Network
# The threshold is the cut-off for edge construction from similarity.

get.homo.graph = function(sim, threshold) {
  ind = sim>threshold
  res = ind * sim
  return(res)
}

get.homo.neighbour = function(net) {
  return(apply(net, 1, function(x) sum(x>0)))
}

get.knn = function(net, k) {
  nr = nrow(net)
  res = matrix(0, nr, k)
  for (i in 1:nr) {
    ind = order(net[i,], decreasing = TRUE)[1:k]
    res[i,] = ind * as.numeric(net[i,ind] > 0)
  }
  return(res)
}

get.knn.sim = function(knn, sim) {
  nr = nrow(knn)
  nc = ncol(knn)
  res = matrix(0, nr, nc)
  for (i in 1:nr) {
    tmp = sim[i, knn[i,]]
    m = nc - length(tmp)
    if (m>0)
      tmp = c(tmp, rep(0, m))
    res[i, ] = tmp
  }
  return(res)
}

get.ave.knn = function(knn, identity) {
  nr = nrow(knn)
  nc = ncol(identity)
  res = matrix(0, nr, nc)
  for (i in 1:nr) {
    res[i,] = colMeans(identity[knn[i,],,drop = FALSE])
  }
  return(res)
}

get.w.ave.knn = function(knn, sim, identity) {
  nr = nrow(knn)
  nc = ncol(identity)
  res = matrix(0, nr, nc)
  for (i in 1:nr) {
    mat = identity[knn[i,],,drop = FALSE] * sim[i, knn[i,]]
    res[i,] = colMeans(mat)
  }
  return(res)
}

# Helper to generate graph obj from net
get.graph.from.net = function(net) {
  n = nrow(net)
  m = sum(net > 0)
  edges = matrix(0, 2, m)
  wt = rep(0,m)
  l = 1
  for (i in 1:n) {
    for (j in 1:n) {
      if (net[i,j] > 0) {
        edges[,l] = c(i,j)
        wt[l] = -log(net[i,j])
        l = l+1
      }
    }
  }
  wt = wt + 1e-7
  edges = as.vector(edges)
  graph = make_graph(edges, directed = FALSE)
  graph = set.edge.attribute(graph, "weight", value=wt)
  return(graph)
}

# TBD
# igraph::{closeness, betweenness, evcent}
get.centrality = function(graph) {
  cl = igraph::closeness(graph)
  bt = igraph::betweenness(graph, directed = FALSE)
  ev = igraph::evcent(graph, scale = FALSE)
  return(list(cl, bt, ev))
}

# TBD
# igraph::page.rank
get.pagerank = function(graph) {
  pagerank = igraph::page_rank(graph, directed = FALSE)
  return(pagerank$vector)
}

# # TBD
# # e1071::allShortestPaths
# # Use product of similarity, and use log + sum to avoid underflow
# get.floyd = function(net) {
#   log.dist = -log(net)
#   res = e1071::allShortestPaths(log.dist)
# }

get.homogeneous.feature = function(inds, ratings, sim,  k,
                                   threshold, threshold.identity) {
  # Get identity results
  res.identity = get.identity(inds, ratings, sim, threshold.identity)
  identity = res.identity
  
  # get the cut net
  net = get.homo.graph(sim, threshold)
  
  # number of neighbour
  num.nb = get.homo.neighbour(net)
  
  # k nearest neighbour, returns an n x k matrix
  knn = get.knn(net, k)
  
  # k nearest neighbour's similarity vec
  knn.sim = get.knn.sim(knn, sim)
  colnames(knn.sim) = paste0('knn.sim', 1:ncol(knn.sim))
  
  # average knn identity features
  ave.knn = get.ave.knn(knn, identity)
  colnames(ave.knn) = paste0('ave.knn', 1:ncol(ave.knn))
  
  # Weighted average knn identity features
  w.ave.knn = get.w.ave.knn(knn, sim, identity)
  colnames(w.ave.knn) = paste0('w.ave.knn', 1:ncol(w.ave.knn))
  
  # graph = get.graph.from.net(net)
  # # Centrality
  # centrality = get.centrality(graph)
  # cl = centrality[[1]]
  # bt = centrality[[2]]
  # ev = centrality[[3]]
  # ev = ev$vector
  # 
  # # Pagerank
  # pagerank = get.pagerank(graph)
  cl = 1
  bt = 1
  ev = 1
  pagerank = 1
  
  # browser()
  res = cbind(num.nb, knn.sim, ave.knn, w.ave.knn, 
              cl, bt, ev, pagerank)
  names(res) = c('num.nb', 'k.sim', paste0('k.ave.feat', 1:ncol(ave.knn)),
                 paste0('k.w.ave.feat', 1:ncol(w.ave.knn)), 'cl', 'bt', 'ev', 'pr')
  return(res)
}


# a = get.homogeneous.feature(dataset_triplets[,1], dataset_triplets[,3],
#                             drug_adj_mat, 4, 0.6, 0.6)
