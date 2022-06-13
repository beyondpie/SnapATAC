#' @include snap-class.R
NULL
#' K Nearest Neighbour Graph Construction
#'
#' Constructs a K Nearest Neighbor (SNN) Graph from a snap object. 
#' The k-nearest neighbors of each cell were identified and used to create
#' a KNN graph. 
#' 
#' Using the selected significant principal components (PCs), we next calculated pairwise 
#' Euclidean distance between every two cells, using this distance, we created a k-nearest 
#' neighbor graph in which every cell is represented as a node and edges are drawn between 
#' cells within k nearest neighbors. Edge weight between any two cells can be refined by shared 
#' overlap in their local neighborhoods using Jaccard similarity (snn).
#'
#' @param obj A snap object
#' @param eigs.dims A vector of the dimensions to use in construction of the KNN graph.
#' @param weight.by.lambda Weight the cell embeddings by the sd of each PC
#' @param k K for the k-nearest neighbor algorithm.
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN.
#' default of 0.0 implies exact nearest neighbor search
#' @param save.knn Default is to store the KNN in object@@kmat. Setting
#' to FALSE can be used together with a provided filename to only write the KNN
#' out as an edge file to disk. This is compatible with runCluster.
#' @param filename Write KNN directly to file named here as an edge list
#' compatible with runCluster.
#' @param snn Setting to TRUE can convert KNN graph into a SNN graph.
#' @param snn.prune Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap. Any edges with values less than or 
#' equal to this will be set to 0 and removed from the SNN graph. Essentially 
#' sets the strigency of pruning (0 --- no pruning, 1 --- prune everything).
#' @param method character, "RANN" or "Annoy", the later is faster but approximated
#' method.
#' 
#' @examples
#' data(demo.sp);
#' demo.sp = runKNN(obj=demo.sp, eigs.dims=1:5, k=15, snn=FALSE, save.knn=FALSE);
#' 
#' @importFrom RANN nn2
#' @importFrom igraph similarity graph_from_edgelist E
#' @importFrom Matrix sparseMatrix
#' @importFrom plyr count
#' @importFrom methods as
#' @import Matrix
#' @return Returns the object with object@kmat filled
#' @export
runKNN <- function(obj, eigs.dims, weight.by.lambda, k, nn.eps, save.knn,
                   filename, snn, snn.prune, method) {
  UseMethod("runKNN", obj);
}

#' @export
runKNN.default <- function(
  obj,
  eigs.dims,
  weight.by.lambda = FALSE,
  k = 15,
  nn.eps = 0,
  save.knn = FALSE,
  filename = NULL,
  snn = FALSE,
  snn.prune = 1/15,
  method = "RANN"
){
	
	cat("Epoch: checking input parameters\n", file = stderr())
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap obj")
		}
	}
	
	if(!(isDimReductComplete(obj@smat))){
		stop("dimentionality reduction is not complete, run 'runDimReduct' first")
	}
	
	ncell = nrow(obj);
	nvar = dimReductDim(obj@smat);
	
	if(missing(eigs.dims)){
		stop("eigs.dims is missing")
	}else{
		if(is.null(eigs.dims)){
			eigs.dims=1:nvar;	
		}else{
			if(any(eigs.dims > nvar) ){
				stop("'eigs.dims' exceeds PCA dimentions number");
			}		
		}
	}
	
	if(save.knn){
		if(is.null(filename)){
			stop("save.knn is TRUE but filename is NULL")
		}else{			
			if(!file.create(filename)){
				stop("fail to create filename")
			}			
		}
	}
	
	if(!is.logical(weight.by.lambda)){
		stop("weight.by.lambda must be a logical variable")
	}
	
	data.use = weightDimReduct(obj@smat, eigs.dims, weight.by.lambda);
	
    if (ncell < k) {
      warning("k set larger than number of cells. Setting k to number of cells - 1.")
      k <- ncell - 1
    }
	
	if(is.na(as.integer(k))){
		stop("k must be an integer")
	}else{
		if(k < 10 || k > 50){
			warning("too small or too large k, recommend to set k within range [10 - 50]")
		}
	}
	
	cat("Epoch: computing nearest neighbor graph\n", file = stderr())
	
	# exclude self neibours
  if (method == "RANN") {
    message("Using RANN to get the KNN result.")
    nn.ranked <- nn2(
      data = data.use,
      k = k,
      searchtype = 'standard',
      eps = nn.eps)$nn.idx
  } else {
    message("Using Annoy to get the KNN result.")
    nn.ranked <- AnnoyNN(data = data.use, k = k)$nn.idx
  }

	j <- as.numeric(x = t(x = nn.ranked))
	i <- ((1:length(x = j)) - 1) %/% k + 1	
	edgeList <- data.frame(i, j, 1)
		
	if(snn){
		cat("Epoch: converting knn graph into snn graph\n", file = stderr())	
		g = graph_from_edgelist(as.matrix(edgeList[,c(1,2)]), directed=FALSE);
		adj = as(similarity(g), "sparseMatrix");
		i = adj@i+1;
		j = findInterval(seq(adj@x)-1,adj@p[-1])+1;
		w = adj@x;
		idx = which(w >= snn.prune);
		edgeList = data.frame(i[idx], j[idx], w[idx]);
	}
	
	if(save.knn){
		cat("Epoch: writing resulting graph into a file\n", file = stderr())
		writeEdgeListToFile(edgeList, filename);
		obj@graph = newKgraph(file=filename, k=k, snn=snn, snn.prune=snn.prune);
	}else{
		kmat = Matrix(0, ncell, ncell, sparse=TRUE);
		kmat = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3]);
		obj@graph = newKgraph(mat=kmat, k=k, snn=snn, snn.prune=snn.prune);
	}
	gc();
	return(obj);
} 



####################
## From Seurat
####################
# Run annoy
#
#' @param data Data to build the index with
#' @param query A set of data to be queried against data
#' @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
#' "hamming"
#' @param n.trees More trees gives higher precision when querying
#' @param k Number of neighbors
#' @param search.k During the query it will inspect up to search_k nodes which
#' gives you a run-time tradeoff between better accuracy and speed.
#' @param include.distance Include the corresponding distances
#' @param index optional index object, will be recomputed if not provided
AnnoyNN <- function(data,
                    query = data,
                    metric = "euclidean",
                    n.trees = 50,
                    k,
                    search.k = -1,
                    include.distance = TRUE,
                    index = NULL
) {
  idx <- index %||% AnnoyBuildIndex(
    data = data,
    metric = metric,
    n.trees = n.trees)
  nn <- AnnoySearch(
    index = idx,
    query = query,
    k = k,
    search.k = search.k,
    include.distance = include.distance)
  nn$idx <- idx
  nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
  return(nn)
}

# Build the annoy index
#
#' @param data Data to build the index with
#' @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
#' "hamming"
#' @param n.trees More trees gives higher precision when querying
#
#' @importFrom RcppAnnoy AnnoyEuclidean AnnoyAngular AnnoyManhattan AnnoyHamming
AnnoyBuildIndex <- function(data, metric = "euclidean", n.trees = 50) {
  f <- ncol(x = data)
  a <- switch(
    EXPR = metric,
    "euclidean" =  new(Class = RcppAnnoy::AnnoyEuclidean, f),
    "cosine" = new(Class = RcppAnnoy::AnnoyAngular, f),
    "manhattan" = new(Class = RcppAnnoy::AnnoyManhattan, f),
    "hamming" = new(Class = RcppAnnoy::AnnoyHamming, f),
    stop ("Invalid metric")
  )
  for (ii in seq(nrow(x = data))) {
    a$addItem(ii - 1, data[ii, ])
  }
  a$build(n.trees)
  return(a)
}

# Search an Annoy approximate nearest neighbor index
#
# @param Annoy index, built with AnnoyBuildIndex
# @param query A set of data to be queried against the index
# @param k Number of neighbors
# @param search.k During the query it will inspect up to search_k nodes which
# gives you a run-time tradeoff between better accuracy and speed.
# @param include.distance Include the corresponding distances in the result
#
# @return A list with 'nn.idx' (for each element in 'query', the index of the
# nearest k elements in the index) and 'nn.dists' (the distances of the nearest
# k elements)
#
#' @importFrom future plan
#' @importFrom future.apply future_lapply
AnnoySearch <- function(index, query, k, search.k = -1, include.distance = TRUE) {
  n <- nrow(x = query)
  idx <- matrix(nrow = n,  ncol = k)
  dist <- matrix(nrow = n, ncol = k)
  convert <- methods::is(index, "Rcpp_AnnoyAngular")
  if (!inherits(x = plan(), what = "multicore")) {
    oplan <- plan(strategy = "sequential")
    on.exit(plan(oplan), add = TRUE)
  }
  res <- future_lapply(X = 1:n, FUN = function(x) {
    res <- index$getNNsByVectorList(query[x, ], k, search.k, include.distance)
    # Convert from Angular to Cosine distance
    if (convert) {
      res$dist <- 0.5 * (res$dist * res$dist)
    }
    list(res$item + 1, res$distance)
  })
  for (i in 1:n) {
    idx[i, ] <- res[[i]][[1]]
    if (include.distance) {
      dist[i, ] <- res[[i]][[2]]
    }
  }
  return(list(nn.idx = idx, nn.dists = dist))
}

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}
