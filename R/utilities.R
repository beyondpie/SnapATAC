# Find centroid of each cluster
findCentrod <- function(x, y){
	x.ls = split(data.frame(x),y);
	centroid.ls = lapply(split(data.frame(x),y), function(xx) apply(xx, 2, median))
	centroid.df = data.frame(do.call(rbind, centroid.ls))
	centroid.df$Y = names(centroid.ls);
	
	return(centroid.df);		
}


findNegCells <- function(obj, idx.pos, method=c("knn", "random", "other")){
	method=match.arg(method);
	mat.use = getGraph(obj@graph);
	ncell = nrow(obj);
	ncell.pos = length(idx.pos);
	ncell.neg = min(length(idx.pos), ncell - ncell.pos);
	
	if(method == "knn"){
		idx.neg = setdiff(seq(ncell), idx.pos);
		dx = order(Matrix::colSums(mat.use[idx.pos,idx.neg]), decreasing=TRUE)[1:ncell.neg];
		idx.neg = idx.neg[dx]
	}else if(method == "random"){
		idx.neg = setdiff(seq(ncell), idx.pos);
		nn.num = min(length(idx.pos), length(idx.neg));
		set.seed(10);
		idx.neg = sample(idx.neg, nn.num);	
	}else if(method == "other"){
		idx.neg = setdiff(seq(ncell), idx.pos);			
	}
	return(idx.neg);
}

pvalue2fdr <- function(p1, p2){
	fdr.ls <- lapply(as.list(seq(0.001, 1, by=0.001)), function(p_i){
		(length(which(p2 <= p_i))) / (length(which(p1 <= p_i)))
	})
	fdr.tab = data.frame(p=seq(0.001, 1, by=0.001), fdr=do.call(c, fdr.ls));	
	return(fdr.tab);
}


#' @importFrom parallel mclapply
splitBmat <- function(mat, id.ls, num.cores=1, tmp.folder){
	if(missing(tmp.folder)){
		stop("tmp.folder is missing")
	}else{
		if(!dir.exists(tmp.folder)){
			stop("tmp.folder does not exist");			
		}
	}
	fileList = mclapply(as.list(seq(id.ls)), function(i){
		file.name = tempfile(fileext=".rds", tmpdir=tmp.folder);
		zz <- file(description=file.name, "w"); 
		saveRDS(mat[id.ls[[i]],], file=file.name);		
		close(zz);
		file.name
	}, mc.cores=num.cores);
	return(fileList);
}


runJaccard2 <- function(
	obj1,
	obj2,
	input.mat=c("bmat", "pmat", "gmat")
){	
	calJaccard <- function(X_i, X_j){
		A = Matrix::tcrossprod(X_i, X_j);
		bi = Matrix::rowSums(X_i);
		bj = Matrix::rowSums(X_j);
		jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));
		rm(A);
		rm(bi);
		rm(bj);
		gc();	
		return(jmat);				
	}
	
	input.mat = match.arg(input.mat);	
	if(input.mat == "bmat"){
		data.use1 = obj1@bmat;
		data.use2 = obj2@bmat;
	}else if(input.mat == "pmat"){
		data.use1 = obj1@pmat;
		data.use2 = obj2@pmat;
	}else if(input.mat == "gmat"){
		data.use1 = obj1@gmat;
		data.use2 = obj2@gmat;
	}
	
	obj1@jmat@jmat = calJaccard(data.use1, data.use2);
	obj1@jmat@p1 = Matrix::rowMeans(data.use1);
	obj1@jmat@p2 = Matrix::rowMeans(data.use2);
	obj1@jmat@norm = FALSE;
	obj1@jmat@nmat = matrix(0,0,0);
	obj1@jmat@input.mat = input.mat;
	return(obj1);
}

#' Eig decomposition.
#' @param M Matrix
#' @param n_eighs integer
#' @param sym bool, is M symmetric, default is isSymmetric(M)
#' @param method character, "arpack" or "RSpectra" for decomposition, default is arpack
#' 
#' @export
eig_decomp <- function(M, n_eigs, sym = isSymmetric(M), method = "arpack") {
	n <- nrow(M)
	f <- function(x, A = NULL) as.matrix(A %*% x)
	wh <- if (sym) 'LA' else 'LM'
  if(!sym) {
    warning("Input matrix is not symmetric, double check if this is true.")
  }
	#constraints: n >= ncv > nev
  if (method == "arpack") {
    message("Use igraph::arpack for eig decomposition.")
    ar <- igraph::arpack(f, extra = M, sym = sym, options = list(
      which = wh, n = n, ncv = min(n, 4*n_eigs), nev = n_eigs + 1))
  } else {
    message("Use RSpectra::eigs for eig decomposition.")
    if (sym) {
      ar <- RSpectra::eigs_sym(A = M, k = n_eigs + 1, which = "LA")
    } else {
      ar <- RSpectra::eigs(A = M, k = n_eigs + 1, which = "LM")
    }
  }
	if (!sym) {
		ar$vectors <- Re(ar$vectors)
		ar$values  <- Re(ar$values)
	}
	if (length(dim(ar$vectors)) == 0L) {
		ar$vectors <- matrix(ar$vectors, ncol = 1L)
  }
	return(ar)
}

.normOVE <- function(p1, p2){
    pp = tcrossprod(p1, p2);
	ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
	ee = pp/(ss - pp)
	return(ee)	
}

trainRegression <- function(obj){
	row.covs = log(Matrix::rowSums(obj@bmat)+1,10);		
	row.covs.dens <- density(x = row.covs, bw = 'nrd', adjust = 1)
	sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
	idx.ds <- sort(sample(x = seq(row.covs), size = min(1000, length(row.covs)), prob = sampling_prob));
	jmat.tr = obj@jmat@jmat[idx.ds,idx.ds];
	b1.tr = obj@jmat@p1[idx.ds];
	b2.tr = obj@jmat@p2[idx.ds];
	# calculate the expected jaccard index matrix given the read depth
	emat.tr = .normOVE(b1.tr, b2.tr);
	# estimate the global scaling factor
	data = data.frame(x=emat.tr[upper.tri(emat.tr)], y=jmat.tr[upper.tri(jmat.tr)])	
	model <- lm(y ~ x + I(x^2), data);
	beta0 = as.numeric(model$coefficients)[1]
	beta1 = as.numeric(model$coefficients)[2]
	beta2 = as.numeric(model$coefficients)[3]
	obj@regModel = c(beta0, beta1, beta2);
	rm(jmat.tr);
	rm(emat.tr);
	rm(data);
	rm(model);
	rm(row.covs);
	return(obj)
}

normJaccard <- function(obj, beta0, beta1, beta2){
	b1.te = obj@jmat@p1;
	b2.te = obj@jmat@p2;
	jmat.te = obj@jmat@jmat;
	emat.te = .normOVE(b1.te, b2.te);
	preds = beta0 + beta1 * emat.te + beta2 * (emat.te ** 2);
	nmat.te = jmat.te/preds;
	obj@jmat@nmat = nmat.te;
	obj@jmat@norm = TRUE;
	rm(jmat.te);
	rm(emat.te);
	rm(nmat.te);
	rm(preds);
	return(obj);
}

runEigDecomp <- function(obj, num.eigs, method = "arpack"){
	nmat <- obj@jmat@nmat
	diag(nmat) <- 0
	norm_p1 <- nmat
	d_norm1 <- Matrix::rowSums(norm_p1)
	d_rot1 <- Diagonal(x = d_norm1 ^ -.5)
	transitions <- as.matrix(d_rot1 %*% norm_p1 %*% d_rot1)
	diag(transitions) <- 0
	eig_transitions <- eig_decomp(transitions, num.eigs, method = method)
	obj@smat@dmat <- eig_transitions$vectors[,2:(num.eigs+1)]
	obj@smat@sdev <- eig_transitions$value[2:(num.eigs+1)]
	return(obj)
}

runEigDecompExd <- function(obj1, obj2){
	# obj1 is reference cells
	# obj2 is query cells
	nmat1 = obj1@jmat@nmat;
	diag(nmat1) = 0;
	norm_p1 = nmat1
	d_norm1 <- Matrix::rowSums(norm_p1);
	d_rot1 <- Diagonal(x = d_norm1 ^ -.5);
	trans_p = obj2@jmat@nmat;
	norm_p = trans_p
	d_new <- Matrix::rowSums(trans_p, na.rm = TRUE)
	d_norm_new <- Matrix::rowSums(norm_p)
	d_rot_new <- Diagonal(x = d_norm_new ^ -.5)
	M_new <- d_rot_new %*% norm_p %*% d_rot1
	obj2@smat@dmat = as.matrix(t(t(M_new %*% obj1@smat@dmat) / obj1@smat@sdev));
	obj2@smat@sdev = obj1@smat@sdev;
	rm(d_norm1);
	rm(d_rot1);
	rm(trans_p);
	rm(d_norm_new);
	rm(M_new);
	gc();
	return(obj2);
}


#' Check if all the snap files exist.
#' @param snapFiles list of string
#' @return bool, TRUE if all of them exist; FALSE if not.
#' @export
checkSnapFilesExist <- function(snapFiles) {
  fileList <- as.list(unique(snapFiles))
  message("check if all the snap files exist.")
  if(any(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)){
    idx <- which(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)
    message(fileList[idx], " doest not exist.")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Check if all the files are snap files.
#' @param snapFiles list of string
#' @param do.par bool, if do in parallel, default is FALSE.
#' @param num.cores integer, cores to be used, default is 1
#' @return bool, TRUE if all of them are snap files; FALSE if no
#' @export
checkFilesSnapFiles <- function(snapFiles, do.par = FALSE, num.cores = 1) {
  fileList <- as.list(unique(snapFiles))
  message("check is all the files are snap files.")
  if(do.par) {
    isSnaps <- unlist(parallel::mclapply(fileList, isSnapFile, mc.cores = num.cores))
  } else {
    isSnaps <- unlist(lapply(fileList, isSnapFile))
  }
  if(any(isSnaps == FALSE)){
    idx <- which(isSnaps == FALSE)
    message(fileList[idx], " does not exist.")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Check if all snap files have a given session.
#' @param snapFiles list of string
#' @param session character
#' @param do.par bool, if do in parallel, default is FALSE.
#' @param num.cores integer, cores to be used, default is 1
#' @return bool, TRUE if all snap files have the session; FALSE if no.
#' @export
checkSnapFileSession <- function(snapFiles, session, do.par = FALSE, num.cores = 1) {
  fileList <- as.list(unique(snapFiles))
  message("check if session ", session, " exist.")
  if(do.par) {
    hasAMs <- unlist(
      parallel::mclapply(fileList,function(x){"AM" %in% rhdf5::h5ls(x, recursive=1)$name},
                         mc.cores = num.cores))
  } else {
    hasAMs <- unlist(lapply(fileList, function(x) {"AM" %in% rhdf5::h5ls(x, recursive=1)$name}))
  }
	if(any(hasAMs == FALSE)){
    idx <- which(hasAMs == FALSE)
    message(fileList[idx], " do not have ", session, "session.")
    return(FALSE)
	} else {
    return(TRUE)
  }
}
