#' Explode dimension of a 2d matrix to 3d array
#'
#' @description
#' Turns a 2d matrix *A* of size *s* x *c* into a 3d array by replicating
#' rows of the matrix *A* either vertically or horizontally.
#'
#' This function takes a matrix *A*
#' ```
#' A = [[a_11, ..., a_1c],
#'      [... , ..., ... ]
#'      [a_s1, ..., a_sc]]
#' ```
#'
#' transforms each row `r` to form a separate 2d matrix *Ar*
#' based on the `byrow` parameters, and returns a 3d array
#' `[A1, ..., As]`.
#'
#' @param A
#' a `matrix`.
#'
#' @param byrow
#' states if the rows of the matrix should be replicated
#' _vertically_ (`byrow=FALSE`)
#' or _horizontally_ (`byrow=TRUE`).
#' Note the default behaviour is to replicate _horizontally_.
#'
#' _vertical_ replication means each *Ar* looks like
#' ```
#' A_r = [[a_r1, ..., a_r1],
#'        [... , ..., ... ],
#'        [a_rc, ..., a_rc]],
#' ```
#' while _horizontal_ replication means each *Ar* looks like
#' ```
#' A_r = [[a_r1, ..., a_rc],
#'        [... , ..., ... ],
#'        [a_r1, ..., a_rc]].
#' ```
#'
#' @return a 3d `array` with appropriately replicated rows.
#'
#' @examples
#' A <- matrix(c(1, 2, 3, 4), ncol=2)
#' # > A
#' #      [,1] [,2]
#' # [1,]    1    3
#' # [2,]    2    4
#' B <- explode_dimension(A)
#' # > B
#' # , , 1
#' #
#' #      [,1] [,2]
#' # [1,]    1    1
#' # [2,]    3    3
#' #
#' # , , 2
#' #
#' #      [,1] [,2]
#' # [1,]    2    2
#' # [2,]    4    4
explode_dimension <- function(A, byrow=FALSE) {
  stopifnot(is.matrix(A))
  apply(A, 1,
        function(x) matrix(x, ncol(A), ncol(A), byrow=byrow), simplify=FALSE) |>
    unlist() |>
    array(dim=c(ncol(A), ncol(A), nrow(A)))
}

update_cluster_centers <-
  function(
    U,
    X,
    alpha=NULL,
    F_=NULL
  ) {
    if (is.null(alpha)) {
      V <-
        t(sweep(
          t(X) %*% U^2,
          2,
          colSums(U^2),
          "/"
        ))
    } else {
      UF <- alpha * (U-F_)^2
      i_indices <- which(rowSums(F_) != 0)
      j_indices <- 1:nrow(F_)
      h_indices <- setdiff(j_indices, i_indices)
      UF[h_indices,] <- 0.
      Phi <- U^2 + UF

      V <- t(sweep(
        t(X) %*% Phi,
        2,
        colSums(Phi),
        "/"
      ))
    }

    return(V)
  }

#' Calculating distances between two matrices
#'
#' A base R approach to calculate distances between
#' rows of two matrices A and B that are of the same ncol.
#'
#' TODO it should be moved to tests to compare with other, externally
#' imported methods because this function turned out to be very slow.
#'
#' @param X
#' a matrix of dimension (N, p)
#'
#' @param V
#' a matrix of dimension (C, p)
#'
#' @returns a matrix of dimension (N, C)
#'
calculate_distances <-
  function(
    X,
    V
  ) {
    process_distance <- function(x, y) {
      output <- as.matrix(dist(rbind(x, y)))
      output[2:nrow(output), 1]
    }

    D_ <- apply(X, 1, function(x, y) process_distance(x, y), V, simplify=FALSE)
    D <- t(do.call(cbind, D_))

    return(D)
  }

#' Estimated U matrix with memberships
#'
#' The trick is to use 3d arrays for efficient calculations.
#' Let's take a D^2 matrix
#' ```
#' D62 = [[5, 145],
#'        [25, 85],
#'        [61, 41]]
#' ```
#'
#' `D_nominator` will be a 3d array of dimensions (2, 2, 3):
#' it contains 3 matrices 2x2.
#' Let's take a single 2x2 matrix from the array `D_nominator[,,1]`
#' ```
#'      [,1] [,2]
#' [1,]    5    5
#' [2,]  145  145
#' ```
#' It is effectively
#' ```
#'           [,1]      [,2]
#' [1,]    d^2_11    d^2_11
#' [2,]    d^2_12    d^2_12
#' ```
#' Let's now take `D_denominator[,,1]`
#' ```
#'      [,1] [,2]
#' [1,]    5  145
#' [2,]    5  145
#' ```
#' It is effectively
#' ```
#'           [,1]      [,2]
#' [1,]    d^2_11    d^2_12
#' [2,]    d^2_11    d^2_12
#' ```
#' If we now divide `D_nominator[,,1]` by `D_denominator[,,1]`, we get
#' ```
#'          [,1]     [,2]
#' [1,]    5/5    5/145
#' [2,]    145/5  145/145
#' ```
#' which is
#'                  [,1]             [,2]
#' [1,]    d^2_11/d^2_11    d^2_11/d^2_12
#' [2,]    d^2_12/d^2_11    d^2_12/d^2_12
#' The `rowSums` of the above will give us a vector
#' ```
#' [1] e_11 e_12
#' ```
#' *almost* the first row of evidence matrix, which is the inverse of the above
#' i.e. a vector `[1] 1/e_11 1/e_12`.
#'
#' The vectorized operations on entire `D_nominator` and `D_denominator`
#' follow the above logic.
#'
#' @param X
#' a matrix *X* of dimension (N, p) containing predictor variables.
#'
#' @param V
#' a prototypes matrix of dimension (C, p)
#'
#' @param F_
#' the supervision  binary matrix of the same dimension as *U*.
#'
#' @param alpha
#' the scaling factor, a floating point > 0.
#'
#' @param fun.distances
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
update_memberships <-
  function(
    X,
    V,
    F_,
    alpha,
    fun.distances, metric
  ) {
    D <- fun.distances(X, V, metric)^2

    D_nominator <- explode_dimension(D)
    D_denominator <- explode_dimension(D, byrow=TRUE)

    E_reciprocal <- t(apply(D_nominator/D_denominator, c(3), rowSums))
    E <- 1/E_reciprocal

    if (is.null(alpha)) {
      return(E)
    } else {
      i_indices <- which(rowSums(F_) != 0)
      M <- matrix(1, nrow(F_), ncol(F_))
      M[i_indices, ] <- 1/(1+alpha)

      ALB = F_*(alpha/(1+alpha))

      return(M*E + ALB)
    }
  }

#' Semi-Supervised Fuzzy C-Means model.
#'
#' @description
#' If *alpha* and *F_* are not supplied (their default values are `NULL`),
#' then a regular unsupervised Fuzzy C-Means algorithm is fitted.
#'
#' @param X
#' a matrix *X* with predictor variables.
#'
#' @param C
#' a number of clusters to find.
#'
#' @param U
#' optionally: a first memberships matrix to initialize the algorithm.
#' Used mainly for reproducibility to compare calculations with other packages
#' (e.g. in Python).
#'
#' @param fun.distances
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
#' @param alpha
#' the scaling factor, a floating point > 0.
#'
#' @param F_
#' the supervision  binary matrix of the same dimension as *U*.
#'
#' @export
#'
#' @examples
#' # simulate 50 obs from N2((5, 5)) and 50 obs from N2((7, 7))
#' library(MASS)
#' X <- rbind(
#'   MASS::mvrnorm(50, mu=c(5, 8), Sigma=matrix(c(3, 0, 0, 3), ncol=2)),
#'   MASS::mvrnorm(50, mu=c(7, 10), Sigma=matrix(c(3, 0, 0, 3), ncol=2))
#'   )
#'
#' # simulate supervision for 10% of each class
#' F_ <- matrix(0, nrow=100, ncol=2)
#' F_[sample(1:50, 10), 1] <- 1
#' F_[sample(51:100, 10), 2] <- 1
#'
#' model <- SSFCM(X=X, C=2)
#' model.ss <- SSFCM(X=X, C=2, alpha=1, F_=F_)
#'
#' acc.unsupervised.1 <- sum(apply(model$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
#' acc.unsupervised.2 <- sum(apply(model$U, 1, which.max) == c(rep(2, 50), rep(1, 50)))
#' acc.supervised <- sum(apply(model.ss$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
#'
SSFCM <- function(
    X,
    C,
    U=NULL,
    max_iter=200,
    conv_criterion=1e-4,
    fun.distances=rdist::cdist,
    alpha=NULL,
    F_=NULL
) {
  # random U if not supplied
  if (is.null(U)) {
    U <- matrix(runif(nrow(X)*C), ncol=C)
  }

  # normalize U
  U <- t(apply(U, 1, function(x) x / sum(x)))

  counter = 0

  # calculations loop
  for (iter in 1:max_iter) {
    counter <- counter + 1
    U_previous_iter <- U

    V <- update_cluster_centers(
      U=U_previous_iter,
      X=X,
      alpha=alpha,
      F_=F_)

    U <- update_memberships(
      X=X,
      V=V,
      F_=F_,
      alpha=alpha,
      fun.distances=fun.distances)

    conv_iter <- base::norm(U - U_previous_iter, type="F")

    if (conv_iter < conv_criterion) {
      break
    }
  }
  return(list(U=U,V=V, counter=counter))
}



SSFCM_1 <- function(
    X,
    C,
    U=NULL,
    max_iter=200,
    conv_criterion=1e-4,
    fun.distances=rdist::cdist,
    alpha=NULL,
    F_=NULL,
    metric= "euclidean")
  {
  # random U if not supplied
  #print(paste0('C', str(C)))
  
  if (is.null(U)) {
    U <- matrix(runif(nrow(X)*C), ncol=C)
  }
  
  # normalize U
  U <- t(apply(U, 1, function(x) x / sum(x)))
  
  counter = 0
  
  # calculations loop
  for (iter in 1:max_iter) {
    counter <- counter + 1
    U_previous_iter <- U
    
    V <- update_cluster_centers(
      U=U_previous_iter,
      X=X,
      alpha=alpha,
      F_=F_)
    
    U <- update_memberships(
      X=X,
      V=V,
      F_=F_,
      alpha=alpha,
      fun.distances=fun.distances, metric="euclidean")
    
    conv_iter <- base::norm(U - U_previous_iter, type="F")
    
    if (conv_iter < conv_criterion) {
      break
    }
  }
  return(list(U=U,V=V, counter=counter))
}


##########3

######################
  
  calculate_reconstruction_error <- function(X, U,V) {
    N_t <- nrow(X)
    K <- ncol(U)
    n <- ncol(X)
    x_hat <- matrix(0,nrow=N_t,ncol=ncol(V))
    q = sqrt(sum(X)^2)
    X$pred <- max.col(U)
    ability<-matrix(0,nrow=N_t,ncol=K)
    for (j in 1:N_t) {
      x_hat[j,] <- rowSums(U[j, ]^2 * t(V))/sum(U[j,]^2)
      for (k in 1:K) {
        if (X[j,'pred']==k){
          ability[j,k] <- 1/q * sqrt(sum((X[j, 1:2] - x_hat[j, ])^2))
        }
      }
    }
    errors <- list()
    max_column <- which(ability == max(ability), arr.ind = TRUE)[, "col"]
    errors$error <- max(ability, na.rm = TRUE)
    errors$lowest_ability <- max_column
    return(errors)
  } 

splitting <- function(X){
  #print(dim(X))
  clust<-FKM(X,k=2)
  return(list(U=clust$U,V=clust$H))
}
  
DISSFCM <- function(chunks, C, F_s,selected_columns, error=0.5,max_iter = 200, conv_criterion = 1e-4, fun.distances = rdist::cdist) {
  T <- length(chunks)
  reconstruction_error <- list()
  reconstruction_error[[1]] <- c(error)
  U_s <- list()
  V_s <- list()
  U_new_s <- list()
  for (t in 1:T){
    chunk <- chunks[[t]]
    if (dim(chunk)[1]>0) {
      if (t==1){
        result <- SSFCM_1(chunk[,selected_columns], C = C,F_=F_s[[t]], metric="euclidean")}
      else if (t>1){ 
        non_zero_columns <- colSums(F_s[[t]]) > 0
        if (all(!non_zero_columns)) {
          #C <- ncol(F_s[[t]])
          F_filtered <- F_s[[t]]
        } else {
         # C <- ncol(F_s[[t]][, non_zero_columns])
          F_filtered <- F_s[[t]][, non_zero_columns]
        }
        
        if (!is.null(ncol(F_filtered))) {
          result <- SSFCM_1(chunk[, selected_columns], C = ncol(F_filtered), F_ = F_filtered, metric = 'euclidean')
        } else {
          print("C jest NULL")
        }
      }
      U <- result$U
      V <- result$V
      chunk$pred <- max.col(U)
      U_s[[t]] <- U
      V_s[[t]] <- V
      check <- checking_1(chunks,C,U_s,V_s, selected_columns,reconstruction_error,t,F_s)
      reconstruction_error <- check$error
      U_s[[t]]<-check$U
      V_s[[t]] <- check$V
      mapping<- mapClusters(F_s=F_s,U_s=U_s,t=t,t2=t+1)
      U_new_s[[t]] <- mapping$U_new
      F_s[[t]] <- mapping$F_new
      } else {
      break
    }
  }
  return(list(U=U_new_s,V=V_s,reconstruction_error=reconstruction_error))
}

checking_1 <- function(chunks, C, U_s, V_s,selected_columns,reconstruction_error,t,F_s){
  U <- U_s[[t]]
  V <- V_s[[t]]
  chunk <- chunks[[t]]
  chunk$pred <- max.col(U)
  reconstruction <- calculate_reconstruction_error(chunk[,selected_columns], U, V)
  reconstruction_error[[t+1]] <- c(reconstruction$error)
  iteration <- 1
  while (reconstruction_error[[t+1]][[iteration]] > tail(reconstruction_error[[t]],1)){
    if (iteration > 1){
      result <- SSFCM_1(chunk[,selected_columns], C = ncol(U), F_=F_new,metric="euclidean")
      chunks$pred <- max.col(result$U)
      U = result$U
      V = result$V
    }
    k_star <- reconstruction$lowest_ability
    print(paste("w chunku",t,"Split cluster", k_star))
    in_cluster <- which(chunk$pred==k_star)
    if (length(in_cluster)<2){
      paste0('in chunk',t,'in cluster',k_star,'we dont have enough observation')
      break
    }
    split_1 <- splitting(chunk[in_cluster,selected_columns])
    U_1 <- matrix(0,ncol=2,nrow=nrow(U))
    V_1 <- split_1$V
    U_1[in_cluster,] <- split_1$U
    U  <- cbind(U[,-k_star], U_1 * U[,k_star])
    V  <- rbind(V[ -k_star,], V_1)
    iteration <- iteration +1
    U_s[[t]] <- U
    V_s[[t]] <- V
    mapping<- mapClusters(F_s=F_s,U_s=U_s,t=t,t2=t)
    F_new <- mapping$F_new
    reconstruction <- calculate_reconstruction_error(chunk[,selected_columns], U, V)
    reconstruction_error[[t+1]] <- append(reconstruction_error[[t+1]],reconstruction$error)
    if (iteration >= 10) {
      break 
    }
  }
  return(list(U=U, V=V,error=reconstruction_error))
}
  

mapClusters <- function(F_s,U_s,t,t2){
  n_classes <- ncol(F_s[[t]])
  n_clusters <- ncol(U_s[[t]])
  #print(n_classes)
  #print(n_clusters)
  F_ <- F_s[[t]]
  U <- U_s[[t]]
  match_matrix <- matrix(0, nrow = n_clusters, ncol = n_classes)
  for (i in 1:n_clusters) {
    for (j in 1:n_classes) {
      match_matrix[i, j] <- sum(F_[, j] == 1 & max.col(U) == i)
    }
  }
  best_match <- apply(match_matrix, 1, which.max)
  best_match
  U_new <- rep(0, nrow(U))
  n_obs <- nrow(F_)
  for (i in 1:nrow(U)) {
    cluster_index <- which.max(U[i, ])
    U_new[i] <- best_match[cluster_index]
  }
  if (length(F_s) >= t2 && !is.null(F_s[[t2]])) {
    if (is.null(F_s[[t2]])) {
      F_new <- NULL
    } else {
      F_new <- matrix(0, nrow = nrow(F_s[[t2]]), ncol = n_clusters)
      klasa <- c()
      skupienie <- list()
      for (i in 1:nrow(F_s[[t2]])) {
        if (all(F_s[[t2]][i, ] == 0)) {
          klasa[i] <- NA
          skupienie[[i]] <- NA
        } else {
          klasa[i] <- max.col(F_s[[t2]], ties.method = "first")[i]
          skupienie[[i]] <- which(best_match == klasa[i])
          F_new[i, skupienie[[i]]] <- 1
        }
      }
    }
  } else {
    F_new <- NULL
  }
  
  return(list(U_new=U_new, F_new = F_new))
}

assfcm <- function(chunks, C, F_s,selected_columns,threshold,error=0.5, metric='euclidean', max_iter = 200, conv_criterion = 1e-4, fun.distances = rdist::cdist) {
  T <- length(chunks)
  reconstruction_error <- list()
  reconstruction_error[[1]] <- error
  U_s <- list()
  V_s <- list()
  U_new_s <- list()
  for (t in 1:T){
    chunk <- chunks[[t]]
    if (dim(chunk)[1]>0) {
      if (t==1){
        result <- SSFCM_1(chunk[,selected_columns], C =ncol(F_s[[t]]),F_=F_s[[t]],metric='euclidean')}
      else if (t>1){ 
        non_zero_columns <- colSums(F_s[[t]]) > 0
        if (all(!non_zero_columns)) {
          C <- ncol(F_s[[t]])
          F_filtered <- F_s[[t]]
        } else {
          C <- ncol(F_s[[t]][, non_zero_columns])
          F_filtered <- F_s[[t]][, non_zero_columns]
        }
        if (!is.null(ncol(F_filtered))) {
          #print(paste0('C', C))
          result <- SSFCM_1(chunk[, selected_columns], C = ncol(F_filtered), F_ = F_filtered, metric = 'euclidean')
        } else {
          print("C jest NULL")
        }
      }
      U <- result$U
      V <- result$V
      chunk$pred <- max.col(U)
      U_s[[t]] <- U
      V_s[[t]] <- V
      check <- checking_2(chunks,U_s,V_s, selected_columns,reconstruction_error=reconstruction_error,t=t,F_s=F_s,threshold=threshold)
      reconstruction_error[t+1] <- check$error
      U_s[[t]]<-check$U
      V_s[[t]] <- check$V
      max_distance <- check$max_distance
      dim(F_s[[t]])
      mapping<- mapClusters(F_s=F_s,U_s,t,t2=t+1)
      U_new_s[[t]] <- mapping$U_new
      F_s[[t]] <- mapping$F_new
    } else {
      break
    }
  }
  return(list(U=U_new_s,V=V_s,reconstruction_error=reconstruction_error,max_distance=max_distance))
}

checking_2 <- function(chunks, U_s, V_s, selected_columns, reconstruction_error, t,F_s, threshold, metric = 'euclidean', fun.distances = rdist::cdist) {
  chunk <- chunks[[t]]
  #print(paste("F na początku checking",dim(F_s)))
  V <- V_s[[t]]
  U <- U_s[[t]]
  chunk$pred <- max.col(U)
  reconstruction <- calculate_reconstruction_error(chunk[, selected_columns], U, V)
  reconstruction_error[[t + 1]] <- c(reconstruction$error)
  iteration <- 1
  max_distance <- 0
  for (chunk1 in chunks) {
    data_points <- chunk1[, selected_columns] 
    for (i in 1:(nrow(data_points) - 1)) {
      for (j in (i + 1):nrow(data_points)) {
        distance <- fun.distances(data_points[i, ], data_points[j, ], metric)
        if (distance > max_distance) {
          max_distance <- distance
        }
      }
    }
  }
  print(paste('max distance',max_distance))
  
  k_stars <- c()
  threshold_1 <- threshold * max_distance
  print(paste('threshold',threshold_1))
  while (reconstruction_error[[t+1]][[iteration]] > tail(reconstruction_error[[t]],1)) {
    if (iteration > 1) {
      result <- SSFCM_1(chunk[, selected_columns], C = ncol(U), F_=F_new, metric = 'euclidean')
      chunk$pred <- max.col(result$U)
      U <- result$U
      V <- result$V
    }
    
    k_star <- reconstruction$lowest_ability
    if (k_star %in% k_stars) {
      print('k star się powtarza')
      break } else {
      k_stars[iteration] <- k_star
      print(paste("w chunku", t, "Split cluster", k_star))
      in_cluster <- which(chunk$pred == k_star)
      split_1 <- splitting(chunk[in_cluster, selected_columns])
      U_1 <- matrix(0, ncol = 2, nrow = nrow(U))
      V_1 <- split_1$V
      U_1[in_cluster, ] <- split_1$U
      U <- cbind(U[, -k_star], U_1 * U[, k_star])
      U_s[[t]]<- U
      V <- rbind(V[-k_star, ], V_1)
      rownames(V) <- 1:nrow(V)
    
      while (nrow(V) > 2) {
        #print(nrow(V))
        min_pair_info <- find_min_distance_pair(V, metric)
        min_pair <- min_pair_info$pair
        min_distance <- min_pair_info$distance
        print(paste('min distance',min_distance))
        if (min_distance >= threshold_1) {
          break
        }
        
        i <- min_pair[1]
        j <- min_pair[2]
        
        if (i <= nrow(V) && j <= nrow(V)) {
          new_V <- (V[i, ] + V[j, ]) / 2
          V <- rbind(V[-c(i, j), ], new_V)
          U <- cbind(U[, -c(i, j)], U[, i] + U[, j])
          cat("Usunięto parę:", i, j, "\n")
          #print(dim(V))
          #print(dim(U))
        } else {
          print(paste("Błąd: i =", i, "lub j =", j, "jest poza zakresem liczby wierszy V"))
        }
    }
    
    
    U_s[[t]] <- U
    V_s[[t]] <- V
    iteration <- iteration + 1
    reconstruction <- calculate_reconstruction_error(chunk[, selected_columns], U, V)
    reconstruction_error[[t + 1]] <- append(reconstruction_error[[t+1]],reconstruction$error)
    mapping<- mapClusters(F_s,U_s,t,t)
    F_new <- mapping$F_new
    }
  }
  
  return(list(U = U, V = V, error = reconstruction_error[t + 1], max_distance=max_distance))
}

find_min_distance_pair <- function(V, metric = 'euclidean', fun.distances = rdist::cdist) {
  min_distance <- Inf
  min_pair <- NULL
  for (i in 1:(nrow(V) - 1)) {
    for (j in (i + 1):nrow(V)) {
      distance <- fun.distances(t(matrix(V[i, ])), t(matrix(V[j, ])), metric)
      if (distance < min_distance) {
        min_distance <- distance
        min_pair <- c(i, j)
      }
    }
  }
  return(list(pair = min_pair, distance = min_distance))
}





