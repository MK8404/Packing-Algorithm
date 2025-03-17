# matrix_generation.R
# ======================================
# Functions for generating different types of matrices
# ======================================

#' Generate a band matrix in sparse representation
#' @param l Integer, the bandwidth of the matrix
#' @param d Integer, the size of the matrix
#' @return A list representing the sparse band matrix
bamat <- function(l, d) {
  if (!is.numeric(l) || l %% 1 != 0 || l <= 0 || !is.numeric(d) || d %% 1 != 0 || d <= 0) {
    stop("Error: The parameters 'l' and 'd' must be natural numbers.")
  }  
  if (l>=d-1) {
    warning("Since the condition l>=d-1 holds, all elements of the generated matrix are non-zero")
  }
  A <- list()
  
  for (j in 1:d) {
    # Create a list of indices representing the band structure
    A <- append(A, list(max(1,j-l):min(d,j+l)))
  }
  return(A)
}


#' Convert sparse representation to matrix form
#' @param A List representing the sparse band matrix
#' @return Matrix representation of the List A
sp2mat <- function(A) {

  n <- length(A)
  B <- matrix(0, ncol = n, nrow = n)
  for (j in 1:n) {
    # Populate the matrix using the List A
    B[j, A[[j]]] <- 1
  }
  return(B)
}


#' Convert matrix form to sparse representation
#' @param B Matrix
#' @return Sparse representation as a list
mat2sp <- function(B) {
  A <- list()
  n <- dim(B)[1]
  for (j in 1:n) {
    # Find nonzero elements in each row
    A <- append(A, list(which(B[j, ] == 1)))
  }
  return(A)
}


#' Permute a sparse matrix
#' @param A List, sparse representation of a matrix
#' @param per Vector, permutation of indices
#' @return Permuted sparse matrix
permute <- function(A, per) {
  B <- list()
  n <- length(A)
  iper <- 1:n
  for (i in 1:n) iper[i] <- which(per == i)  #the inverse permutation so iper[i] is telling from where i came 
  for (j in 1:n) {
    # Apply permutation to the sparse representation
    B <- append(B, list(iper[A[[per[j]]]])) 
  }
  return(B)
}


#' Generate a block tridiagonal matrices matrix
#' @param k Integer,  breadth (block size)
#' @param N Integer, height 
#' @return A block matrix with block tridiagonal matrices structure
bbmat <- function(k, N) {
  A <- bamat(1, 2^(N+1) - 1)  
  A <- sp2mat(A)              
  A[1, nrow(A)] <- A[nrow(A), 1] <- 0  # Remove wrap-around connections
  D <- matrix(1, k, k)        
  A <- A %x% D                # Kronecker product to expand structure
  return(A)
}


#' Generate a permuted band matrix
#' @param d Integer, matrix size
#' @param l Integer, bandwidth
#' @param p Numeric, probability for random sparsification
#' @return List with sparse matrix and permutation vector
f_pbm <- function(d, l, p) {
  A <- bamat(l, d)   
  B <- sp2mat(A)     
  
  C <- diag(d)
  for (i in 2:d) {
    # Introduce random sparsity
    C[1:(i-1), i] <- C[i, 1:(i-1)] <- sample(c(0,1), i-1, replace = TRUE, prob = c(1-p, p))
  }
  
  B <- B * C
  A <- mat2sp(B)  
  
  per <- sample(1:d)
  A <- permute(A, per)
  
  return(list(A = A, per = order(per)))
}

#' Generate a permuted block band tridiagonal matrix 
#' @param k Integer,  breadth (block size)
#' @param N Integer, height 
#' @param l Integer, bandwidth
#' @param p Numeric, probability for sparsification
#' @param Id Logical, whether to use identity permutation
#' @return List with sparse matrix and permutation vector
f_pbbm <- function(k, N, l, p, Id = FALSE) {
  A1 <- bbmat(k, N)  
  d <- nrow(A1)      
  A <- bamat(l, d)   
  B <- sp2mat(A)     
  
  C <- diag(d)
  for (i in 2:d) {
    # Introduce random sparsity
    C[1:(i-1), i] <- C[i, 1:(i-1)] <- sample(c(0,1), i-1, replace = TRUE, prob = c(1-p, p))
  }
  
  B <- B * C
  A1 <- A1 * B  
  A2 <- mat2sp(A1)  
  
  per <- if (Id) 1:d else sample(1:d)
  A2 <- permute(A2, per)
  
  return(list(A = A2, per = order(per)))
}


#' Generate a dyadic matrix
#' @param k Integer, breadth 
#' @param N Integer, height
#' @return Dyadic matrix

dyadmat <- function(k, N) {
  if (N < 1 || !is.numeric(N) || N%%1 != 0 || k < 1 || !is.numeric(k) || k%%1 != 0) {
    stop("N and k must be integers greater than 0")
  }
  M <- 1
  d1 <- 1
  if (N == 1) {
    return(matrix(1, k, k))
  }
  if (N >= 2) {
    for (i in 1:(N - 1)) {
      # Construct dyadic pattern recursively
      A <- matrix(0, 2 * d1 + 1, 2 * d1 + 1)
      A[1:d1, 1:d1] <- A[(d1 + 2):(dim(A)[1]), (d1 + 2):(dim(A)[1])] <- M
      A[d1 + 1, ] <- A[, d1 + 1] <- 1
      M <- A
      d1 <- dim(M)[1]
    }
    D <- matrix(1, k, k)
    A <- A %x% D  # Expand using Kronecker product
  }
  return(A)
}


#' Generate a permuted dyadic band matrix (DBM)
#' @param k Integer, breadth 
#' @param N Integer, height
#' @param l Integer, bandwidth
#' @param p Numeric, probability for sparsification
#' @param Id Logical, whether to use identity permutation
#' @return List with sparse matrix and permutation vector
f_dbm <- function(k, N, l, p, Id = FALSE) {
  A1 <- dyadmat(k, N)  # Generate dyadic matrix
  d <- nrow(A1)        # Get matrix size
  A <- bamat(l , d)     # Generate band matrix in sparse representation
  B <- sp2mat(A)       # Convert sparse representation to matrix form
  
  # Introduce additional zeros in the band matrix
  C <- diag(d)
  for (i in 2:d) {
    C[1:(i-1), i] <- C[i, 1:(i-1)] <- sample(c(0,1), i-1, replace = TRUE, prob = c(1-p, p))
  }
  
  B <- B * C          # Apply random sparsification
  A1 <- A1 * B        # Element-wise multiplication with dyadic matrix
  A2 <- mat2sp(A1)    # Convert back to sparse representation
  
  # Determine whether to use identity permutation or random permutation
  per <- if (Id) 1:d else sample(1:d)
  A2 <- permute(A2, per)
  
  return(list(A = A2, per = order(per)))
}

