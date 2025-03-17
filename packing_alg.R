library(data.table)
library(Rcpp)

# packing_alg.R
# ======================================
# Functions for determining the optimal permutation
# ======================================


power = function(A,s){
  B = A
  for(i in 1:length(A)){
    B[[i]] = neighborhood_s_cpp(i,A,s)[,1]
  }
  return(B)
}




#' Breadth-First Search (BFS) algorithm
#' @param u Integer, the starting vertex
#' @param A List, adjacency list representing the graph
#' @return Integer vector, distances from the starting vertex to all other vertices
BSF_algorithm <- function(u, A) {
  distances <- rep(-length(A)-1, length(A))  # Initialize distances with a large negative value
  distances[u] <- 0  # Distance to the starting vertex is 0
  
  queue <- as.list(u)  # Initialize a queue with the starting vertex u
  
  while (length(queue) > 0) {
    v <- queue[[1]]  # Dequeue the first vertex from the queue
    queue <- queue[-1]  # Remove the dequeued vertex
    
    for (neighbor in A[[v]]) {  # Explore neighbors of v
      if (distances[neighbor] == -length(A)-1) {
        distances[neighbor] <- distances[v] + 1  # Update the distance
        queue <- append(queue, neighbor)  # Enqueue the neighbor
      }
    }
  }
  return(distances)
}


#' Align two sets of coordinates
#' @param a Data table, first set of coordinates
#' @param b Data table, second set of coordinates
#' @return Data table, aligned coordinates
align = function(a,b){
  joined = as.matrix(merge(a,b,by = "vertex"))
  
  if(cov(joined)[2,4]<0){ 
    a[,2] = -1*a[,2]
    joined[,2] = -1*joined[,2]
  }
  
  a[,2] = a[,2] - mean(joined[,2]) + mean(joined[,4])
  return(a)
}

#' Find the final permutation of vertices
#' @param A List, adjacency list representing the graph
#' @param v Integer, starting vertex (optional)
#' @return Integer vector, permutation of vertices
find_per_final = function(A, v=0){
  n_A = length(A)
  
  lengths <- vapply(A, length, integer(1))
  max_length <- max(lengths)
  max_indices <- which(lengths == max_length)
  if (length(max_indices) == 1) {
    v = max_indices
  } else {
    v = sample(max_indices, 1)
  }
  
  graph_dist = data.table(vertex = 1:n_A, g_dist = BSF_algorithm(v,A))
  setkey(graph_dist, g_dist)
  s_max = max(graph_dist[,2])
  
  v_nei_s = sort(A[[v]])
  v_B = dist_mat_rcpp(A[v_nei_s])
  v_res = as.vector(cmdscale(v_B,k=1))
  
  init_ord = data.table(vertex = v_nei_s, 
                        pos = v_res, 
                        Layer = rep(1,length(v_nei_s)) , 
                        key = "vertex")
  setkey(init_ord,vertex)
  
  if(s_max>1){
    for(s in 2:s_max){
      X = graph_dist[g_dist  == s]
      Y = graph_dist[g_dist  == s-1]
      setkey(Y, vertex)
      
      while(nrow(X)>0){
        v = X[sample(1:nrow(X), 1),'vertex']
        
        v_nei_s = A[[as.numeric(v)]]
        Y_sub = Y[Y$vertex %in% v_nei_s]
        setkey(Y_sub, g_dist)
        
        lengths_Y_sub = lengths[as.numeric(Y_sub$vertex)]
        max_length_Y_sub <- max(lengths_Y_sub)
        max_indices <- which(lengths_Y_sub == max_length_Y_sub)
        if (length(max_indices) == 1) {
          w = as.numeric(Y_sub[max_indices,'vertex'])
        } else {
          w = as.numeric(Y_sub[sample(max_indices, 1),'vertex'])
        }
        
        w_nei_s = sort(A[[w]])
        w_B = dist_mat_rcpp(A[w_nei_s])
        w_res = as.vector(cmdscale(w_B,k=1))
        
        v_ord = data.table(vertex = w_nei_s, 
                           pos = w_res, 
                           Layer = rep(s,length(w_nei_s)) , 
                           key = "vertex")
        v_ord = align(v_ord,init_ord)    
        init_ord = rbind(init_ord,v_ord)
        setkey(init_ord,vertex)
        
        X = X[! vertex  %in% as.matrix(w_nei_s)]
        
      }
    }
  }
  
  result = as.matrix(init_ord[,mean(pos),by = vertex])
  result = result[order(result[,2]),]
  
  per = order(result[,1])
  per = order(per)
  return(per)
}





#' Compute a distance matrix using Rcpp
#' @param B List, adjacency list representing the graph
#' @return Numeric matrix, computed distance matrix
Rcpp::cppFunction('
NumericMatrix dist_mat_rcpp(List B) {
  int n = B.size();
  NumericMatrix M(n, n);

  for (int i = 0; i < n - 1; i++) {
    IntegerVector vec_i = B[i];
    for (int j = i + 1; j < n; j++) {
      IntegerVector vec_j = B[j];
      
      int len_i = vec_i.size();
      int len_j = vec_j.size();
      
      std::set<int> intersect_set;
      for (int x : vec_i) {
        if (std::find(vec_j.begin(), vec_j.end(), x) != vec_j.end()) {
          intersect_set.insert(x);
        }
      }
      int intersect_size = intersect_set.size();
      
      double dist = (len_i / 2.0) + (len_j / 2.0) - intersect_size;
      M(i, j) = dist;
      M(j, i) = dist;
    }
  }
  return M;
}
')



Rcpp::cppFunction('
#include <queue>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix neighborhood_s_cpp(int u, const List& A, int s) {
  // Breadth-First Search (BFS) to find the s-step neighborhood of row u
  
  std::queue<std::pair<int, int>> bfs_queue;  // Queue to process row indices
  bfs_queue.push({u, 0});                     // Start with row u 
  
  std::unordered_set<int> visited;  // Set to track visited row indices
  visited.insert(u);
  
  std::vector<std::pair<int, int>> result;  // Store pairs (row index, distance)
  
  while (!bfs_queue.empty()) {
    auto current = bfs_queue.front();  // Get the first row index in the queue
    bfs_queue.pop();
    int v = current.first;   // Current row index
    int depth = current.second;  // Current distance
    
    result.push_back({v, depth});  // Store the row index and its distance
    
    // If limit distance is not reached, explore neighbors
    if (depth < s) {
      IntegerVector neighbors = A[v - 1];  // Adjust index (R is 1-based)
      
      for (int neighbor : neighbors) {
        // Add neighbor to the queue if it has not been visited yet
        if (visited.find(neighbor) == visited.end()) {
          visited.insert(neighbor);
          bfs_queue.push({neighbor, depth + 1});
        }
      }
    }
  }
  
  // Convert result vector into an R IntegerMatrix
  IntegerMatrix result_matrix(result.size(), 2);
  for (size_t i = 0; i < result.size(); ++i) {
    result_matrix(i, 0) = result[i].first;   // row index
    result_matrix(i, 1) = result[i].second;  // Distance from the starting row index
  }
  
  return result_matrix;
}')




