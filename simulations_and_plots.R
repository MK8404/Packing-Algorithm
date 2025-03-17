require(ggplot2)
source("matrix_generation.r")
source("packing_alg.r")
rotate <- function(x) t(apply(x, 2, rev))

#========================================================================
#========================================================================
#Fig 1. The consequences of increasing the neighborhood size 
# for finding the optimal permutation:
#========================================================================
#========================================================================
#Permuted full band matrix
A = f_pbm(d = 100,l = 20,p = 1) 

#The first row in the Fig. 1. Plots of band matrices \Sigma(s), s = 1,2,3,4
image(rotate(sp2mat(permute(A$A, A$per))),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

image(rotate(sp2mat(permute(power(A$A,2), A$per))),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

image(rotate(sp2mat(permute(power(A$A,3), A$per))),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

image(rotate(sp2mat(permute(power(A$A,4), A$per))),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

#The second row in Fig.1. The results of permutations obtained from the packing algorithm 
#applied to randomly permuted matrices in the first row.

per = find_per_final(A = A$A,0)
C = permute(A$A,per)
image(rotate(sp2mat(C)),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

per = find_per_final(A = power(A$A,2),0)
C = permute(A$A,per)
image(rotate(sp2mat(C)),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

per = find_per_final(A = power(A$A,3),0)
C = permute(A$A,per)
image(rotate(sp2mat(C)),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

per = find_per_final(A = power(A$A,4),0)
C = permute(A$A,per)
image(rotate(sp2mat(C)),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)



#========================================================================
#========================================================================
# Simulations for permuted band matrices
#========================================================================
#========================================================================
n_rep =  100 # number of repetitions
d = 1000  # number of rows (columns) in the matrix
l_values = c(10,40) # half bandwidth of a full band matrix
p_values = c(0.25,0.5,0.75) # sparsity parameter
s_values = c(1,2,3) # neighborhood order
results1 = data.frame(d = 0,l = 0,p = 0, s = 0, stat = 0, m_half_width = 0,half_bw = 0)
results2 = data.frame(d = 0,l = 0,p = 0, s = 0, F_s = 0)


for(l in l_values){
  A_ref = f_pbm(d = d,l = l,p = 1) #reference full band matrix
  for(p in p_values){
    for(i in 1:n_rep){
      A     = f_pbm(d = d,l = l,p = p) # initial band matrix
      #If the initial matrix is block-diagonal 
      #with at least two blocks, we regenerate the matrix,
      #until we get matrix with one block:
      while(length(A$A) != sum(BSF_algorithm(1,A$A)>=0)){  
        A     = f_pbm(d = d,l = l, p = p)
      }
      
      for(s in s_values){
        AA = power(A$A,s) #AA is a matrix representing row connections in s-th order neighborhoods 
        per = find_per_final(AA,0) # permutation returned by the packing algorithm
        
        # We permute the rows and columns of the matrix A$A 
        #according to the permutation returned by the algorithm:
        C = permute(A$A,per) 
        #the half-widths of D_i's under 'per':
        l_seq = c()
        for(j in 1:length(C)){
          l_seq[j] = max(abs(C[[j]]-j))
        }
        results1 = rbind(results1, c(d,l,p,s,1,mean(l_seq), max(l_seq)))
        
        # We permute the rows and columns of the matrix A$A 
        #according to the initial permutation A$per:
        C = permute(A$A,A$per)
        #the half-widths of D_i's under 'A$per' (original permutation)
        l_seq = c()
        for(j in 1:length(C)){
          l_seq[j] = max(abs(C[[j]]-j))
        }
        results1 = rbind(results1, c(d,l,p,s,2,mean(l_seq), max(l_seq)))
        
        #the level of (band) filling F(s):
        F_s = sum(sp2mat(AA) !=0)/(d*(2*s*l+1)-l*s*(l*s+1))
        results2 = rbind(results2, c(d,l,p,s,F_s))
        
      }
    }
  }
  
}
results1 = results1[-1,] #removing the initial (zero-filled) row 
#Determination of the mean across repetitions for different parameter value combinations:
summar1 = aggregate(results1, by = list(results1$l,results1$p,results1$s,results1$stat), mean)

#Fig.2. Results for the average half-width (\Vert \boldsymbol l^{(\pi)} \Vert_1/d) 
qplot(p, m_half_width, data = summar1, colour = as.factor(stat), geom = c("point", "line")) + 
  facet_wrap(s ~ l, scales = "free_y", ncol = 2, 
             labeller = labeller(s = function(x) paste("s =", x), 
                                 l = function(x) paste("\u03BB =", x), .multi_line = FALSE)) + 
  labs(y = "") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, hjust = 0.5, vjust = 0.5, angle = 0),  # Larger font size for panel labels
    axis.text = element_text(size = 12),  # Larger font size for axis text
    axis.title = element_text(size = 14),  # Larger font size for axis titles
    plot.title = element_text(size = 16, face = "bold")  # Larger font size for the plot title
  ) +
  geom_hline(aes(yintercept = l), linetype = "dashed", color = "black") +
  geom_line(linewidth = 1.0) +  # Thicker lines
  geom_point(size = 2)  # Larger points


#Fig.2. Results for the half-bandwidth:
qplot(p, half_bw, data = summar1, colour = as.factor(stat), geom = c("point", "line")) + 
  facet_wrap(s ~ l, scales = "free_y", ncol = 2, 
             labeller = labeller(s = function(x) paste("s =", x), 
                                 l = function(x) paste("\u03BB =", x), .multi_line = FALSE)) + 
  labs(y = "")  +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, hjust = 0.5, vjust = 0.5, angle = 0),  # Larger font size for panel labels
    axis.text = element_text(size = 12),  # Larger font size for axis text
    axis.title = element_text(size = 14),  # Larger font size for axis titles
    plot.title = element_text(size = 16, face = "bold")  # Larger font size for the plot title
  ) +
  geom_line(size = 1.0) +  # Thicker lines
  geom_point(size = 2)  # Larger points



results2 = results2[-1,] #removing the initial (zero-filled) row 

#Determination of the mean across repetitions for different parameter value combinations:
summar2 = aggregate(results2, by = list(results2$l,results2$p,results2$s), mean)
summar2$s = factor(summar2$s)


#Fig. 3 (left panel). Average results for the level of (band) filling F(s):
qplot(p, F_s, data = summar2[summar2$l == 10,],  
      linetype = s, color = s, geom = c("point", "line")) +  
  labs(y = "") +  
  scale_y_continuous(limits = c(NA, 1)) +  
  theme(  
    legend.position = "none",  
    strip.text = element_text(size = 14, hjust = 0.5, vjust = 0.5, angle = 0),  # Larger font size for panel labels  
    axis.text = element_text(size = 18),  # Larger font size for axis text  
    axis.title = element_text(size = 18),  # Larger font size for axis titles  
    plot.title = element_text(size = 16, face = "bold")  # Larger font size for the plot title  
  ) +  
  geom_line(size = 1.0) +  # Thicker lines  
  geom_point(size = 2)  # Larger points  





#========================================================================
#========================================================================
# Simulations for block tridiagonal matrices
#========================================================================
#========================================================================

n_rep =  100 # number of repetitions
N_values = 4 # height 
k_values = c(10,40) # breadth (block size)
p_values = c(0.25,0.5,0.75,1) # sparsity parameter
s_values = c(1,2,3)  # neighborhood order
results1 = data.frame(N = 0,k = 0,p = 0, s = 0, stat = 0 , m_half_width = 0)
results2 = data.frame(N = 0,k = 0,p = 0, s = 0, G_s = 0)


for(N in N_values){
  for(k in k_values){
    #reference block tridiagonal matrix,
    #the value of the l parameter is chosen in such a way 
    #as not to zero the elements in the block tridiagonal structure    
    A_ref = f_pbbm(k = k,N = N,l = 2*k,p = 1,Id = T)
    for(p in p_values){
      for(i in 1:n_rep){
        A     = f_pbbm(k = k,N = N,l = 2*k, p = p, Id = F)# initial block tridiagonal matrix
        #If the initial matrix is block-diagonal 
        #with at least two blocks, we regenerate the matrix,
        #until we get matrix with one block:
        while(length(A$A) != sum(BSF_algorithm(1,A$A)>=0)){
          A     = f_pbbm(k = k,N = N,l = 2*k, p = p, Id = F)
        }
        for(s in s_values){
          
          AA = power(A$A,s) #AA is a matrix representing row connections in s-th order neighborhoods 
          per = find_per_final(AA,0) # permutation returned by the packing algorithm
          
          # We permute the rows and columns of the matrix A$A 
          #according to the original permutation 'A$per':
          C = permute(A$A, A$per)
          #the half-widths of D_i's under 'A$per' (original permutation)
          l_seq = c()
          for(j in 1:length(C)){
            l_seq[j] = max(abs(C[[j]]-j))
          }
          results1 = rbind(results1, c(N,k,p,s,2, sum(l_seq)/(k*2^N-1)))

          # We permute the rows and columns of the matrix A$A 
          #according to the permutation returned by the algorithm:
          C = permute(A$A,per) 
          #the half-widths of D_i's under 'per':
          l_seq = c()
          for(j in 1:length(C)){
            l_seq[j] = max(abs(C[[j]]-j))
          }
          results1 = rbind(results1, c(N,k,p,s,1, sum(l_seq)/(k*2^N-1)))
          
          # The proportion of nonzero elements in the matrix, reordered by the packing algorithm,
          #within the tridiagonal block structure of the reference matrix:
          B = 1 - sum(sp2mat(C)>sp2mat(A_ref$A))/sum(sp2mat(C))
          results2 = rbind(results2, c(N,k,p,s,B))
          
        }
      }
    }
  }
}


results1 = results1[-1,] #removing the initial (zero-filled) row

#Determination of the mean across repetitions for different parameter value combinations:
summar1 = aggregate(results1, by = list(results1$N,results1$k,results1$p,results1$s,results1$stat), mean)

#Fig.4. Results for the average half-width (\Vert \boldsymbol l^{(\pi)} \Vert_1/d) 
qplot(p, m_half_width, data = summar1[summar1$p < 1,], colour = as.factor(stat), geom = c("point", "line")) +  
  facet_wrap(s ~ k, ncol = 2, scales = "free_y",  
             labeller = labeller(s = function(x) paste("s =", x),  
                                 k = function(x) paste("k =", x), .multi_line = FALSE)) +   
  labs(y = "") +  
  theme(  
    legend.position = "none",  
    strip.text = element_text(size = 14, hjust = 0.5, vjust = 0.5, angle = 0),  # Larger font size for panel labels  
    axis.text = element_text(size = 12),  # Larger font size for axis text  
    axis.title = element_text(size = 14),  # Larger font size for axis titles  
    plot.title = element_text(size = 16, face = "bold")  # Larger font size for the plot title  
  ) +  
  geom_line(size = 1.0) +  # Thicker lines  
  geom_point(size = 2) +  # Larger points  
  geom_hline(data = unique(summar1[, c("s", "k")]),  
             aes(yintercept = summar1[summar1$p == 1 & summar1$stat == 1, 'm_half_width']),  
             linetype = "dashed", color = "black")  # Dashed reference line  

  
  
results2 = results2[-1,] #removing the initial (zero-filled) row
#Determination of the mean across repetitions for different parameter value combinations:
summar2 = aggregate(results2, by = list(results2$N,results2$k,results2$p,results2$s), mean)
summar2$k = factor(summar2$k)
summar2$s = factor(summar2$s)


#Fig 3.(middle and right panel). The average proportion of nonzero elements in the matrices,
# reordered by the packing algorithm, within the tridiagonal block structures of reference matrices:
qplot(p, G_s, data = summar2[summar2$p < 1,], color = s, linetype = s, geom = c("point", "line")) +  
  facet_wrap(~k, ncol = 2, scales = "free_y",  
             labeller = labeller(k = function(x) paste("k =", x), .multi_line = FALSE)) +  
  labs(y = "") +  
  theme(  
    legend.position = "right",  
    legend.text = element_text(size = 16),   # Larger font size for legend text  
    legend.title = element_text(size = 18, face = "bold"),  # Larger font size for legend title  
    strip.text = element_text(size = 18, hjust = 0.5, vjust = 0.5, angle = 0),  # Larger font size for panel labels  
    axis.text = element_text(size = 18),  # Larger font size for axis text  
    axis.title = element_text(size = 18),  # Larger font size for axis titles  
    plot.title = element_text(size = 18, face = "bold")  # Larger font size for the plot title  
  ) +  
  geom_line(size = 1.0) +  # Thicker lines  
  geom_point(size = 2) +  # Larger points  
  scale_y_continuous(limits = c(NA, NA))  # Flexible y-axis limits  


#===================================================================
#===================================================================
#Fig.5. Banded dyadic matrix and dyadic matrix
#===================================================================
#===================================================================
#Banded dyadic matrix
A = f_dbm(k = 10,N = 5,l = 60, p = 0.5, Id = F)

# The first plot in the first row of Fig.5, original non-permuted banded dyadic matrix:
image(rotate(sp2mat(permute(      A$A   , A$per))),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

# The second plot in the first row of Fig.5, second power of original non-permuted banded dyadic matrix:
image(rotate(sp2mat(permute(power(A$A,2), A$per))),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)


# We permute the rows and columns of the matrix A$A 
#according to the permutation returned by the algorithm applied to neighborhoods of order 1:
per = find_per_final(A$A,0)
C = permute(A$A,per)
# The first plot in the second row of Fig.5, matrix recovered by the packing algorithm 
image(rotate(sp2mat(C)),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

# We permute the rows and columns of the matrix A$A 
#according to the permutation returned by the algorithm applied to neighborhoods of order 2:
per = find_per_final(power(A$A,2),0)
C = permute(A$A,per)
# The second plot in the second row of Fig.5, matrix recovered by the packing algorithm 
image(rotate(sp2mat(C)),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)




#Dyadic matrix
A = f_dbm(k = 10,N = 5,l = 160, p = 0.45, Id = F)

# The third plot in the first row of Fig.5, original non-permuted dyadic matrix:
image(rotate(sp2mat(permute(      A$A   , A$per))),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

# The fourth plot in the first row of Fig.5, second power of original non-permuted dyadic matrix:
image(rotate(sp2mat(permute(power(A$A,2), A$per))),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

# We permute the rows and columns of the matrix A$A 
#according to the permutation returned by the algorithm applied to neighborhoods of order 1:
per = find_per_final(A$A,0)
C = permute(A$A,per)
# The third plot in the second row of Fig.5, matrix recovered by the packing algorithm 
image(rotate(sp2mat(C)),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)
# We permute the rows and columns of the matrix A$A 
#according to the permutation returned by the algorithm applied to neighborhoods of order 2:
per = find_per_final(power(A$A,2),0)
C = permute(A$A,per)
# The fourth plot in the second row of Fig.5, matrix recovered by the packing algorithm 
image(rotate(sp2mat(C)),axes = FALSE, xlab = "", ylab = "", col = c("white", "black"))
box(col = "black", lwd = 2)

