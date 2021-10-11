##############################################
######Brandon Woosuk Park#####################
######George Mason University#################
##############################################

#This file includes R codes for numerical studies for spectral clustering


setwd('C:\\Users\\ars_w\\Desktop\\Programming Code\\Codes_Submission')

##Source code
source('automatic functions.r')

library(MASS)
library(kernlab)
library(anocva)
library(matrixcalc)

set.seed(111)

p=9 #Number of covariates (X) , Increase the number of covariates
n=100 #Number of samples
k=3 #Number of clusters
nsim = 5000 #Number of simulations



beta = c(1,-1,2,1,1,2,-1,2,-1) #regresson coefficients

mu = matrix(0,p,1) #Mean vector for X
sigma = matrix(0,p,p) #Covariance matrix of X
for (i in 1:p){
  for (j in 1:p){
    if (beta[i]==beta[j]) sigma[i,j]=.5
    else sigma[i,j]=.5
  }
}
diag(sigma)=1

spec_cluster = matrix(0,nsim,p) #matrix storing simulation results


for (m in 1:nsim){
  x = mvrnorm(n,mu,sigma) #Generate X from multivariate normal
  e = rnorm(n) #error follows a standard normal distribution
  
  
  y = x%*%beta+e #Generate Y
  

  ###Spectral Clustering based on R-squared
  
  # r2_total = summary(lm(y~x-1))$r.squared
  # sim = matrix(0,p,p)
  # for (i in 1:p){
  #   for (j in 1:p){
  #     if (i!=j){
  #       r2_red = summary(lm(y~x[,-c(i,j)]-1))$r.squared
  #       sim[i,j] = r2_total-r2_red
  #     }
  #     else sim[i,j]=0
  #   }
  # }
  # 
  # summary(lm(y~x[,-c(3,6)]-1))$r.squared
  ###differences as similarity between two covariates
  sim = sim_r2(x,y) # Similarity matrix
  D = diag(apply(sim,1,sum)) #Degree matrix
  L = D-sim #Laplacian
  #L = solve(D)%*%L #Normalized Laplacian
  ei=eigen(L, symmetric="True") #Compute eigenvalues and corresponding eigen vectors
  Z1 = ei$vectors[,(p-k+1):p] #the k eigenvectors corresponding with k smallest eigenvalues
  spec_cluster[m,]=kmeans(Z1, centers=k, nstart=25)$cluster #Use K-means clustering

  print(m)
}

spec_cluster[is.na(spec_cluster)]=0

###All possibilities for true clusters
# cluster_ver1=c(1,3,1,1,1,1,2,1,3)
# cluster_ver2=c(2,3,2,2,2,2,1,2,3)
# cluster_ver3=c(2,1,2,2,2,2,3,2,1)
# cluster_ver4=c(3,1,3,3,3,3,2,3,1)
# cluster_ver5=c(3,2,3,3,3,3,1,3,2)
# cluster_ver6=c(1,2,1,1,1,1,3,1,2)

cluster_ver1=c(1,1,1,2,3,1,1,1,1)
cluster_ver2=c(1,1,1,3,2,1,1,1,1)
cluster_ver3=c(2,2,2,3,1,2,2,2,2)
cluster_ver4=c(2,2,2,1,3,2,2,2,2)
cluster_ver5=c(3,3,3,1,2,3,3,3,3)
cluster_ver6=c(3,3,3,2,1,3,3,3,3)

count1=0

for (m in 1:nsim){
  if (all(spec_cluster[m,]==cluster_ver1) || all(spec_cluster[m,]==cluster_ver2) || all(spec_cluster[m,]==cluster_ver3) || 
      all(spec_cluster[m,]==cluster_ver4) || all(spec_cluster[m,]==cluster_ver5) || all(spec_cluster[m,]==cluster_ver6))  count1= count1+1
}

count1/nsim



