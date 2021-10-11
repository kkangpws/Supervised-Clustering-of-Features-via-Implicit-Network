##########################################################
####Simulation ###########################################
####Very large p##########################################
####SCAD##################################################
##########################################################
setwd('C:\\Users\\ars_w\\Desktop\\Programming Code\\Codes_Submission')

###Library
library(MASS)
library(ncvreg)
library(glmnet)
library(lars)
library(FitAR)
library(itsmr)
library(e1071)
library(ppcor)
library(sarima)
library(xlsx)

source('automatic functions.r')

####Basic settings
set.seed(111) #Set seed
nsim = 5000 #Number of simulations
n=100 #Number of observations
p=9 # Number of parameters
#mu = matrix(c(1,1,1,-1,-1,-1,2,2,2),p,1) # Vector of means
mu = matrix(0,p,1)
sigma = matrix(0,p,p) #Covariance matrix
for (i in 1:p){
  for (j in 1:p){
    if (beta[i]==beta[j]) sigma[i,j]=.5
    else sigma[i,j]=0
  }
}
diag(sigma)=1

#sigma[1:3,1:3] = sigma[4:6,4:6] = sigma[7:9,7:9]=.5
#diag(sigma)=1


beta = matrix(0,p,1)
beta[1:3]=1
beta[4:6]=-1
beta[7:9]=2

c=3 #Number of clusters


deg_cluster = matrix(0,nsim, p) #Clusters based on Euclidean dist
deg_cluster1 = matrix(0,nsim,p) #Clusters based on correlation among X
deg_cluster2 = matrix(0,nsim,p) #Clusters based on partial correlation
kmeans_dist = matrix(0,nsim,p)
kmeans_corr = matrix(0,nsim,p)

for (m in 1:nsim){
  
  
  x = mvrnorm(n,mu=mu, Sigma=sigma)
  y = x%*%beta + rnorm(n)
  
  
  
  #######Degree centrality based on the Euclidean distance (no Y)
  d = matrix(0,p,p)
  
  for (i in 1:p){
    for (j in 1:p) d[i,j] = dist(t(cbind(x[,i],x[,j])), method="euclidean")
  }
    
  degree_dist = apply(d,2,sum)
  
  ##K-means clustering based on Euclidean distance
  kmeans_dist[m,]=kmeans(d,3)$cluster
  
  #####Degree centrality based on the correlation among X (no Y)
  corr = cor(x)
  diag(corr)=0
  degree_corr = apply(corr,2,sum)
  
  ##K-means clustering based on correlation
  kmeans_corr[m,] = kmeans(corr,3)$cluster
  
  
  # 
  # ######Degree centrality based on partial correlation of Y and X
  # degree_par = pcor.degree_new(x, y,seq(1:p))[,2]
  # 
  # 
  # #####Bootstrap to estimate the covariance of degree centrality
  # B = 500
  # degree_dist_boot = matrix(0,B,p)
  # degree_corr_boot = matrix(0,B,p)
  # degree_par_boot = matrix(0,B,p)
  # for (b in 1:B){
  #   s = sample(seq(1:n),n,replace=TRUE)
  #   x_boot = x[s,]
  #   y_boot = y[s]
  #   ###Euclidean distance
  #   d_boot = matrix(0,p,p)
  #   for (i in 1:p){
  #     for (j in 1:p) d_boot[i,j] = dist(t(cbind(x_boot[,i],x_boot[,j])), method="euclidean")
  #   }
  #   deg = apply(d_boot,2,sum)
  #   degree_dist_boot[b,] = deg
  #   
  #   ###Correlation
  #   corr_boot = cor(x_boot)
  #   diag(corr_boot)=0
  #   deg_corr = apply(corr_boot,2,sum)
  #   degree_corr_boot[b,]=deg_corr
  #   
  #   ###Partial correlations
  #   deg_par = pcor.degree_new(x_boot,y_boot,seq(1:p))[2]
  #   for (i in 1:p) degree_par_boot[b,i] = deg_par[i,1]
  # }
  # 
  # deg_dist_cov = cov(degree_dist_boot)
  # deg_corr_cov = cov(degree_corr_boot)
  # deg_par_cov = cov(degree_par_boot)
  # 
  # #source('automatic functions.r')
  # 
  # cluster_dist=sequential_max(n, p,.05, c,degree_dist, deg_dist_cov,seq(1:9),0)
  # for (i in 1:c){
  #   for (j in 1:p) {
  #     if (cluster_dist[i,j+1]!=0) deg_cluster[m,j]=cluster_dist[i,1]
  #   }
  # }
  # 
  # cluster_corr=sequential_min(n, p,.05, c,degree_corr, deg_corr_cov,seq(1:9),0)
  # for (i in 1:c){
  #   for (j in 1:p) {
  #     if (cluster_corr[i,j+1]!=0) deg_cluster1[m,j]=cluster_corr[i,1]
  #   }
  # }
  # 
  # cluster_par=sequential_min(n, p,.05, c,degree_par, deg_par_cov,seq(1:9),0)
  # for (i in 1:c){
  #   for (j in 1:p) {
  #     if (cluster_par[i,j+1]!=0) deg_cluster2[m,j]=cluster_par[i,1]
  #   }
  # }
  # 

    # for (i in 1:c){
    #   cluster_close_multi[(b-1)*c+i,1]=i
    #   for (j in 2:p_sel){
    #     cluster_close_multi[(b-1)*c+i,j] = cluster_close_sequential[i,j]
    #   }
    # }
   
  #if (m==3) break
  print(m)
}


###All possibilities for true clusters
cluster_ver1=c(1,1,1,2,2,2,3,3,3)
cluster_ver2=c(1,1,1,3,3,3,2,2,2)
cluster_ver3=c(2,2,2,1,1,1,3,3,3)
cluster_ver4=c(2,2,2,3,3,3,1,1,1)
cluster_ver5=c(3,3,3,2,2,2,1,1,1)
cluster_ver6=c(3,3,3,1,1,1,2,2,2)

count_k=0
count_k1=0
count1=0
count2=0
count3=0

for (m in 1:nsim){
  if (all(kmeans_dist[m,]==cluster_ver1) || all(kmeans_dist[m,]==cluster_ver2) || all(kmeans_dist[m,]==cluster_ver3) || 
     all(kmeans_dist[m,]==cluster_ver4) || all(kmeans_dist[m,]==cluster_ver5) || all(kmeans_dist[m,]==cluster_ver6))  count_k= count_k+1
}

for (m in 1:nsim){
  if (all(kmeans_corr[m,]==cluster_ver1) || all(kmeans_corr[m,]==cluster_ver2) || all(kmeans_corr[m,]==cluster_ver3) || 
      all(kmeans_corr[m,]==cluster_ver4) || all(kmeans_corr[m,]==cluster_ver5) || all(kmeans_corr[m,]==cluster_ver6))  count_k1= count_k1+1
}

for (m in 1:nsim){
  if (all(deg_cluster[m,]==cluster_ver1) || all(deg_cluster[m,]==cluster_ver2) || all(deg_cluster[m,]==cluster_ver3) || 
      all(deg_cluster[m,]==cluster_ver4) || all(deg_cluster[m,]==cluster_ver5) || all(deg_cluster[m,]==cluster_ver6))  count1= count1+1
}

for (m in 1:nsim){
  if (all(deg_cluster1[m,]==cluster_ver1) || all(deg_cluster1[m,]==cluster_ver2) || all(deg_cluster1[m,]==cluster_ver3) || 
      all(deg_cluster1[m,]==cluster_ver4) || all(deg_cluster1[m,]==cluster_ver5) || all(deg_cluster1[m,]==cluster_ver6))  count2= count2+1
}

for (m in 1:nsim){
  if (all(deg_cluster2[m,]==cluster_ver1) || all(deg_cluster2[m,]==cluster_ver2) || all(deg_cluster2[m,]==cluster_ver3) || 
      all(deg_cluster2[m,]==cluster_ver4) || all(deg_cluster2[m,]==cluster_ver5) || all(deg_cluster2[m,]==cluster_ver6))  count3= count3+1
}

count_k/nsim
count_k1/nsim
count1/nsim
count2/nsim
count3/nsim
