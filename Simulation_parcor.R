##########################################################
####Simulation ###########################################
####Degree Centrality#####################################
####Using partial correlations############################
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
nsim = 100 #Number of simulations
n=100 #Number of observations
p=9 # Number of parameters
beta = c(1,-1,2,1,1,2,-1,2,-1) #regresson coefficients


mu = matrix(0,p,1) #Mean vector for X
sigma = matrix(0,p,p) #Covariance matrix of X
for (i in 1:p){
  for (j in 1:p){
    if (beta[i]==beta[j]) sigma[i,j]=.5 #Within-cluster correlation
    else sigma[i,j]=.5 #Between-cluster correlation
  }
}
diag(sigma)=1

c=3 #Number of clusters

deg_cluster = matrix(0,nsim,p) #Clusters based on partial correlation


for (m in 1:nsim){
  
  #Generate X and Y
  x = mvrnorm(n,mu=mu, Sigma=sigma)
  y = x%*%beta + rnorm(n)
  
  
  ######Degree centrality based on partial correlation of Y and X
  degree_par = pcor.degree_new(x, y,seq(1:p))[,2]
  
  
  #####Bootstrap to estimate the covariance of degree centrality
  B = 500
  degree_par_boot = matrix(0,B,p)
  for (b in 1:B){
    s = sample(seq(1:n),n,replace=TRUE)
    x_boot = x[s,]
    y_boot = y[s]

    ###Partial correlations
    deg_par = pcor.degree_new(x_boot,y_boot,seq(1:p))[2]
    for (i in 1:p) degree_par_boot[b,i] = deg_par[i,1]
  }

  deg_par_cov = cov(degree_par_boot)
  
  #source('automatic functions.r')
  
  cluster_par=sequential_min(n, p,.05, c,degree_par, deg_par_cov,seq(1:9),0)
  for (i in 1:c){
    for (j in 1:p) {
      if (cluster_par[i,j+1]!=0) deg_cluster[m,j]=cluster_par[i,1]
    }
  }
  

  
  #if (m==3) break
  print(m)
}


###All possibilities for true clusters
cluster_ver1=c(1,2,3,1,1,3,2,3,2)
cluster_ver2=c(1,3,2,1,1,2,3,2,3)
cluster_ver3=c(2,1,3,2,2,3,1,3,1)
cluster_ver4=c(2,3,1,2,2,1,3,1,3)
cluster_ver5=c(3,2,1,3,3,1,2,1,2)
cluster_ver6=c(3,1,2,3,3,2,1,2,1)

count=0


for (m in 1:nsim){
  if (all(deg_cluster[m,]==cluster_ver1) || all(deg_cluster[m,]==cluster_ver2) || all(deg_cluster[m,]==cluster_ver3) || 
      all(deg_cluster[m,]==cluster_ver4) || all(deg_cluster[m,]==cluster_ver5) || all(deg_cluster[m,]==cluster_ver6))  count= count+1
}

count/nsim
