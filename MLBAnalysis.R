#######################################
#######Brandon Woosuk Park#############
######George Mason University##########
#######################################

#######################################
####Data Analysis######################
#######################################

##This file corresponds with teh Table 5 in Secton 7.

setwd('C:\\Users\\ars_w\\Desktop\\Programming Code\\Codes_Submission')


###Load originial data
mlb = read.csv(file="MLB.csv", header=T)
mlb = as.matrix(mlb[,2:36], row.names=mlb[,1])


##Source file
source("automatic functions.r")

##Needed library
library(Matrix)
library(e1071)
library(igraph)
library(MASS)
library(xtable)
library(matrixcalc)
library(ppcor)
library(glmnet)
library(lars)
library(ncvreg)
library(qgraph)
library(tnet)
library(reshape2)


####Create X and y
n = dim(mlb)[1]
col = dim(mlb)[2]
p = col-1

x = data.matrix(mlb[,1:p])
x = scale(x)
colnames(x) = seq(1:p)
y = data.matrix(mlb[,col])
y = scale(y)

###Multiple Splits####
l=30  # l = number of multiple splits
c = 3 # c = number of clusters desired
cluster_degree_multi = matrix(0,l*c, p+1)
cluster_cluster_multi = matrix(0,l*c, p+1)
var_selected = matrix(0,l,p)


####Multiple splits for cleand and screen procedure


#Set seed
set.seed(1111)
for (b in 1:l){
  ###Split the data into two parts
  s = sample(seq(1:n), n, replace=F)
  s1 = s[1:((n+1)/2)]
  s2 = s[(((n+1)/2)+1):n]
  
  x_train = x[s1,]
  x_test = x[s2,]
  y_train = y[s1]
  y_test = y[s2]
  
  # ####Selecting Variables via LASSO-Screen####
  # cv.out = cv.glmnet(x_train,y_train,alpha=1)
  # bestlam = cv.out$lambda.min
  # out = glmnet(x_train,y_train,alpha=1, lambda=bestlam, intercept= FALSE)
  # lasso.coef = predict(out, type="coefficients", s=bestlam)[2:(p+1),]
  # xidx = which(lasso.coef!=0)
  # lp = length(xidx)
  # if (lp==1) break
  # x_train_lasso = x_train[,xidx]
  # beta_lasso = lasso.coef[xidx]
  # x_lasso_name = as.numeric(colnames(x_train_lasso))
  # 
  # # ###Selecting Variables via MCP####
  # cv.out = cv.ncvreg(x_train,y_train,alpha=1, penalty="MCP")
  # bestlam = cv.out$lambda.min
  # options(warn=-1)
  # out = ncvreg(x_train,y_train,alpha=1,penalty="MCP", lambda= bestlam, intercept= FALSE)
  # mcp.coef = out$beta[-1]
  # xidx = which(mcp.coef!=0)
  # lp = length(xidx)
  # if (lp==1) next
  # x_train_lasso = x_train[,xidx]
  # beta_lasso = mcp.coef[xidx]
  # x_lasso_name = as.numeric(colnames(x_train_lasso))
  
  ####Selecting Variables via SCAD####
  cv.out = cv.ncvreg(x_train,y_train,alpha=1, penalty="SCAD")
  bestlam = cv.out$lambda.min
  options(warn=-1)
  out = ncvreg(x_train,y_train,alpha=1,penalty="SCAD", lambda= bestlam, intercept=FALSE)
  scad.coef = out$beta[-1]
  xidx = which(scad.coef!=0)
  lp = length(xidx)
  if (lp==1) next
  x_train_lasso = x_train[,xidx]
  beta_lasso = scad.coef[xidx]
  x_lasso_name = as.numeric(colnames(x_train_lasso))
  
  ####Clean procedure with t test#####
  ##Function for t score
  var_sel = matrix(0,1,lp)
  colnames(var_sel)=x_lasso_name
  for (i in 1:lp){
    reg = summary(lm(y_train~x_train_lasso[,i]))$coefficients[2,4]
    if (reg<.1) var_sel[i]=1
    else var_sel[i]=0
  }
  xxidx = x_lasso_name[which(var_sel!=0)]
  
  p_sel = length(xxidx)
  x_sel = x_test[,xxidx]
  
  for (i in 1:p_sel) var_selected[b,xxidx[i]]=1

  ##Compute degree centrality
  degree = pcor.degree_new(x_sel, y_test,xxidx)[,2]

  
  #Compute clustering coeffiients
  clustercoef = pcor.cluster_new(x_sel, y_test,xxidx)[,2]
  
  
  ##Bootstrap to estimate the variance
  B=500
  degree_boot = matrix(0,B,p_sel)
  cluster_boot = matrix(0,B,p_sel)
  
  colnames(degree_boot) = xxidx
  colnames(cluster_boot) = xxidx
  
  for (m in 1:B){
    s = sample(seq(1:length(s2)), length(s2),replace=T)
    x_boot = x_sel[s,]
    y_boot = y_test[s]
    deg = pcor.degree_new(x_boot,y_boot,xxidx)[,2]
    clc = pcor.cluster_new(x_boot,y_boot,xxidx)[,2]
    for (j in 1:p_sel){
      degree_boot[m,j]=deg[j]
      cluster_boot[m,j]=clc[j]
    }
    #print(m)
  }
  
  ###Compute covariance matrix based on bootstrap estimates
  degree_boot_cov = cov(degree_boot)
  cluster_boot_cov = cov(cluster_boot)
  
  
  ###Sequential test
  multitest = (p_sel*(p_sel-1))/2

  cluster_degree_sequential=sequential_new(length(s2), p_sel,.1/multitest, c,degree, degree_boot_cov, xxidx,0)
  cluster_cluster_sequential=sequential_min(length(s2), p_sel,.1/multitest, c,clustercoef, cluster_boot_cov,xxidx,0)

  for (i in 1:c){
    cluster_degree_multi[(b-1)*c+i,1]=i
    for (j in 2:p_sel){
    cluster_degree_multi[(b-1)*c+i,j] = cluster_degree_sequential[i,j]
    }
  }


  for (i in 1:c){
    cluster_cluster_multi[(b-1)*c+i,1]=i
    for (j in 2:p_sel){
      cluster_cluster_multi[(b-1)*c+i,j] = cluster_cluster_sequential[i,j]
    }
  }

  
  print(b)
}

####Count the number of times each variable is selected
var_sum = apply(var_selected,2,sum)
var_freq = var_sum/l
selected = which(var_freq >=.5)



###For sequential testing
###Count the number of appearances
count_deg = matrix(0,c,p+1)
count_cluster = matrix(0,c,p+1)

for (i in 1:c){
  count_deg[i,1]=i
  count_cluster[i,1]=i
}

for (b in 1:l){
  for (i in 1:c){
    if (cluster_degree_multi[(b-1)*c+i,1]==i){
      for (j in 2:(p+1)){
        if (cluster_degree_multi[(b-1)*c+i,j]!=0) count_deg[i,cluster_degree_multi[(b-1)*c+i,j]+1] = count_deg[i,cluster_degree_multi[(b-1)*c+i,j]+1]+1 
      }
    }
  }
}


for (b in 1:l){
  for (i in 1:c){
    if (cluster_cluster_multi[(b-1)*c+i,1]==i){
      for (j in 2:(p+1)){
        if (cluster_cluster_multi[(b-1)*c+i,j]!=0) count_cluster[i,cluster_cluster_multi[(b-1)*c+i,j]+1] = count_cluster[i,cluster_cluster_multi[(b-1)*c+i,j]+1]+1 
      }
    }
  }
}

deg_freq = count_deg[,-1]/l
cluster_freq = count_cluster[,-1]/l


for (i in 1:c){
  for (j in 1:p){
    if (deg_freq[i,j]>=.5) deg_freq[i,j]=j
    else deg_freq[i,j]=0
  }
}


for (i in 1:c){
  for (j in 1:p){
    if (cluster_freq[i,j]>=.5) cluster_freq[i,j]=j
    else cluster_freq[i,j]=0
  }
}


selected
deg_freq 
cluster_freq

