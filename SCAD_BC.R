##################################
##########Brandon Park############
#####George Mason University######
##################################


######################################################################
####Data Analysis for Scheetz data ###################################
#####17322 variables and 536 observations#############################
######################################################################

###Variable selection method is SCAD + t test

##This file corresponds with the Table 5 in Section 8.


setwd('C:\\Users\\ars_w\\Desktop\\Programming Code\\Codes_Submission')


load("bcTCGA.Rdata")
###Library
library(MASS)
library(ncvreg)
library(Matrix)
library(e1071)
library(xtable)
library(matrixcalc)
library(ppcor)
library(glmnet)
library(SIS)

###Source
source('automatic functions.R')

#####Data is called by opening workspace
#X = scale(X,center=T,scale=T)
#######################################################################
########Comparison among different methodologies begins here###########
#######################################################################

n=dim(X)[1] #number of observations
p = dim(X)[2] #Number of parameters

####Choose 3000 genes with highest variance
x_var= apply(X,2,var)
x_idx = seq(1:p)
x_sel = X[,which(order(x_var)<=3000)]
x_idx_sel = x_idx[which(order(x_var)<=3000)]
x_sel_names = colnames(X)[x_idx_sel]

####Select q genes with largest absolute correlation with Y
q=500
corr = cor(y, x_sel)
xx = x_sel[,which(order(corr)<=q)]
xx_idx = x_idx_sel[which(order(corr)<=q)]
colnames(xx) = xx_idx
xx_names = colnames(X)[xx_idx]
#####Split data into two parts for Screen & Clean procedure and clustering
set.seed(1111)

###Multiple Splits####
l=20  # l = number of multiple splits
c = 3 # c = number of clusters desired
cluster_degree_multi = matrix(0,l*c, p+1)
cluster_cluster_multi = matrix(0,l*c, p+1)
cluster_degree_multi_group = matrix(0,l*c, p+1)
cluster_cluster_multi_group = matrix(0,l*c, p+1)
var_sc = matrix(0,l,q) #Selected variables of l splits
var_sc_scad = matrix(0,l,q) #selected variables from scad of l splits
for (b in 1:l){
  #Split data into two parts
  s = sample(seq(1:n), n, replace=F)
  s1 = s[1:(n/2)]
  s2 = s[((n/2)+1):n]
  
  x_vs = xx[s1,] #Design matirx for a variable selection
  y_vs = y[s1] # Response variable for a variable selection
  x_cl = xx[s2,] #Design matrix for clustering
  y_cl = y[s2] #Response variable for clustering
  
  ##########Clean and Screen##################
  ####Selecting Variables via scad####
  cv.out = cv.ncvreg(x_vs,y_vs,alpha=1, penalty="SCAD", family="gaussian")
  bestlam = cv.out$lambda.min
  options(warn=-1)
  out = ncvreg(x_vs,y_vs,alpha=1, family="gaussian",penalty="SCAD", lambda= bestlam, intercept= FALSE)
  mcp.coef = out$beta[-1]
  xidx = which(mcp.coef!=0)
  lp = length(xidx)
  if (lp==1) next
  x_vs_mcp = x_vs[,xidx]
  colnames(x_vs_mcp)=xidx
  beta_mcp = mcp.coef[xidx]
  x_mcp_name = colnames(x_vs_mcp)
  
  #for (i in 1:lp) var_mcp_selected[b,xidx[i]]=1
  ####Clean procedure with t test#####
  ##Function for t score
  
  var_sel = matrix(0,1,lp)
  pval = matrix(0,1,lp)
  colnames(var_sel)=x_mcp_name
  for (i in 1:lp){
    reg = summary(lm(y_vs~x_vs_mcp[,i]))$coefficients[2,4]
    pval[i]=reg
    if (reg<.1) var_sel[i]=1
    else var_sel[i]=0
  }
  xxidx = as.numeric(x_mcp_name[which(var_sel!=0)])
  
  p_sel = length(xxidx)
  #x_sel = x_test[,xxidx]
  #colnames(x_sel)=xxidx
  for (i in 1:p_sel) var_sc[b,xxidx[i]]=1
  

  
  if (p_sel<=2) next
  #print(b)
  

  x_sel = x_cl[,xxidx]
  
  #Compute degree centrality
  degree = pcor.degree_new(x_sel, y_cl,xxidx)[,2]
  
  #for (i in 1:p_sel) degree_cent[m,var_sc_sel[i]]=degree[i]
  
  
  #Compute clustering coefficients
  clustercoef = pcor.cluster_new(x_sel, y_cl,xxidx)[,2]
  
  #for (i in 1:p_sel) clustering_coef[m,var_sc_sel[i]]=clustercoef[i]
  
  
  ##Bootstrap to estimate the variance
  B=500
  degree_boot = matrix(0,B,p_sel)
  cluster_boot = matrix(0,B,p_sel)
  
  for (k in 1:B){
    s_boot = sample(seq(1:length(s2)), length(s2),replace=T)
    x_boot = x_sel[s_boot,]
    y_boot = y_cl[s_boot]
    deg = pcor.degree_new(x_boot,y_boot,xxidx)[,2]
    clc = pcor.cluster_new(x_boot,y_boot,xxidx)[,2]
    for (j in 1:p_sel){
      degree_boot[k,j]=deg[j]
      cluster_boot[k,j]=clc[j]
    }
    #print(m)
  }
  degree_boot_cov = cov(degree_boot)
  cluster_boot_cov = cov(cluster_boot)
  
  ###Sequential test
  multitest = (p_sel*(p_sel-1))/2
  #Source file
  #source('automatic functions.r')
  #c = 3 # c = number of clusters desired
  cluster_degree_sequential=sequential_min(length(s2), p_sel,.05, c,degree, degree_boot_cov,xxidx,0)
  cluster_cluster_sequential=sequential_new(length(s2), p_sel,.05, c,clustercoef, cluster_boot_cov,xxidx,0)
  
  
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
var_sum = apply(var_sc,2,sum)
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


###Assign a variable to a cluster if it appears more than or equal to 50% in a cluster
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



d_len = matrix(0,1,c)
cc_len = matrix(0,1,c)

for (i in 1:c){
  d_len[i]=length(which(deg_freq[i,]!=0))
  cc_len[i]=length(which(cluster_freq[i,]!=0))
}


#store finalized clusters
deg_seq_final = matrix(0,c,max(d_len))
cluster_seq_final = matrix(0,c,max(cc_len))
for (i in 1:c){
  for (j in 1:d_len[i]) deg_seq_final[i,j] = which(deg_freq[i,]!=0)[j]
}

for (i in 1:c){
  for (j in 1:cc_len[i]) cluster_seq_final[i,j] = which(cluster_freq[i,]!=0)[j]
}

xx_names[selected]
deg_seq_final
cluster_seq_final 

