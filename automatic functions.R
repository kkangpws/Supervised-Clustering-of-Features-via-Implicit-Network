##############################################
######Brandon Woosuk Park#####################
######George Mason University#################
##############################################


#This file includes all R functions that were used for numerical studies and data analyses



##Find mode in a vector
Mode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}




##Compute the norm of a vector
norm <- function(x)
{
  x_sq = sum(x^2)
  return(sqrt(x_sq))
}




##################################
####Weighted Degree Centrality####
##################################

#Following two R functions compute degree centrality of a weight matrix. 
#A weighted matrix is constructed based on the partial correlations
#between the covariate and the response variable

pcor.degree <- function(x,y,b)
{

p = dim(x)[2]
colname = b

######Table for the correlation#######
table.pcor  = matrix(0,nrow=1,ncol=p)
for (i in 1:p){
	table.pcor[i] = pcor.test(x[,i],y,x[,-i])$estimate}

####Make the partial correlation be in between 0 and 1
table = matrix(0,nrow=1,ncol=p)
for (i in 1:p) table[i] = (1+table.pcor[i])/2

#####Adjacency Matrix with Cor Cent####
A = matrix(0,nrow=p,ncol=p)
colnames(A) = colname
rownames(A) = colname
for (j in 1:p){
	for (u in 1:p){
	if (j != u){
		A[j,u]=A[j,u]+sqrt((table[1,j])^2+(table[1,u])^2)} #Weight is the square root of summation of two squared partial correlations
		else{A[j,u]=A[j,u]}}}

num = matrix(colname,p,1)
d=apply(A,1,sum)
return(as.data.frame(list(Node =num, Weighted_Degree = d)))
}



#A different function of partial correlation is used to construct the weight network
pcor.degree_new <- function(x,y,b)
{
  
  p = dim(x)[2]
  colname = b

  ######Table for the correlation#######
  table.pcor  = matrix(0,nrow=1,ncol=p)
  for (i in 1:p){
    table.pcor[i] = pcor.test(x[,i],y,x[,-i])$estimate}
  
  ####Make the partial correlation be in between 0 and 1
  table = matrix(0,nrow=1,ncol=p)
  for (i in 1:p) table[i] = (1+table.pcor[i])/2
  
  for (i in 1:p) table[i] = sqrt(2*(1-table[i]^2))
  
  #####Adjacency Matrix with Cor Cent####
  A = matrix(0,nrow=p,ncol=p)
  colnames(A) = colname
  rownames(A) = colname
  for (j in 1:p){
    for (u in 1:p){
      if (j != u){
        A[j,u]=A[j,u]+table[1,j]+table[1,u]
        }
      else{
        A[j,u]=A[j,u]
      }
    }
  }
  
  num = matrix(colname,p,1)
  d=apply(A,1,sum)
  return(as.data.frame(list(Node =num, Weighted_Degree = d)))
}


#Following codes describe how to compute degree centrality
#of a weighted matrix that is constructed based on the difference
#of R^2 between the full and reduced models
rsqr.degree <- function(x,y,b)
{
  
  p = ncol(x)
  colname = b
  
  ###Construct Similarity Matrix
  A = sim_r2(x,y)
  num = matrix(colname,p,1)
  
  #Compute degree centrality
  d=apply(A,1,sum)
  return(as.data.frame(list(Node =num, Weighted_Degree = d)))
}

###Similarity matrix using r-squared#####
sim_r2 <- function(x,y){
  p = ncol(x)
  sim = matrix(0,p,p)
  r2_total = summary(lm(y~x-1))$r.squared
  for (i in 1:p){
    for (j in 1:p){
      if (i!=j){
        r2_red = summary(lm(y~x[,-c(i,j)]-1))$r.squared
        sim[i,j] = r2_total-r2_red
      }
      else sim[i,j]=0
    }
  }
  return(sim)
}		



###########################################
###Weighted Closeness Centrality###########
###########################################
#Following two R functions compute closeness centrality of a weight matrix. 
#A weighted matrix is constructed based on the partial correlations
#between the covariate and the response variable

pcor.closeness <- function(x,y,b)
{
n = dim(x)[1]
p = dim(x)[2]
data = cbind(y,x)
ss = sample(seq(1:n), n/2, replace=F)
data.train = data[ss,]
data.test = data[-ss,]
y.train = data.train[,1]
x.train = data.train[,2:(p+1)]
y.test = data.test[,1]
x.test = data.test[,2:(p+1)]
colname = b

######Table for the correlation#######
table  = matrix(0,nrow=1,ncol=p)
for (i in 1:p){
	table[i] = pcor.test(x.train[,i],y.train,x.train[,-i])$estimate}

####Make the partial correlation be in between 0 and 1
for (i in 1:p)	table[i] = (1+table[i])/2	
	
#####Adjacency Matrix with Cor Cent####
A = matrix(0,nrow=p,ncol=p)
colnames(A) = colname
rownames(A) = colname
for (j in 1:p){
	for (u in 1:p){
	if (j != u){
		A[j,u]=A[j,u]+sqrt((table[1,j])^2+(table[1,u])^2)}
		else{A[j,u]=A[j,u]}}}

z=allShortestPaths(1/A)$middlepoints


###With the testing data###
######Table for the correlation#######
table.test  = matrix(0,nrow=1,ncol=p)
for (i in 1:p){
	table.test[i] = pcor.test(x.test[,i],y.test,x.test[,-i])$estimate}

for (i in 1:p) table.test[i] = (1+table.test[i])/2
#####Adjacency Matrix with Cor Cent####
A.test = matrix(0,nrow=p,ncol=p)

for (j in 1:p){
	for (u in 1:p){
	if (j != u){
		A.test[j,u]=A.test[j,u]+sqrt((table.test[1,j])^2+(table.test[1,u])^2)}
		else{A.test[j,u]=A.test[j,u]}}}

A.dist = 1/A.test
A.dist[A.dist==Inf]=0
dist = matrix(0,p,p)
if (is.null(z)){
	dist = A.dist
}
else{
for (i in 1:p){
	for (j in 1:p){
		d = extractPath(z,i,j)
		dl = length(d)
		for (k in 1:(dl-1)){
				dist[i,j] = dist[i,j]+ A.dist[d[k],d[k+1]]
		}
	}
}
}
num = matrix(colname,p,1)
cc=1/apply(dist,1,sum)
return(as.data.frame(list(Node =num, Weighted_Closeness = cc)))
}


#A different function of partial correlation is used to construct the weight network
pcor.closeness_new <- function(x,y,b)
{
  n = dim(x)[1]
  p = dim(x)[2]
  data = cbind(y,x)
  ss = sample(seq(1:n), n/2, replace=F)
  data.train = data[ss,]
  data.test = data[-ss,]
  y.train = data.train[,1]
  x.train = data.train[,2:(p+1)]
  y.test = data.test[,1]
  x.test = data.test[,2:(p+1)]
  colname = b

  ######Table for the correlation#######
  table  = matrix(0,nrow=1,ncol=p)
  for (i in 1:p){
    table[i] = pcor.test(x.train[,i],y.train,x.train[,-i])$estimate}
  
  ####Make the partial correlation be in between 0 and 1
  for (i in 1:p)	table[i] = (1+table[i])/2	
  
  for (i in 1:p) table[i] = sqrt(2*(1-table[i]^2))
  
  #####Adjacency Matrix with Cor Cent####
  A = matrix(0,nrow=p,ncol=p)
  colnames(A) = colname
  rownames(A) = colname
  for (j in 1:p){
    for (u in 1:p){
      if (j != u){
        A[j,u]=A[j,u]+table[1,j]+table[1,u]
      }
      else{
        A[j,u]=A[j,u]
      }
    }
  }
  
  z=allShortestPaths(1/A)$middlepoints
  
  
  ###With the testing data###
  ######Table for the correlation#######
  table.test  = matrix(0,nrow=1,ncol=p)
  for (i in 1:p){
    table.test[i] = pcor.test(x.test[,i],y.test,x.test[,-i])$estimate}
  
  for (i in 1:p) table.test[i]=(1+table.test[i])/2
  
  for (i in 1:p) table.test[i] = sqrt(2*(1-table.test[i]^2))

  
  #####Adjacency Matrix with Cor Cent####
  A.test = matrix(0,nrow=p,ncol=p)
  
  for (j in 1:p){
    for (u in 1:p){
      if (j != u){
        A.test[j,u]=A.test[j,u]+table.test[1,j]+table.test[1,u]}
      else{A.test[j,u]=A.test[j,u]}}}
  

  dist = matrix(0,p,p)
  if (is.null(z)){
    dist = A.test
  }
  else{
    for (i in 1:p){
      for (j in 1:p){
        d = extractPath(z,i,j)
        dl = length(d)
        for (k in 1:(dl-1)){
          dist[i,j] = dist[i,j]+ A.test[d[k],d[k+1]]
        }
      }
    }
  }
  num = matrix(colname,p,1)
  cc=1/apply(dist,1,sum)
  return(as.data.frame(list(Node =num, Weighted_Closeness = cc)))
}


################################
###Weighted Clustering Coef#####
################################
#Following two R functions compute clustering coefficient of a weight matrix. 
#A weighted matrix is constructed based on the partial correlations
#between the covariate and the response variable

pcor.cluster <- function(x,y,b)
{

p = dim(x)[2]
colname = b

######Table for the correlation#######
table  = matrix(0,nrow=1,ncol=p)
for (i in 1:p){
	table[i] = pcor.test(x[,i],y,x[,-i])$estimate}

####Make the partial correlation be in between 0 and 1
for (i in 1:p)	table[i] = (1+table[i])/2	
	
#####Adjacency Matrix with Cor Cent####
A = matrix(0,nrow=p,ncol=p)
colnames(A) = colname
rownames(A) = colname
for (j in 1:p){
	for (u in 1:p){
	if (j != u){
		A[j,u]=A[j,u]+sqrt((table[1,j])^2+(table[1,u])^2)}
		else{A[j,u]=A[j,u]}}}
c = matrix(0,p,1)
for (i in 1:p){
	c[i] = sum(A[-i,-i])/(2*(p-1)*(p-2))}
num = matrix(colname,p,1)
return(as.data.frame(list(Node =num, Weighted_Clustering = c)))
}


#A different function of partial correlation is used to construct the weight network
pcor.cluster_new <- function(x,y,b)
{
  
  p = dim(x)[2]
  colname = b
  #colnames(x) = colname
  ######Table for the correlation#######
  table  = matrix(0,nrow=1,ncol=p)
  for (i in 1:p){
    table[i] = pcor.test(x[,i],y,x[,-i])$estimate}
  
  ####Make the partial correlation be in between 0 and 1
  for (i in 1:p)	table[i] = (1+table[i])/2	
  
  ##New Metric
  for (i in 1:p) table[i] = sqrt(2*(1-table[i]^2))
  
  #####Adjacency Matrix with Cor Cent####
  A = matrix(0,nrow=p,ncol=p)
  colnames(A) = colname
  rownames(A) = colname
  for (j in 1:p){
    for (u in 1:p){
      if (j != u){
        A[j,u]=A[j,u]+table[1,j]+table[1,u]}
      else{A[j,u]=A[j,u]}}}
  c = matrix(0,p,1)
  for (i in 1:p){
    c[i] = sum(A[-i,-i])/(2*(p-1)*(p-2))}
  num = matrix(colname,p,1)
  return(as.data.frame(list(Node =num, Weighted_Clustering = c)))
}




################################################
##########Sequential Test#######################
################################################


#Following R functions describe the proposed clustering algorithm,
#the sequential testing algorithm, to identify clusters.
#Since we need a starting value to implement this algorithm,
#we begin with either maximum or minumum Network wide metric values



#Sequential Testing Method using the Maximum
sequential_max <- function(n,p,alpha,c,degree,degree_cov, idx, tau)
{
  df = n+n-2
  #tscore = qt(alpha/2,df)
  #print(tscore)
  temp = degree
  temp_2 = degree
  #print(temp)
  cluster = matrix(0,c,p+1)
  #print(cluster)
  for (i in 1:c){
    cluster[i,1]=i
    max_idx = which(degree==max(temp))
    temp_len = length(which(temp!=0))
    tscore = qt(alpha/(temp_len),df)
    #print(max_idx)
    #print(tscore)
    cluster[i,max_idx+1] = idx[max_idx]
    temp_2[max_idx]=0
    
    for (j in 1:p){
      #print(j)
      max_idx_2 = which(temp_2== max(temp_2))
      #print(max_idx_2)
      t = tstat(temp[max_idx], temp_2[max_idx_2], degree_cov[max_idx,max_idx], degree_cov[max_idx_2,max_idx_2], degree_cov[max_idx,max_idx_2],tau)
      #print(t)
      if (abs(t)<= abs(tscore)){
        cluster[i,max_idx_2+1] = idx[max_idx_2]
        temp_2[max_idx_2]=0
      }
      #print(temp_2)
      #print(sum(temp_2))
      if (sum(temp_2)==0) break
      if (abs(t)> abs(tscore)) break
    }
    #print(temp_2)
    temp = temp_2
    
    #print(temp)
    #print(h)
    if (sum(temp)==0){
      print("No more cluster!")
      break
    }
    else next
  }
  cluster_idx = matrix(0,1,p)
  for (i in 1:p) cluster_idx[i] = which(cluster[,-1]==i, arr.ind=TRUE)[1]
  return(cluster_idx)
}




#Sequential Testing Method using the Minimum
sequential_min <- function(n,p,alpha,c,degree,degree_cov, idx, tau)
{
  df = n+n-2
  temp = degree
  temp_2 = degree
  cluster = matrix(0,c,p+1)
  for (i in 1:c){
    cluster[i,1]=i
  	min_idx = which(degree==min(temp[temp!=0]))
  	temp_len = length(which(temp!=0))
  	tscore = qt(alpha/(2*temp_len),df)
  	cluster[i,min_idx+1] = idx[min_idx]
  	temp_2[min_idx]=0
  	if (sum(temp_2)==0) break
  	for (j in 1:p){
    	min_idx_2 = which(temp_2== min(temp_2[temp_2!=0]))
    	t = tstat(temp[min_idx], temp_2[min_idx_2], degree_cov[min_idx,min_idx], degree_cov[min_idx_2,min_idx_2], degree_cov[min_idx,min_idx_2], tau)
  		if (abs(t)<= abs(tscore)){
	  		cluster[i,min_idx_2+1] = idx[min_idx_2]
	  		temp_2[min_idx_2]=0
  		}
	  	if (abs(t)> abs(tscore)) break
	  	if (sum(temp_2)==0) break
	}
	temp = temp_2
    if (sum(temp)==0){
      print("No more cluster!")
      break
    }
    else next
  }
  return(cluster)
}


########################################
####Sequential Testing with maximum#####
########################################
sequential_new <- function(n,p,alpha,c,degree,degree_cov, idx, tau)
{
  df = n+n-2
  temp = degree
  cluster = matrix(0,c,p+1)
  for (i in 1:c){
    cluster[i,1]=i
    max_idx = which(degree==max(temp[temp!=0]))
	temp_len = length(which(temp!=0))
	tscore = qt(alpha/(2*temp_len),df)
    #print(max_idx)
    h = matrix(0,1,p)
    g = matrix(0,1,p)
    for (j in 1:p){
      #print(j)
      t = tstat(temp[max_idx], temp[j], degree_cov[max_idx,max_idx], degree_cov[j,j], degree_cov[max_idx,j], tau)
      #print(t)
      if (abs(t)<abs(tscore)){
        g[j]=j
        h[j]=temp[j]
      }
      else{
        g[j]=0
        h[j]=0
      }
    }
    for (j in 1:p)
		if (g[j]!=0) cluster[i,j+1]=idx[g[j]]
    #print(cluster[i,])
    for (j in 1:p) temp[j]=temp[j]-h[j]
    #print(temp)
    #print(h)
    if (sum(temp)==0){
      print("No more cluster!")
      break
    }
    else next
  }
  return(cluster)
}

#R function for the test statistic
tstat <- function(x,y,varx,vary,varxy, tau)
{
  top = abs(x-y)-tau
  bot = varx+vary-(2*varxy)
  if (bot==0) t=0
  else t = top/(bot^.5)
  return (t)
}
