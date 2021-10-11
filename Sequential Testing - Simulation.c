/*Brandon Woosuk Park
George Mason University*/

/*Main C Code for simulations to detect communities based on
network wide metrics using sequential testing method
Here we focus on the degree centrality*/

//This files corresponds with the Table 3 in the manuscript.


//Simulation Procedure
/*
1. Generate samples
2. Using multiple splits (for the time being there are 20 splits), select the active set of variables based on the frequency
3. With the selected variables, use bootstrap technique to derive the distribution of the network wide metric of selected variables
4. Detect the communities based on Tau and test the difference between clusters
*/


#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"
#include "NRCRoutines.c"
#include "26function.c"
#include "mcpfunction.c"
#include "clusterfunction.c"

int main(){
	int n,p,i,j,b,l,c,m,nsim,k;
	double *beta;
	int **degcluster,**clustercoef;
	
	set_seed(41,42);
	nsim=10; //Number of simulation studies
	n=200; //Number of observations
	p=20; //Number of parameters
	c=3; //Desired number of clusters
	//Create true betas
	beta=dvector(1,p);
	for (i=1; i<= 3;i++) beta[i]=1;	
	for (i=4; i<= 6;i++) beta[i]=-1;
	for (i=7; i<= 9;i++) beta[i]=2;
	for (i=10; i<= p;i++) beta[i]=0;

	//Arrays to store clusters
	degcluster = imatrix(1,nsim*c,1,p+2);
	clustercoef = imatrix(1,nsim*c,1,p+2);
	for (i=1;i<=nsim*c;i++){
		for (j=1;j<=p+2;j++){
			degcluster[i][j]= 0;
		}
	}
	

	for (i=1;i<=nsim*c;i++){
		for (j=1;j<=p+2;j++){
			clustercoef[i][j]= 0;
		}
	}
	
	
	//Simulation begins
	for (m=1;m<=nsim;m++){
		printf("%d \n", m);
	//Generate Samples
	double **x, *y, *eps, mu, **x_scale;
	x = dmatrix(1,n,1,p);
	x_scale = dmatrix(1,n,1,p);
	y = dvector(1,n);
	eps = dvector(1,n);
	mu=0;
	
	//Generating X
	for (i=1; i<=n;i++){
		for (j=1;j<=p;j++){
			x[i][j]=rnorm(0,1)+mu;
		}
	}

	//Generating epsilons
	for (i=1; i<=n;i++){
		eps[i]=rnorm(0,1);
	}

	//Generating Y
	double sum;
	for (i=1;i<=n;i++){
		sum=0;
		for (j=1;j<=p;j++){
			sum = sum+ x[i][j]*beta[j];
		}
		y[i]=sum + eps[i];
	}

	//Scale X
	scale(n,p,x,x_scale);
	
	//Multiple splits
	l=1; //Number of Multiple splits
	int **dcluster,**clusterc,**closec;
	dcluster = imatrix(1,l*c,1,p+1);
	clusterc = imatrix(1,l*c,1,p+1);
	closec = imatrix(1,l*c,1,p+1);
	
	//Set 'variables' and 'degree_split' equal to 0
	for (i=1;i<=l*c;i++){
		for (j=1;j<=p+1;j++) dcluster[i][j]=0;
	}

	
	for (i=1;i<=l*c;i++){
		for (j=1;j<=p+1;j++) clusterc[i][j]=0;
	}
	
	for (b=1;b<=l;b++){
		//Split the samples into two parts
		int nn,*idx1, *idx2;
		nn=n/2;
		double *y_test, *y_train, **x_test, **x_train;
		idx1 = ivector(1,nn);
		idx2 = ivector(1,nn);
		rand_split(n,idx1,idx2);

		y_test = dvector(1,nn);
		y_train = dvector(1,nn);
		x_test = dmatrix(1,nn,1,p);
		x_train = dmatrix(1,nn,1,p);
		
		for (i=1;i<=nn;i++){
		y_train[i] = y[idx1[i]];
		y_test[i] = y[idx2[i]];
		}
		
		for (i=1;i<=nn;i++){
			for (j=1;j<=p;j++){
				x_train[i][j]=x_scale[idx1[i]][j];
				x_test[i][j]=x_scale[idx2[i]][j];
			}
		}
		
		//Variable selection vis MCP with k-fold cross-validation
		int k=10;
		int *var_chosen;
		double *beta_mcp;
		var_chosen = ivector(1,p);
		beta_mcp = dvector(1,p);
		
		var_selection(k,n/2,p,2,0,50,4,1,30,1000,pow(10,-8),y_train, x_train, var_chosen, beta_mcp);
		
		
		//Reduce the x_test by choosing only selected variables
		int len_var=0;
		for (i=1;i<=p;i++){
		if (var_chosen[i]!=0) len_var = len_var + 1;
		else continue;
		}
		pushZerosToEnd(p,var_chosen);

		
		double **x_test_mcp, *beta_test_mcp;
		x_test_mcp= dmatrix(1,nn,1,len_var);
		beta_test_mcp = dvector(1,len_var);
		for (i=1;i<=nn;i++){
			for (j=1;j<=len_var;j++){
				x_test_mcp[i][j] = x_test[i][var_chosen[j]];
			}
		}
		for (i=1;i<=len_var;i++) beta_test_mcp[i]=beta_mcp[var_chosen[i]];

		
		//Do t-test for further variable selection
		double *beta_t;
		int *sel_t;
		beta_t = dvector(1,len_var);
		sel_t = ivector(1,len_var);
		varselection_t(nn,len_var,.15, beta_test_mcp,0,y_test, x_test_mcp,
		var_chosen, beta_t, sel_t);
		

		//Reduce the x_test by choosing only selected variables
		int len_var_t=0;
		for (i=1;i<=len_var;i++){
			if (sel_t[i]!=0) len_var_t = len_var_t + 1;
			else continue;
		}
		pushZerosToEnd(len_var,sel_t);

		
		//select X's from variable selection
		double **x_test_t;
		x_test_t = dmatrix(1,nn,1,len_var_t);
	
		for (i=1;i<=nn;i++){
			for (j=1;j<=len_var_t;j++){
				x_test_t[i][j]=x_test[i][sel_t[j]];
			}
		}

		//Matrix X only with selected variables
		double **data;
		data = dmatrix(1,nn,1,len_var_t+1);
	
		for (i=1;i<=nn;i++){
			data[i][1]=y_test[i];
			for (j=2;j<=len_var_t+1;j++){
				data[i][j] = x_test_t[i][j-1];
			}
		}

		//Compute the covariance matrix of p+1 variables
		double *d1, *d2, samvar, **var, **varinv;
		d1 = dvector(1,nn);
		d2 = dvector(1,nn);
		var = dmatrix(1,len_var_t+1,1,len_var_t+1);
		for (j=1;j<=len_var_t+1;j++){
			for (k=1;k<=len_var_t+1;k++){
				for (i=1;i<=nn;i++){
					d1[i]=data[i][j];
					d2[i]=data[i][k];
				}
			samplecov(nn,d1,d2,&samvar);
			var[j][k]=samvar;
			}
		}
		
		//Calculate the precision matrix
		varinv = dmatrix(1,len_var_t+1,1,len_var_t+1);
		dmatrix_inv(var,varinv,len_var_t+1);
		
		//Compute the vector of the partial correlation coefficients
		double *parcor;
		parcor = dvector(1,len_var_t);
		partialcorr(len_var_t,varinv,parcor);

				
		//Make partial correlations between 0 and 1
		for (i=1;i<=len_var_t;i++) parcor[i]= (1+parcor[i])/2;
		
		//New defined weight between Y and X_i
		for (i=1;i<=len_var_t;i++) parcor[i]= sqrt(2*(1-pow(parcor[i],2)));

		//Compute the weight matrix
		double **weight;
		weight = dmatrix(1,len_var_t,1,len_var_t);
		weight_matrix_new(len_var_t,parcor,weight);
	
		//Compute the degree centrality
		double *degree_cent, degavg;
		degree_cent = dvector(1,len_var_t);
		degree_weight(len_var_t,weight,degree_cent);

		
		//Compute the clustering coefficient
		double *clustering_coef, clusteravg;
		clustering_coef = dvector(1,len_var_t);
		clustering_coef_weight(len_var_t,weight,clustering_coef);
		

		//Use Bootstrap to compute the mean and variance of each degree centrality
		//Use only selected variables from variable selection
		double **x_boot, **data_boot;
		x_boot = dmatrix(1,n,1,len_var_t);
		for (i=1;i<=n;i++){
			for (j=1;j<=len_var_t;j++){
				x_boot[i][j]=x[i][sel_t[j]];
			}
		}
	
		data_boot = dmatrix(1,n,1,len_var_t+1);
		for (i=1;i<=n;i++){
			data_boot[i][1]= y[i];
			for (j=2;j<=len_var_t+1;j++){
				data_boot[i][j]=x_boot[i][j-1];
			}
		}

		int nbt;
		double **boot_degree, **boot_cluster, 
		**boot_parcor,bootavgdeg, bootavgcluster, 
		bootavgclose, bootavgdegvar, bootavgclustervar, bootavgclosevar;
		nbt=500;
		boot_parcor = dmatrix(1,nbt,1,len_var_t);
		boot_degree = dmatrix(1,nbt,1,len_var_t);
		boot_cluster = dmatrix(1,nbt,1,len_var_t);
		boot_cent_new(n,nbt,len_var_t,data_boot,boot_parcor,boot_degree,
		boot_cluster);

		//Distribution of degree centrality
		double **deg_cov, *deg_avg_vec;
		deg_avg_vec = dvector(1,len_var_t);
		deg_cov = dmatrix(1,len_var_t,1,len_var_t);
		avg_vec(nbt,len_var_t,boot_degree,deg_avg_vec);	
		cov_mat(nbt,len_var_t,boot_degree,deg_cov);

		
		//Distribution of clustering coefficient
		double **cc_cov, *cc_avg_vec;
		cc_avg_vec = dvector(1,len_var_t);
		cc_cov = dmatrix(1,len_var_t,1,len_var_t);
		avg_vec(nbt,len_var_t,boot_cluster,cc_avg_vec);	
		cov_mat(nbt,len_var_t,boot_cluster,cc_cov);
	

		//Define the clusters using sequential testing
		int **degree_cluster;
		degree_cluster = imatrix(1,c,1,len_var_t+1);
		for (i=1;i<=c;i++){
			for (j=1;j<=len_var_t+1;j++){
				degree_cluster[i][j]=0;
			}
		}

		
		//Define the clusters using sequential testing
		int **cluster_coef;
		cluster_coef = imatrix(1,c,1,len_var_t+1);
		for (i=1;i<=c;i++){
			for (j=1;j<=len_var_t+1;j++){
				cluster_coef[i][j]=0;
			}
		}
		
		int multi_test; //Bonferroni Adjustment
		multi_test =  (len_var_t*(len_var_t-1))/2;

		//Process of the sequential tests
		diff_test_boot_new(n,len_var_t,c,.1/(double) multi_test,degree_cent,deg_cov,degree_cluster);
		diff_test_boot_new(n,len_var_t,c,.1/(double) multi_test,clustering_coef,cc_cov,cluster_coef);

		for (i=1;i<=c;i++){
			for (j=1;j<=len_var_t+1;j++){
				dcluster[(b-1)*c+i][j]=degree_cluster[i][j];
			}
		}

		for (i=1;i<=c;i++){
			for (j=1;j<=len_var_t+1;j++){
				clusterc[(b-1)*c+i][j]=cluster_coef[i][j];
			}
		}

		free_ivector(idx1,1,nn);
		free_ivector(idx2,1,nn);
		free_dvector(y_test,1,nn);
		free_dvector(y_train,1,nn);
		free_dmatrix(x_test,1,nn,1,p);
		free_dmatrix(x_train,1,nn,1,p);
		free_ivector(var_chosen,1,p);
		free_dvector(beta_mcp,1,p);
		free_dmatrix(x_test_mcp,1,nn,1,len_var);
		free_dvector(beta_test_mcp,1,len_var);
		free_dvector(beta_t,1,len_var);
		free_ivector(sel_t,1,len_var);
		free_dmatrix(x_test_t,1,nn,1,len_var_t);
		free_dmatrix(data,1,nn,1,len_var_t+1);
		free_dvector(d1,1,nn);
		free_dvector(d2,1,nn);
		free_dmatrix(var,1,len_var_t+1,1,len_var_t+1);
		free_dmatrix(varinv,1,len_var_t+1,1,len_var_t+1);
		free_dvector(parcor,1,len_var_t);
		free_dmatrix(weight,1,len_var_t,1,len_var_t);
		free_dvector(degree_cent,1,len_var_t);
		free_dvector(clustering_coef,1,len_var_t);
		free_dmatrix(x_boot,1,n,1,len_var_t);
		free_dmatrix(data_boot,1,n,1,len_var_t+1);
		free_dmatrix(boot_parcor,1,nbt,1,len_var_t);
		free_dmatrix(boot_degree,1,nbt,1,len_var_t);
		free_dmatrix(boot_cluster,1,nbt,1,len_var_t);
		free_dvector(deg_avg_vec,1,len_var_t);
		free_dmatrix(deg_cov,1,len_var_t,1,len_var_t);
		free_dvector(cc_avg_vec,1,len_var_t);
		free_dmatrix(cc_cov,1,len_var_t,1,len_var_t);
		free_imatrix(degree_cluster,1,c,1,len_var_t+1);
		free_imatrix(cluster_coef,1,c,1,len_var_t+1);
	}

	//Create a matrix to count the number of apprearnces
	double **count_mat_deg, **count_mat_cluster;
	count_mat_deg = dmatrix(1,c,1,p+1);
	count_mat_cluster = dmatrix(1,c,1,p+1);
	for (i=1;i<=c;i++){
		count_mat_deg[i][1]=i;
		count_mat_cluster[i][1]=i;
		for (j=2;j<=p+1;j++){
			count_mat_deg[i][j]=0;
			count_mat_cluster[i][j]=0;
		} 
	}
	

	for (b=1;b<=l;b++){
		for (i=1;i<=c;i++){
			if (dcluster[(b-1)*c + i][1]==i){
				for (j=2;j<=p+1;j++){
					if 	(dcluster[(b-1)*c + i][j]!=0) count_mat_deg[i][dcluster[(b-1)*c + i][j]+1]= count_mat_deg[i][dcluster[(b-1)*c + i][j]+1]+1;
				}
			}
		}
	}
	
	
	for (b=1;b<=l;b++){
		for (i=1;i<=c;i++){
			if (clusterc[(b-1)*c + i][1]==i){
				for (j=2;j<=p+1;j++){
					if 	(clusterc[(b-1)*c + i][j]!=0) count_mat_cluster[i][clusterc[(b-1)*c + i][j]+1]= count_mat_cluster[i][clusterc[(b-1)*c + i][j]+1]+1;
				}
			}
		}
	}
	

	for (i=1;i<=c;i++){
		for (j=2;j<=p+1;j++){
			count_mat_deg[i][j] = count_mat_deg[i][j]/l;
			count_mat_cluster[i][j] = count_mat_cluster[i][j]/l;
		}
	}
	

	//Select variables that appear the most in the first, second, third cluster, respectively
	for (i=1;i<=c;i++){
		for (j=2;j<=p+1;j++){
			if (count_mat_deg[i][j]>= .5) count_mat_deg[i][j]=j-1;
			else count_mat_deg[i][j]=0;
		}
	}

	for (i=1;i<=c;i++){
		for (j=2;j<=p+1;j++){
			if (count_mat_cluster[i][j]>= .5) count_mat_cluster[i][j]=j-1;
			else count_mat_cluster[i][j]=0;
		}
	}
	

	//Push zeroes to the end
	for (i=1;i<=c;i++){
		int *temp;
		temp = ivector(1,p+1);
		for (j=1;j<=p+1;j++) temp[j]=count_mat_deg[i][j];
		pushZerosToEnd(p+1,temp);
		for (j=1;j<=p+1;j++) count_mat_deg[i][j]=temp[j];
		
		free_ivector(temp,1,p+1);
	}

	for (i=1;i<=c;i++){
		int *temp;
		temp = ivector(1,p+1);
		for (j=1;j<=p+1;j++) temp[j]=count_mat_cluster[i][j];
		pushZerosToEnd(p+1,temp);
		for (j=1;j<=p+1;j++) count_mat_cluster[i][j]=temp[j];
		
		free_ivector(temp,1,p+1);
	}
	
	
	for (i=1;i<=c;i++){
		degcluster[(m-1)*c+i][1]=m;
		for (j=2;j<=p+2;j++){
			degcluster[(m-1)*c+i][j]=count_mat_deg[i][j-1];
		}
	}

	
	for (i=1;i<=c;i++){
		clustercoef[(m-1)*c+i][1]=m;
		for (j=2;j<=p+2;j++){
			clustercoef[(m-1)*c+i][j]=count_mat_cluster[i][j-1];
		}
	}
	
	free_dvector(eps,1,n);
	free_dmatrix(x,1,n,1,p);
	free_dvector(y,1,n);
	free_dmatrix(x_scale,1,n,1,p);
	free_imatrix(dcluster,1,l*c,1,p+1);
	free_imatrix(clusterc,1,l*c,1,p+1);
	free_dmatrix(count_mat_deg,1,c,1,p+1);
	free_dmatrix(count_mat_cluster,1,c,1,p+1);
	}

	//Create output files
	FILE *fp;
	fp = fopen("cluster_degree.csv","w");
	
	for(j=1;j<=nsim*c;j++){
		for (i=1;i<=p+2;i++){
			fprintf(fp,"%d \t", degcluster[j][i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	
	FILE *fp2;
	fp2 = fopen("cluster_cc.csv","w");
	
	for(j=1;j<=nsim*c;j++){
		for (i=1;i<=p+2;i++){
			fprintf(fp2,"%d \t", clustercoef[j][i]);
		}
		fprintf(fp2,"\n");
	}
	fclose(fp2);
	
	free_dvector(beta,1,p);
	free_imatrix(degcluster,1,nsim*c,1,p+2);
	free_imatrix(clustercoef,1,nsim*c,1,p+2);

}



