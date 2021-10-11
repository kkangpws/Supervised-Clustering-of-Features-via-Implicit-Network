/*Brandon Woosuk Park
George Mason University*/

/*Functions for Simulation Studies in C*/

#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/*Sort a vector*/
void sorts(int n,double *x, int *y)
{
    int i,j;
	double *d1;
	d1 = dvector(1,n);
	
	for (i=1;i<=n;i++) d1[i] = x[i];
	sort(n,d1);
	for (i=1;i<=n;i++){
		for (j=1;j<=n;j++){
			if (d1[i]==x[j]) y[i]=j;
		}
	}
	free_dvector(d1,1,n);
}


/*Compute the average of a vector
Given n and x, the function returns the average of a vector*/
void average(int n, double *x, double *avg)
{
	int i;
	double sum=0;
	for (i=1;i<=n;i++){
		sum= sum + x[i];
	}
	
	*avg = sum/n;
}

/*Compute the sample covariance of two n x 1 vectors
Given n, x, and y, the function computes the sample covariance
If x=y, then the function computes the sample variance*/
void samplecov(int n, double *x, double *y, double *samvar)
{
	int i,j;
	double avg1, avg2,*d1,*d2, sum;
	d1 = dvector(1,n);
	d2 = dvector(1,n);
	average(n,x,&avg1);
	average(n,y,&avg2);
	
	for (i=1;i<=n;i++){
		d1[i] = x[i]-avg1;
		d2[i] = y[i]-avg2;
	}
	sum=0;
	for (j=1;j<=n;j++){
		sum = sum + (d1[j]*d2[j]);
	}
	
	*samvar = sum /(n-1);

	free_dvector(d1,1,n);
	free_dvector(d2,1,n);
}

/*Compute the sample variance of a vector*/
void samplevar(int n, double *x, double *var)
{
	int i;
	double avg, *d, sum;
	d = dvector(1,n);
	average(n,x,&avg);
	sum=0;
	for (i=1;i<=n;i++){
		d[i]=x[i]-avg;
	}
	for (i=1;i<=n;i++)
		sum = sum + pow(d[i],2);

	*var = sum/(n-1);
	
	free_dvector(d,1,n);
}

//Compute the correlation between two vectors
void cor(int n, double *x, double *y, double *corr)
{
	int i;
	double avgx, avgy;
	average(n,x,&avgx);
	average(n,y,&avgy);
	double sum1=0,sum2=0,sum3=0;
	for (i=1;i<=n;i++){
		sum1 = sum1 + (x[i]-avgx)*(y[i]-avgy);
		sum2 = sum2 + pow(x[i]-avgx,2);
		sum3 = sum3 + pow(y[i]-avgy,2);
	}
	
	*corr = sum1/(pow(sum2,.5)*pow(sum3,.5));

}


//Scale a matrix so that each column of the matrix has mean 0 and variance 1
void scale(int n,int p, double **x, double **x_scale)
{
	int i,j;
	double avg, var, *d;
	d = dvector(1,n);
	for (i=1;i<=p;i++){
		for (j=1;j<=n;j++){
			d[j] = x[j][i];
		}
		average(n,d,&avg);
		samplevar(n, d, &var);
		for (j=1;j<=n;j++){
			x_scale[j][i] = (d[j]-avg)/pow(var,.5);
		}
		avg=0;
		var=0;
	}
	free_dvector(d,1,n);
}

//Compute the t statistic for two sample t test
void ttest(int n,double beta, double beta_0, double *y, double *x, double *t)
{
	int i, df;
	double avgx,*y_fit, sum1=0, sum2=0, se;
	df = n-2;
	average(n,x,&avgx);
	y_fit = dvector(1,n);
	
	for (i=1;i<=n;i++) y_fit[i] = x[i]*beta;
	
	for (i=1;i<=n;i++){
		sum1 = sum1 + pow(y[i]-y_fit[i],2);
		sum2 = sum2 + pow(x[i]-avgx,2);
	}

	se =pow((1/(double)df)*sum1,.5)/pow(sum2,.5);
	*t = (beta-beta_0)/se;
	free_dvector(y_fit,1,n);
}

//Variable Selection using T-test
void varselection_t(int n, int p, double alpha, double *beta,
double beta_0, double *y, double **x, int *chosen_scad,
double *beta_t, int *sel_t)
{
	int i,j,df;
	double tscore, *dt, *t,tt;
	df = n-2;
	tscore=qt(alpha/2, df, 1, 0);
	dt = dvector(1,n);
	t = dvector(1,p);
	
	for (j=1;j<=p;j++){
		for (i=1;i<=n;i++){
			dt[i]=x[i][j];
		}
		ttest(n,beta[j],beta_0,y,dt,&tt);
		t[j]=tt;
	}

	for (i=1;i<= p;i++){
		if (fabs(t[i])> fabs(tscore)){
			sel_t[i]=chosen_scad[i];
			beta_t[i] = beta[i];
		}
		else {
			sel_t[i]=0;
			beta_t[i]=0;
		}
	}
	
	free_dvector(dt,1,n);
	free_dvector(t,1,p);
}


/*Compute the partial correlation between one variable
and all the other variables from the precision matrix.
For example, in our case, the function gives the
sample partial correlation between Y and Xs. However, the variable
for which the partial correlation is computed must be in
1st column of the precision matrix*/
void partialcorr(int p, double **x, double *b)
{
	int i;
	for (i=2;i<= p+1;i++){
		b[i-1] = - x[1][i]/sqrt(x[1][1]*x[i][i]);
	}
}

/*Computes the weight matrix based on the partial correlation.
This is the implicit network we are constructing.*/
void weight_matrix(int p, double *b, double **x)
{
	int i,j;
	for (i=1;i<=p;i++){
		for (j=1;j<=p;j++){
			if (i==j){
				x[i][j]=0;
			}
			else{
				x[i][j]=sqrt(pow(b[i],2)+pow(b[j],2));
			}
		}
	}
}

/*Computes the weight matrix based on the new-defined partial correlation.
This is the implicit network we are constructing.*/
void weight_matrix_new(int p, double *b, double **x)
{
	int i,j;
	for (i=1;i<=p;i++){
		for (j=1;j<=p;j++){
			if (i==j){
				x[i][j]=0;
			}
			else{
				x[i][j]=b[i]+b[j];
			}
		}
	}
}

/*Compute the p x p  covariance matrix of n x p matrix*/
void cov_mat(int n, int p, double **x, double **cov_x)
{
	int i, j,k;
	double *d1, *d2, samvar;
	d1 = dvector(1,n);
	d2 = dvector(1,n);
	for (j=1;j<=p;j++){
		for (k=1;k<=p;k++){
			for (i=1;i<=n;i++){
				d1[i]=x[i][j];
				d2[i]=x[i][k];
			}
		samplecov(n,d1,d2,&samvar);
		cov_x[j][k]=samvar;
		}
	}
	
	free_dvector(d1,1,n);
	free_dvector(d2,1,n);
}

//Summation of all entries in a p x p matrix
void cov_mat_sum(int p,double **cov_x, double *var_sum)
{
	int i,j;
	double sum=0;
	for (i=1;i<=p;i++){
		for (j=1;j<=p;j++){
			sum = sum + cov_x[i][j];
		}
	}
	*var_sum = sum/(p*p);
}

/*Compute the sample degree centrality based on the weights.
The degree centrality of a node is just the sum of the weight
in a row or a column.*/
//Compute degree centrality if we have weight matrix
void degree_weight(int p, double **weight, double *d)
{
	int i,j,k;
	double sum;
	for (i=1;i<=p;i++){
		sum=0;
		for (j=1;j<=p;j++){
			sum = sum + weight[i][j];
		}
		d[i]=sum;
	}
}

//Compute degree centrality if only data matrix is given
void degree(int n, int p, double **x, double *d)
{
	int i,j,k;
	//Find the Covaraince of p+1 variables
	double *d1, *d2, samvar, **var, **varinv;
	d1 = dvector(1,n);
	d2 = dvector(1,n);
	var = dmatrix(1,p+1,1,p+1);
	for (j=1;j<=p+1;j++){
		for (k=1;k<=p+1;k++){
			for (i=1;i<=n;i++){
				d1[i]=x[i][j];
				d2[i]=x[i][k];
			}
		samplecov(n,d1,d2,&samvar);
		var[j][k]=samvar;
		}
	}
	

	//Calculate the precision matrix
	varinv = dmatrix(1,p+1,1,p+1);
	dmatrix_inv(var,varinv,p+1);

	//Compute the vector of the partial correlation coefficients
	double *parcor;
	parcor = dvector(1,p);
	partialcorr(p,varinv,parcor);

	//Compute the weight matrix
	double **weight;
	weight = dmatrix(1,p,1,p);
	weight_matrix(p,parcor,weight);

	double sum;
	for (i=1;i<=p;i++){
		sum=0;
		for (j=1;j<=p;j++){
			sum = sum + weight[i][j];
		}
		d[i]=sum;
	}
	
	free_dvector(d1,1,n);
	free_dvector(d2,1,n);
	free_dmatrix(var,1,p+1,1,p+1);
	free_dmatrix(varinv,1,p+1,1,p+1);
	free_dvector(parcor,1,p);
	free_dmatrix(weight,1,p,1,p);
}

/*Compute the sample clustering coefficient based on the weights.*/
void clustering_coef(int n, int p, double **x, double *c)
{
	
	int i,j,k;
	//Find the Covaraince of p+1 variables
	double *d1, *d2, samvar, **var, **varinv;
	d1 = dvector(1,n);
	d2 = dvector(1,n);
	var = dmatrix(1,p+1,1,p+1);
	for (j=1;j<=p+1;j++){
		for (k=1;k<=p+1;k++){
			for (i=1;i<=n;i++){
				d1[i]=x[i][j];
				d2[i]=x[i][k];
			}
		samplecov(n,d1,d2,&samvar);
		var[j][k]=samvar;
		}
	}
	
	//Calculate the precision matrix
	varinv = dmatrix(1,p+1,1,p+1);
	dmatrix_inv(var,varinv,p+1);
	
	//Compute the vector of the partial correlation coefficients
	double *parcor;
	parcor = dvector(1,p);
	partialcorr(p,varinv,parcor);
	
	//Compute the weight matrix
	double **weight;
	weight = dmatrix(1,p,1,p);
	weight_matrix(p,parcor,weight);

	double sum,l;
	l = (p-1)*(p-2);
	
	for (i=1;i<=p;i++){
		sum=0;
		for (j=1;j<=p;j++){
			for (k=1;k<=p;k++){
				if (i==j || k==i){
					continue;
				}
				else{
					sum = sum + weight[j][k];
				}
			}
		}
		c[i]=sum/(2*l);
	}
	free_dvector(d1,1,n);
	free_dvector(d2,1,n);
	free_dmatrix(var,1,p+1,1,p+1);
	free_dmatrix(varinv,1,p+1,1,p+1);
	free_dvector(parcor,1,p);
	free_dmatrix(weight,1,p,1,p);
}

//Compute clustering coefficient if we have weight matrix
void clustering_coef_weight(int p, double **weight, double *c)
{
	int i,j,k;
	double sum,l;
	l = (p-1)*(p-2);
	for (i=1;i<=p;i++){
		sum=0;
		for (j=1;j<=p;j++){
			for (k=1;k<=p;k++){
				if (i==j || k==i){
					continue;
				}
				else{
					sum = sum + weight[j][k];
				}
			}
		}
		c[i]=sum/(2*l);
	}
}


/*Now we try to use bootstrap method to find the sample distribution
of the network metrics.*/

//Index number generator
void bt(int n, int bn, long *idum, int *index)
{
    int i;
    
    for (i=1; i<=bn; i++)
        index[i] = (int)ceil(ran2(idum)*n);
}

//Compute bootstrap network wide metrics
void boot_cent(int n, int nbt, int p, double **x,
double **boot_parcor, double **boot_degree,
double **boot_cluster)
{
	int i,j,k,l,m, bn;
	bn=n;
	long idum=(-99);
	
	for (i=1; i<=nbt;i++){
		int *index;
		double **boot_sample;
		index = ivector(1,bn);
		boot_sample = dmatrix(1,bn,1,p+1);
		//Index vector and bootstrap sample vector
		bt(n, bn, &idum, index);
		
		for (j=1;j<=bn;j++){
			for (k=1; k<=p+1;k++){
				boot_sample[j][k] = x[index[j]][k];
			}
		}
		
		//Bootstrap Partial Correlation
		double *d1, *d2, samvar, **var, **varinv;
		d1 = dvector(1,bn);
		d2 = dvector(1,bn);
		var = dmatrix(1,p+1,1,p+1);
		for (j=1;j<=p+1;j++){
			for (k=1;k<=p+1;k++){
				for (l=1;l<=bn;l++){
					d1[l]=boot_sample[l][j];
					d2[l]=boot_sample[l][k];
				}
			samplecov(bn,d1,d2,&samvar);
			var[j][k]=samvar;
			}
		}

		//Calculate the precision matrix
		varinv = dmatrix(1,p+1,1,p+1);
		dmatrix_inv(var,varinv,p+1);

		//Compute the vector of the partial correlation coefficients
		double *parcor;
		parcor = dvector(1,p);
		partialcorr(p,varinv,parcor);
		for (k=1;k<=p;k++) parcor[k]= (1+parcor[k])/2;
		for (k=1;k<=p;k++) boot_parcor[i][k] = parcor[k];

		//Compute the weight matrix
		double **weight;
		weight = dmatrix(1,p,1,p);
		weight_matrix(p,parcor,weight);

		//Bootstrap degree centrality
		double *b_degree, b_degree_avg;
		b_degree = dvector(1,p);
		degree_weight(p,weight,b_degree);
		for (j=1;j<=p;j++) boot_degree[i][j] = b_degree[j];
		
		
		//Compute the clustering coefficients
		double *b_cluster, b_cluster_avg;
		b_cluster = dvector(1,p);
		clustering_coef_weight(p,weight,b_cluster);
		for (j=1;j<=p;j++) boot_cluster[i][j] = b_cluster[j];


		free_dvector(d1,1,bn);
		free_dvector(d2,1,bn);
		free_dmatrix(var,1,p+1,1,p+1);
		free_dmatrix(varinv,1,p+1,1,p+1);
		free_dvector(parcor,1,p);
		free_dmatrix(weight,1,p,1,p);
		free_dvector(b_degree,1,p);
		free_dvector(b_cluster,1,p);
		free_ivector(index,1,bn);
		free_dmatrix(boot_sample,1,bn,1,p+1);
		
		}
	
}

//Bootstrap and network wide metrics using new metric
void boot_cent_new(int n, int nbt, int p, double **x,
double **boot_parcor, double **boot_degree,
double **boot_cluster)
{
	int i,j,k,l,m, bn;
	bn=n;
	long idum=(-99);
	
	for (i=1; i<=nbt;i++){
		int *index;
		double **boot_sample;
		index = ivector(1,bn);
		boot_sample = dmatrix(1,bn,1,p+1);
		//Index vector and bootstrap sample vector
		bt(n, bn, &idum, index);
		
		for (j=1;j<=bn;j++){
			for (k=1; k<=p+1;k++){
				boot_sample[j][k] = x[index[j]][k];
			}
		}
		
		//Bootstrap Partial Correlation
		double *d1, *d2, samvar, **var, **varinv;
		d1 = dvector(1,bn);
		d2 = dvector(1,bn);
		var = dmatrix(1,p+1,1,p+1);
		for (j=1;j<=p+1;j++){
			for (k=1;k<=p+1;k++){
				for (l=1;l<=bn;l++){
					d1[l]=boot_sample[l][j];
					d2[l]=boot_sample[l][k];
				}
			samplecov(bn,d1,d2,&samvar);
			var[j][k]=samvar;
			}
		}

		//Calculate the precision matrix
		varinv = dmatrix(1,p+1,1,p+1);
		dmatrix_inv(var,varinv,p+1);

		//Compute the vector of the partial correlation coefficients
		double *parcor;
		parcor = dvector(1,p);
		partialcorr(p,varinv,parcor);
		//for (k=1;k<=p;k++) parcor[k]= (1+parcor[k])/2;
		
		
		
		//for (k=1;k<=p;k++) boot_parcor[i][k] = parcor[k];
		
		for (k=1;k<=p;k++) parcor[k]= sqrt(2*(1-pow(parcor[k],2)));

		for (k=1;k<=p;k++) boot_parcor[i][k] = parcor[k];
		//Compute the weight matrix
		double **weight;
		weight = dmatrix(1,p,1,p);
		weight_matrix_new(p,parcor,weight);
		
		
		//Bootstrap degree centrality
		double *b_degree, b_degree_avg;
		b_degree = dvector(1,p);
		degree_weight(p,weight,b_degree);
		for (j=1;j<=p;j++) boot_degree[i][j] = b_degree[j];
		
		
		//Compute the clustering coefficients
		double *b_cluster, b_cluster_avg;
		b_cluster = dvector(1,p);
		clustering_coef_weight(p,weight,b_cluster);
		for (j=1;j<=p;j++) boot_cluster[i][j] = b_cluster[j];
		
	
		
		free_dvector(d1,1,bn);
		free_dvector(d2,1,bn);
		free_dmatrix(var,1,p+1,1,p+1);
		free_dmatrix(varinv,1,p+1,1,p+1);
		free_dvector(parcor,1,p);
		free_dmatrix(weight,1,p,1,p);
		free_dvector(b_degree,1,p);
		free_dvector(b_cluster,1,p);
		free_ivector(index,1,bn);
		free_dmatrix(boot_sample,1,bn,1,p+1);
		
		}
	
}


//Randomly split the data into two parts.

void rand_split(int n, int *index1, int *index2)
{
	//const int n; //amount of random numbers that need to be generated
	//const int n; //maximum value (of course, this must be at least the same as AMOUNT;
	int i;
	int *value; //array to store the random numbers in
	double *d;
	value = ivector(1,n);
	d = dvector(1,n);
	//srand(1234); //always seed your RNG before using it
	for (i=1;i<=n;i++) d[i]=runif(0,1);
	
	sorts(n,d,value);
	
	for (i=1;i<=n/2;i++){
		index1[i]=value[i];
	}

	for (i=n/2 + 1;i<=n;i++){
		index2[i-(n/2)]=value[i];
	}

	free_dvector(d,1,n);
	free_ivector(value,1,n);
}