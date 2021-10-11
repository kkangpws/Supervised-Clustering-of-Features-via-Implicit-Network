/*Brandon Woosuk Park
George Mason University*/

/*Function for MCP regularization*/

//Compute the maximum difference between two vectors
void max_diff(int n, double *x, double *y, double *max_dif)
{
	int i,j;
	double max, *d;
	d=dvector(1,n);
	for (i=1;i<=n;i++) d[i]=fabs(x[i]-y[i]);
	
	max = d[1];
	for (i=2;i<=n;i++){
		if (d[i]>max) max = d[i];
	}
	
	*max_dif = max;
	
	free_dvector(d,1,n);
}

//routine to obtain the OLS from a linear regression model
void ols(int n,int p,double *y, double **x, double *beta)
{
	int i,j;
	double **a, *b;
	
	a = dmatrix(1,p,1,p);
	b = dvector(1,p);
	
	//a = x'x
	dmatrix_tXX(n,p, x,a);
	//b= x'y
	dmatrix_tXY(n,p,x,y,b);
	
	//beta - (x'x)^(-1)x'y
	// apply the routine solving linear equations
	solve_linear_equation(a,b,p,beta);
	
	free_dmatrix(a,1,p,1,p);
	free_dvector(b,1,p);
}

//S function used in MCP
double s(double lambda, double z){
	if (z > lambda) return z-lambda;
	if (z < -lambda) return z+lambda;
	if (fabs(z)<= lambda) return 0;
}

//MCP penalty function
void mcp_penalty(double lambda, double gamma, double z, double *b_mcp)
{
	if (fabs(z)<= gamma*lambda) *b_mcp= s(lambda,z)/(1-(1/gamma));
	if (fabs(z)> gamma*lambda) *b_mcp= z;
}

//Main MCP function
void mcp(int n, int p,double gamma, double lambda, int max_iter,
double tol, double **x, double *y, double *b_mcp)
{
	int i,j,k;
	double *b_ols;
	b_ols=dvector(1,p);
	ols(n,p,y,x,b_ols);

	
	double *r;
	r = dvector(1,n);
	
	for (i=1;i<=n;i++){
		double sum=0;
		for (j=1;j<=p;j++){
			sum = sum + x[i][j]*b_ols[j];
		}
		r[i]=y[i]-sum;
	}

	double *b_init, *b_end;
	b_init = dvector(1,p);
	b_end = dvector(1,p);
	
	for (i=1;i<=p;i++) b_init[i]=b_ols[i];
	
	for (i=1;i<=max_iter;i++){
	for (j=1;j<=p;j++){
		double z_j=0, z_sum=0,*x_j; 
		x_j = dvector(1,n);
		for (k=1;k<=n;k++) x_j[k]=x[k][j];
		
		for (k=1;k<=n;k++) z_sum= z_sum + x_j[k]*r[k];
		
		z_j = z_sum/n + b_init[j];
		
		double b_update;
		mcp_penalty(lambda,gamma, z_j, &b_update);
		b_end[j]=b_update;
		
		for (k=1;k<=n;k++) r[k]=r[k]-(b_end[j]-b_init[j])*x_j[k];
		
		free_dvector(x_j,1,n);
	}

	double max;
	max_diff(p,b_init,b_end,&max);
	if (max < tol){
		for (k=1;k<=p;k++){
			b_mcp[k]=b_end[k];
		}
		break;
	}
	for (k=1;k<=p;k++) b_init[k]=b_end[k];	
	

	}
	free_dvector(b_ols,1,p);
	free_dvector(r,1,n);
	free_dvector(b_init,1,p);
	free_dvector(b_end,1,p);
}


/*Compute the mean squared error from the MCP coefficients*/
void mse(int n, int p, double *beta, double *y, double **x, double *merror)
{
	int i,j;
	double *y_fit, sum, *error;
	
	y_fit = dvector(1,n);
	error = dvector(1,n);
	//Compute the fitted y
	for (i=1;i<=n;i++){
		sum=0;
		for (j=1;j<=p;j++){
			sum = sum+ x[i][j]*beta[j];
		}
		y_fit[i]=sum;
	}

	
	for (i=1;i<=n;i++) error[i] = pow(y_fit[i]-y[i],2);

	double s=0;
	for (i=1;i<=n;i++){
		s = s + error[i];
	}
	
	*merror = s/(double)n;
	
	free_dvector(y_fit,1,n);
	free_dvector(error,1,n);
}

/*Random split the data into k parts*/
//Generate random numbers
void rand_num(int n, int *value)
{
	int i,j;
	double *d;
	d = dvector(1,n);
	//srand(1234); //always seed your RNG before using it
	for (i=1;i<=n;i++) d[i]=runif(0,1);
	
	sorts(n,d,value);
}

//Split data randomly into k groups
void k_split(int k,int n, int l,int p,int *value,double *y, double **x,
double *y_ktrain, double *y_ktest, double **x_ktrain, double **x_ktest)
{
	if ((n%k) != 0){
		printf("Data cannot be evenly divided into k parts!");
		exit(0);
	}
	else{
	int i,j,t;
	//Generate y
	if (l==1){
		for (i=1;i<=n/k;i++) y_ktest[i]= y[value[i]];
		for (i=1;i<=n-(n/k);i++) y_ktrain[i]=y[value[i+(n/k)]];
	}
	if (l==k){
		for (i=1;i<=n/k;i++) y_ktest[i]=y[value[i+(l-1)*(n/k)]];
		for (i=1;i<=n-(n/k);i++) y_ktrain[i]=y[value[i]];
	}
	else{
		for (j=1;j<=n-(n/k);j++){
			if (j <= (l-1)*(n/k)){
				y_ktrain[j]=y[value[j]];
			}
			else{
				y_ktrain[j]=y[value[j+(n/k)]];
			}
		}
		for (j=1;j<= n/k;j++){
				y_ktest[j]=y[value[j+(l-1)*(n/k)]];
			}
	}

	
	//Generate X
	for (t=1;t<=p;t++){
	if (l==1){
		for (i=1;i<=n/k;i++) x_ktest[i][t]= x[value[i]][t];
		for (i=1;i<=n-(n/k);i++) x_ktrain[i][t]=x[value[i+(n/k)]][t];
	}
	if (l==k){
		for (i=1;i<=n/k;i++) x_ktest[i][t]=x[value[i+(l-1)*(n/k)]][t];
		for (i=1;i<=n-(n/k);i++) x_ktrain[i][t]=x[value[i]][t];
	}
	else{
		for (j=1;j<=n-(n/k);j++){
			if (j <= (l-1)*(n/k)){
				x_ktrain[j][t]=x[value[j]][t];
			}
			else{
				x_ktrain[j][t]=x[value[j+(n/k)]][t];
			}
		}
		for (j=1;j<= n/k;j++){
				x_ktest[j][t]=x[value[j+(l-1)*(n/k)]][t];
			}
	}
	}
	

	}
}


//Compute MSE based on k-fold Cross validation 
void cv_mse(int k, int n, int p, double gamma, double lambda, 
int max_iter, double tol,double *y_train, double **x_train,
double *m)
{
	int i,l;
	double *msevec;
	msevec = dvector(1,k);
	//k fold cross validation
	int *value;
	value = ivector(1,n);
	srand(123); 
	rand_num(n,value);

	
	for (l=1;l<=k;l++){
		double *y_cvtrain, *y_cvtest, **x_cvtrain, **x_cvtest;
		y_cvtrain = dvector(1,n - n/k);
		y_cvtest = dvector(1,n/k);
		x_cvtrain = dmatrix(1,n - n/k,1,p);
		x_cvtest = dmatrix(1,n/k,1,p);	
		
		
		k_split(k,n,l, p ,value,y_train, x_train, y_cvtrain, y_cvtest,
		x_cvtrain, x_cvtest); 
		
		//MCP using the training set from lth cross validation
		double *beta_mcp;
		beta_mcp = dvector(1,p);
		mcp(n - n/k,p, gamma, lambda,max_iter,tol, x_cvtrain, y_cvtrain, beta_mcp);

		
		//Compute the MSE from lth cross validation
		double merror;
		mse(n/k,p,beta_mcp,y_cvtest,x_cvtest,&merror);
		
		//Plug the value into the MSE vector
		msevec[l]=merror;
		
		free_dvector(y_cvtrain, 1, n-n/k);
		free_dvector(y_cvtest, 1, n/k);
		free_dmatrix(x_cvtrain, 1, n-n/k,1,p);
		free_dmatrix(x_cvtest, 1, n/k,1,p);
		free_dvector(beta_mcp,1,p);
	}

	double s=0;
	for (i=1;i<=k;i++){
		s = s + msevec[i];
	}

	*m = s/(double)k;

	free_dvector(msevec,1,k);
	free_ivector(value,1,n);
}


//Function to give the grid between an interval
void grid(int n, double init, double end, double *g)
{
	int i;
	double range = fabs(init-end);
	g[1]=init;
	
	for (i=1;i<n;i++) g[i+1]=init - (range*i)/(double)n;

}

//Find the index of the minumum in a vector
void min_index(int n, double *x, int *idx)
{
	int i,index;
	double min;
	min = x[1];
	for (i=2;i<=n;i++){
		if (x[i]<min){
			min=x[i];
			index=i;
		}
	}
	
	*idx = index;
}

//Use k-fold cross-validation to choose the value of lambda and a
void cv_par(int k, int n, int p, double init_lam, double end_lam,
int iter_lam, double init_gamma, double end_gamma, int iter_gamma,
int max_iter, double tol, double *y_train, double **x_train, 
double *sel_lambda, double *sel_gamma)
{
	int i,j;
	
	//Get the grid of lambda
	double *lamb;
	lamb = dvector(1,iter_lam);
	grid(iter_lam,init_lam,end_lam,lamb);
	
	//Get the grid of a
	double *gam;
	gam = dvector(1, iter_gamma);
	grid(iter_gamma, init_gamma, end_gamma, gam);

	//Compute the MSE using k-fold cross-validation for each lambda
	double **mse_vector;
	mse_vector = dmatrix(1,iter_lam*iter_gamma,1,3);
	for (i=1;i<=iter_lam;i++){
		for (j=1;j<=iter_gamma;j++){
			double error=0;
			cv_mse(k,n,p,gam[j],lamb[i],max_iter,tol,y_train,x_train,&error);
			mse_vector[(i-1)*(iter_gamma)+j][1]=lamb[i];
			mse_vector[(i-1)*(iter_gamma)+j][2]=gam[j];
			mse_vector[(i-1)*(iter_gamma)+j][3]=error;
			if (gam[j]==1){
				printf("gamma is 1!");
				break;
			}
		}
	}

	int idx;
	double *mse_vec;
	mse_vec = dvector(1,iter_lam*iter_gamma);
	for (i=1;i<=iter_lam*iter_gamma;i++) mse_vec[i] = mse_vector[i][3];
	min_index(iter_lam*iter_gamma,mse_vec,&idx);
	
	//Compute the lambda and a that gives the lowest MSE
	*sel_lambda=mse_vector[idx][1];
	*sel_gamma = mse_vector[idx][2];
	
	free_dvector(lamb,1,iter_lam);
	free_dvector(gam,1,iter_gamma);
	free_dmatrix(mse_vector,1,(iter_lam)*(iter_gamma),1,3);
	free_dvector(mse_vec,1,(iter_lam)*(iter_gamma));
}

//MCP to do the variable selection via lambda chosen from 
//k-fold cross-validation
void var_selection(int k, int n, int p, double init_lam,
double end_lam, int iter_lam, double init_gamma, double end_gamma,
int iter_gamma, int max_iter, double tol, double *y_train,
double **x_train, int *sel_variable, double *beta_mcp)
{
	int i;
	double chosen_lambda, chosen_gamma;
	cv_par(k,n,p,init_lam,end_lam,iter_lam,init_gamma, end_gamma, iter_gamma,
		max_iter,tol,y_train, x_train, &chosen_lambda, &chosen_gamma);
	
	double *beta;
	beta = dvector(1,p);
	mcp(n, p, chosen_gamma,chosen_lambda,max_iter,tol, x_train, y_train, beta);
	for (i=1;i<=p;i++) beta_mcp[i]=beta[i];
	
	for (i=1;i<=p;i++){
		if (beta[i]!=0) sel_variable[i]=i;
		else sel_variable[i]=0;
	}
	
	free_dvector(beta,1,p);
}


//Push 0s to the end
void pushZerosToEnd(int n, int *x)
{
    int count = 1;  // Count of non-zero elements
	int i;
    // Traverse the array. If element encountered is non-
    // zero, then replace the element at index 'count' 
    // with this element
    for (i = 1; i <= n; i++)
        if (x[i] != 0)
            x[count++] = x[i]; // here count is 
                                   // incremented
 
    // Now all non-zero elements have been shifted to 
    // front and  'count' is set as index of first 0. 
    // Make all elements 0 from count to end.
    while (count <= n)
        x[count++] = 0;
}

//Push 0s to the end
void pushZerosToEnd_double(int n, double *x)
{
    int count = 1;  // Count of non-zero elements
	int i;
    // Traverse the array. If element encountered is non-
    // zero, then replace the element at index 'count' 
    // with this element
    for (i = 1; i <= n; i++)
        if (x[i] != 0)
            x[count++] = x[i]; // here count is 
                                   // incremented
 
    // Now all non-zero elements have been shifted to 
    // front and  'count' is set as index of first 0. 
    // Make all elements 0 from count to end.
    while (count <= n)
        x[count++] = 0;
}


