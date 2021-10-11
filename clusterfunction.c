/*Brandon Woosuk Park
George Mason University*/

//This files corresponds with the Table 2 in section 6.4 in the manuscript.

/*Functions for clustering methods*/


//Find the maximum value from a p x p matrix
void max_matrix(int p , double **x, double *max)
{
	int i,j;
	double maximum;
	maximum=x[1][1];
	for (i=1;i<=p;i++){
		for (j=1;j<=p;j++){
			if (x[i][j]>maximum) maximum = x[i][j];
		}
	}
	*max =maximum;
}

//Find the maximum value from p x 1 vector
void max_vector(int p, double *x, double *max)
{
	int i;
	double maximum;
	maximum=x[1];
	for (i=2;i<=p;i++){
		if (x[i]>maximum) maximum = x[i];
	}	
	*max = maximum;
}

//Find the maximum value location from p x 1 vector
void max_vec_loc(int p, double *x, int *idx)
{
	int i, location;
	double maximum;
	maximum=x[1];
	location=1;
	for (i=2;i<=p;i++){
		if (x[i]>maximum){
			maximum = x[i];
			location = i;
		}
	}
	*idx = location;
}
//Find the nonzero minimum value from p x 1 vector
void min_vector(int p, double *x, double *min)
{
	int i;
	double minimum;
	for (i=1;i<=p;i++){
		if (x[i]!=0){
			minimum=x[i];	
		}
		else continue;
	}

	for (i=1;i<=p;i++){
		if (x[i]==0) continue;
		else{
		if (x[i]<minimum){
			minimum = x[i];
		}
		}
	}
	*min = minimum;
}
//Find the minimum value (except 0) location from p x 1 vector
void min_vec_loc(int p, double *x, int *idx)
{
	int i, location;
	double minimum;
	for (i=1;i<=p;i++){
		if (x[i]!=0){
			minimum=x[i];
			location=i;		
		}
		else continue;
	}

	for (i=1;i<=p;i++){
		if (x[i]==0) continue;
		else{
		if (x[i]<minimum){
			minimum = x[i];
			location = i;
		}
		}
	}
	*idx = location;
}


//Remove zeroes
void remove_zero(int p, int *x, int *new_p)
{
	int i,j,k;
	for (i = 1; i <= p; i++) {
      for (j = i + 1; j <= p;) {
         if (x[j] == x[i] || x[j]==0) {
            for (k = j; k <= p; k++) {
               x[k] = x[k + 1];
            }
            p--;
         } else
            j++;
      }
   }
   *new_p=p;
}
//Sort a vector
void sort_vec(int n, int *number)
{

     /*Sort the given array number , of length n*/     
	int i,j,a;
     for (i = 1; i <= n; ++i)
    {
        for (j = i + 1; j <= n; ++j)
        {
            if (number[i] > number[j])
            {
                a =  number[i];
                number[i] = number[j];
                number[j] = a;
            }
        }
    }
}



//Detect based on the test for the difference 
//Compute the average vector
void avg_vec(int n, int p, double **x, double *x_avgvec)
{
	int i,j;
	double *d;
	d = dvector(1,n);
	for (i=1;i<=p;i++){
		double avg;
		for (j=1;j<=n;j++){
			d[j]=x[j][i];
		}
		average(n,d,&avg);
		x_avgvec[i]=avg;
	}
	free_dvector(d,1,n);
}

//Test the difference based on the mean and variance estimated from bootstrap
void test(int n, double x,double y,double var_x, double var_y, double var_xy, double *t)
{
double top, bot;
top = x-y;
bot = var_x + var_y - 2*var_xy;
if (bot==0) *t=0;
else *t = top/pow(bot,.5); 
}

//Sequential Testing for original
void diff_test_boot(int n,int p,int c,double alpha,
double *degree, double **var, int **cluster)
{
	int i,j,df;
	df = n+n-2;
	double tscore=qt(alpha/2,df,1,0);	
	double *temp;
	temp = dvector(1,p);
	for (i=1;i<=p;i++) temp[i]=degree[i];
	for (i=1;i<=c;i++){
		cluster[i][1]=i;
		double *h;
		int *g,max_idx;
		g = ivector(1,p);
		h = dvector(1,p);
		max_vec_loc(p,temp,&max_idx);
		
		for (j=1;j<=p;j++){
			double t;
			test(n,temp[max_idx],temp[j],var[max_idx][max_idx],var[j][j], var[max_idx][j],&t);
			if (fabs(t)<fabs(tscore)){
				g[j]=j;
				h[j]=temp[j];
			}
			else{
				g[j]=0;
				h[j]=0;
			}
		}
		
		pushZerosToEnd(p,g);
		int new_p;
		remove_zero(p,g,&new_p);
		sort_vec(new_p,g);

		
		for (j=1;j<=new_p;j++) cluster[i][j+1] = g[j];
		for (j=new_p+1;j<=p;j++) cluster[i][j+1]=0;
		
		for (j=1;j<=p;j++) g[j]=0;
		for (j=1;j<=p;j++) temp[j]=temp[j]-h[j];
		
		free_ivector(g,1,p);
		free_dvector(h,1,p);
	}
	free_dvector(temp,1,p);
}

//Sequential Testing for new defined metric
void diff_test_boot_new(int n,int p,int c,double alpha,
double *degree, double **var, int **cluster)
{
	int i,j,df;
	df = n+n-2;
	double tscore=qt(alpha/2,df,1,0);	
	double *temp;
	temp = dvector(1,p);
	for (i=1;i<=p;i++) temp[i]=degree[i];
	for (i=1;i<=c;i++){
		cluster[i][1]=i;
		double *h;
		int *g,min_idx;
		g = ivector(1,p);
		h = dvector(1,p);
		min_vec_loc(p,temp,&min_idx);
		
		for (j=1;j<=p;j++){
			double t;
			test(n,temp[min_idx],temp[j],var[min_idx][min_idx],var[j][j], var[min_idx][j],&t);
			if (fabs(t)<fabs(tscore)){
				g[j]=j;
				h[j]=temp[j];
			}
			else{
				g[j]=0;
				h[j]=0;
			}
		}
		
		pushZerosToEnd(p,g);
		int new_p;
		remove_zero(p,g,&new_p);
		sort_vec(new_p,g);
		
		for (j=1;j<=new_p;j++) cluster[i][j+1] = g[j];
		for (j=new_p+1;j<=p;j++) cluster[i][j+1]=0;
		
		for (j=1;j<=p;j++) g[j]=0;
		for (j=1;j<=p;j++) temp[j]=temp[j]-h[j];
		
		free_ivector(g,1,p);
		free_dvector(h,1,p);
	}
	free_dvector(temp,1,p);
}


