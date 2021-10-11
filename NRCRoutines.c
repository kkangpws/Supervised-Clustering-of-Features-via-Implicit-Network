/*Routines used in C*/

#define NR_END 1
#define FREE_ARG char*
int ludcmp_flag, flag;

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
        fprintf(stderr,"Numerical Recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
        if (!v) nrerror("allocation failure in vector()");
        return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
        if (!v) nrerror("allocation failure in ivector()");
        return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
        unsigned char *v;

        v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
        if (!v) nrerror("allocation failure in cvector()");
        return v-nl+NR_END;
}
unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
        unsigned long *v;

        v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
        if (!v) nrerror("allocation failure in lvector()");
        return v-nl+NR_END;
}       

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        if (!v) nrerror("allocation failure in dvector()");
        return v-nl+NR_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;
        
        /* allocate pointers to rows */
        m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;
        
        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}
int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int **m;
        
        /* allocate pointers to rows */
        m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;


        /* allocate rows and set pointers to them */
        m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;
        
        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
        
        /* return pointer to array of pointers to rows */
        return m;
}
        
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
        long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
        long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
        float **m;
        
        /* allocate array of pointers to rows */
        m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
        if (!m) nrerror("allocation failure in submatrix()");
        m += NR_END;
        m -= newrl;
 
        /* set pointers to rows */
        for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
        
        /* return pointer to array of pointers to rows */
        return m;
}
        
double **subdmatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
        long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
        long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
        double **m;

        /* allocate array of pointers to rows */
        m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
        if (!m) nrerror("allocation failure in submatrix()");
        m += NR_END;
        m -= newrl;

        /* set pointers to rows */
        for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
        long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        float **m;

        /* allocate pointers to rows */
        m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
        if (!m) nrerror("allocation failure in convert_matrix()");
        m += NR_END;
        m -= nrl;

        /* set pointers to rows */
        m[nrl]=a-ncl;
        for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
        /* return pointer to array of pointers to rows */
        return m; 
}       
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
        long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
        float ***t;

        /* allocate pointers to pointers to rows */
        t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
        if (!t) nrerror("allocation failure 1 in f3tensor()");
        t += NR_END;
        t -= nrl;
        
        /* allocate pointers to rows and set pointers to them */
        t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
        if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
        t[nrl] += NR_END;
        t[nrl] -= ncl;
 
        /* allocate rows and set pointers to them */
        t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
        if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
        t[nrl][ncl] += NR_END;
        t[nrl][ncl] -= ndl;
        
        for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
        for(i=nrl+1;i<=nrh;i++) {
                t[i]=t[i-1]+ncol;
                t[i][ncl]=t[i-1][ncl]+ncol*ndep;
                for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
        }

        /* return pointer to array of pointers to rows */
        return t;
}       
void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}
        
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}       
        
void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}
        
void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}       
       
void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{       
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}
        
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}
 
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
        free((FREE_ARG) (b+nrl-NR_END));
}
        
void free_subdmatrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
        free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
        free((FREE_ARG) (b+nrl-NR_END));
}
 
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{       
        free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
        free((FREE_ARG) (t[nrl]+ncl-NR_END));
        free((FREE_ARG) (t+nrl-NR_END));
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5.){2p491&.#@#Q0JS[. */

void lubksb(double **a, int n, int *indx, double *b)
{
        int i,ii=0,ip,j;
        double sum;

        for (i=1;i<=n;i++) {
                ip=indx[i];
                sum=b[ip];
                b[ip]=b[i];
                if (ii)
                        for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
                else if (sum) ii=i;
                b[i]=sum;
        }
        for (i=n;i>=1;i--) {
                sum=b[i];
                for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
                b[i]=sum/a[i][i];
        }
}
void ludcmp(double **a, int n, int *indx, double *d)
{
        int i,imax,j,k;
        double big,dum,sum,temp;
        double *vv, TINY = 1.0e-20;

        vv=dvector(1,n);
        *d=1.0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if ((temp=fabs(a[i][j])) > big) big=temp;
                if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
                vv[i]=1.0/big;
        }
        for (j=1;j<=n;j++) {
                for (i=1;i<j;i++) {
                        sum=a[i][j];
                        for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                }
                big=0.0;
                for (i=j;i<=n;i++) {
                        sum=a[i][j];
                        for (k=1;k<j;k++)
                                sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                        if ( (dum=vv[i]*fabs(sum)) >= big) {
                                big=dum;
                                imax=i;
                        }
                }
                if (j != imax) {
                        for (k=1;k<=n;k++) {
                                dum=a[imax][k];
                                a[imax][k]=a[j][k];
                                a[j][k]=dum;
                        }
                        *d = -(*d);
                        vv[imax]=vv[j];
                }
                indx[j]=imax;
                if (a[j][j] == 0.0) a[j][j]=TINY;
                if (j != n) {
                        dum=1.0/(a[j][j]);
                        for (i=j+1;i<=n;i++) a[i][j] *= dum;
                }
        }
        free_dvector(vv,1,n);
}

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

double ran1(idum)
long *idum;
{
        static long ix1,ix2,ix3;
        static double r[98];
        double temp;
        static int iff=0;
        int j;
        void nrerror();

        if (*idum < 0 || iff == 0) {
                iff=1;
                ix1=(IC1-(*idum)) % M1;
                ix1=(IA1*ix1+IC1) % M1;
                ix2=ix1 % M2;
                ix1=(IA1*ix1+IC1) % M1;
                ix3=ix1 % M3;
                for (j=1;j<=97;j++) {
                        ix1=(IA1*ix1+IC1) % M1;
                        ix2=(IA2*ix2+IC2) % M2;
                        r[j]=(ix1+ix2*RM2)*RM1;
                }
                *idum=1;
        }
        ix1=(IA1*ix1+IC1) % M1;
        ix2=(IA2*ix2+IC2) % M2;
        ix3=(IA3*ix3+IC3) % M3;
        j=1 + ((97*ix3)/M3);
        if (j > 97 || j < 1) nrerror("RAN1: This cannot happen.");
        temp=r[j];
        r[j]=(ix1+ix2*RM2)*RM1;
        return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

double gasdev(idum)
long *idum;
{
        static int iset=0;
        static double gset;
        double fac,r,v1,v2;
        double ran2();

        if  (iset == 0) {
                do {
                        v1=2.0*ran2(idum)-1.0;
                        v2=2.0*ran2(idum)-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

double gamdev(int ia, long *idum)
{
        double ran2(long *idum);
        void nrerror(char error_text[]);
        int j;
        double am,e,s,v1,v2,x,y;

        if (ia < 1) nrerror("Error in routine gamdev");
        if (ia < 6) {
                x=1.0;
                for (j=1;j<=ia;j++) x *= ran2(idum);
                x = -log(x);
        } else {
                do {
                        do {
                                do {
                                        v1=2.0*ran2(idum)-1.0;
                                        v2=2.0*ran2(idum)-1.0;
                                } while (v1*v1+v2*v2 > 1.0);
                                y=v2/v1;
                                am=ia-1;
                                s=sqrt(2.0*am+1.0);
                                x=s*y+am;
                        } while (x <= 0.0);
                        e=(1.0+y*y)*exp(am*log(x/am)-s*y);
                } while (ran2(idum) > e);
        }
        return x;
}

void tdev(int df, int dim, double *trannum, long *idum)
{
	int i,j;
	double Y;

	Y = gamdev(df/2.0,idum)*2; Y = sqrt(Y/df);
	for (i=1; i<=dim; i++)
		trannum[i] = gasdev(idum)/Y;

}

/*********************************************************************
 *   Returns the value ln(Gamma(xx)) for xx>0. Full accuracy is obtained
 *     for xx > 1. For 0<xx<1, the reflection formula can be used first:
 *
 *          Gamma(1-z) = pi/Gamma(z)/sin(pi*z) = pi*z/Gamma(1+z)/sin(pi*z)
 *          *********************************************************************/
double gammln(xx)
double xx;
{
        double x,tmp,ser;
        static double cof[6]={76.18009173,-86.50532033,24.01409822,
                -1.231739516,0.120858003e-2,-0.536382e-5};
        int j;

        x=xx-1.0;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.0;
        for (j=0;j<=5;j++) {
                x += 1.0;
                ser += cof[j]/x;
        }
        return -tmp+log(2.50662827465*ser);
}

#define PI 3.141592654

double poidev(xm,idum)
double xm;
long *idum;
{
        static double sq,alxm,g,oldm=(-1.0);
        double em,t,y;
        double ran2(),gammln();

        if (xm < 12.0) {
                if (xm != oldm) {
                        oldm=xm;
                        g=exp(-xm);
                }
                em = -1;
                t=1.0;
                do {
                        em += 1.0;
                        t *= ran2(idum);
                } while (t > g);
        } else {
                if (xm != oldm) {
                        oldm=xm;
                        sq=sqrt(2.0*xm);
                        alxm=log(xm);
                        g=xm*alxm-gammln(xm+1.0);
                }
                do {
                        do {
                                y=tan(PI*ran2(idum));
                                em=sq*y+xm;
                        } while (em < 0.0);
                        em=floor(em);
                        t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
                } while (ran2(idum) > t);
        }
        return em;
}

#undef PI


#define IM1 2147483563
#define IM2 2147483399   
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32  
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(idum)
long *idum;
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        double temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        } 
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV   
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software 5.){2p491&.#@#Q0JS[. */

void sort(n,ra)
int n;
double *ra;
{
        int l,j,ir,i;
        double rra;

        l=(n >> 1)+1;
        ir=n;
        for (;;) {
                if (l > 1)
                        rra=ra[--l];
                else {
                        rra=ra[ir];
                        ra[ir]=ra[1];
                        if (--ir == 1) {
                                ra[1]=rra;
                                return;
                        }
                }
                i=l;
                j=l << 1;
                while (j <= ir) {
                        if (j < ir && ra[j] < ra[j+1]) ++j;
                        if (rra < ra[j]) {
                                ra[i]=ra[j];
                                j += (i=j);
                        }
                        else j=ir+1;
                }
                ra[i]=rra;
        }
}

/*
double gammln(double xx)
{
        double x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
                0.1208650973866179e-2,-0.5395239384953e-5};
        int j;

        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
}
*/
double fgamma(double x)
{
        int i,j;
        double a1,a2,pi=3.1415926535;

        if (x - (int)(x) == 0)
        {

                a1 = 1; a2 = 1;
                i = (int)(x);
                j = 1;
                while (j<i)
                {
                        a2 = a1*j;
                        a1 = a2;
                        j++;
                }
        }
        else if (2*x - (int)(2*x) == 0)
        {

                a1 = 1; a2 = 1;
                i = (int)(x);
                j = 1;
                while (j<=i)
                {
                        a2 = a1*(j-0.5);
                        a1 = a2;
                        j++;
                }
                a2 *= pow(pi,0.5);
        }
        else a2 = exp(gammln(x));

        return a2;
}

double erfcc(double x)
{
        double t,z,ans;

        z=fabs(x);
        t=1.0/(1.0+0.5*z);
        ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
                t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
                t*(-0.82215223+t*0.17087277)))))))));
        return x >= 0.0 ? ans : 2.0-ans;
}

double Phi(double x)
{
	return 1 - 0.5*erfcc(x/sqrt(2.0));
}
// N(0,1) density function
double fnorm(double x)
{
	double pi = 3.141592654;

	return exp(-x*x/2.0)/sqrt(2*pi);
}

double dvector_percentile(int n, double *a, double percent)
{
        double percentile, diff,*tempa;
        int i,n_percent;
	
	if (n == 1) return(a[1]);
	else 
	{
        tempa = dvector(1,n);

        for (i=1;i<=n;i++)
                tempa[i] = a[i];

        sort(n,tempa);

        n_percent = (int)(n*percent);

        diff = n*percent - n_percent;
        if (diff>0)
                percentile = tempa[n_percent+1];
        else percentile = (tempa[n_percent]+tempa[n_percent+1])/2;

        free_dvector(tempa,1,n);

        return(percentile);
	}
}

void dmatrix_inv(double **A, double **invA, int n)
{
        double d, *col,**tempA;
        int i,j,*indx;

        col=dvector(1,n);
        indx=ivector(1,n);
        tempA=dmatrix(1,n,1,n);

        for (i=1;i<=n ;i++ )
                for (j=1;j<=n;j++)
                        tempA[i][j]=A[i][j];

        ludcmp(tempA,n,indx,&d);
        for (j=1;j<=n;j++)
        {
                for (i=1;i<=n;i++) col[i]=0.0;
                col[j]=1.0;
                lubksb(tempA,n,indx,col);
                for (i=1;i<=n;i++) invA[i][j]=col[i];
        }

        free_dvector(col,1,n);
        free_ivector(indx,1,n);
        free_dmatrix(tempA,1,n,1,n);

        return;
}

void choldc(double **a, int n, double *p)
/*
        given a positive-definite symmetric dmatrix a[1..n][1..n]. this routine constructs Cholesky
        decomposition. A=L * L'. On input, only the upper triangle of a need be given; it is not
        modified. The Cholesky factor L is returned in the lower triangle of a, except for its
        diagonal elements which are return in p[1..n]
*/
{
        int i,j,k,l,m;
        double sum;
        char temp;

        for (i=1;i<=n;i++) {
                for (j=i;j<=n;j++) {
                        for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
                        if (i == j) 
			{
                                if (sum <= 0.0)
                                {
                                        printf("choldc failed\n");
	                                flag = 1;
                                }
                        
                        p[i]=sqrt(sum);
                        } else a[j][i]=sum/p[i];
                }
        }
}

void dmatrix_invsqrt(double **a, int n, double **b)
/*
        b = a^{-1/2}
        a = L*L` (cholesky decomposition) and then b = L^{-1}
*/
{
        int i,j;
        double **tempa,*diag;

        tempa = dmatrix(1,n,1,n); diag = dvector(1,n);

        for (i=1; i<=n; i++)
                for (j=1; j<=n; j++)
                        tempa[i][j] = a[i][j];

        choldc(tempa,n,diag);

        for (i=1; i<=n; i++)
        {
                tempa[i][i] = diag[i];

                for (j=1; j<i; j++)
                        tempa[j][i] = 0;

        }

        dmatrix_inv(tempa,b,n);
        free_dmatrix(tempa,1,n,1,n); free_dvector(diag,1,n);
}

void dmatrix_transpose(double **a, double **b, int m, int n)
{
        int i,j;
        
        for (i=1;i<=m;i++)
                for (j=1;j<=n;j++)
                b[j][i]=a[i][j];  
                
        return; 
}

double dmatrix_determinant(double **a, int n)
{
        int i,j, *indx;
        double **tempa, d;

        if (n > 1) {
        indx = ivector(1,n);
        tempa = dmatrix(1,n,1,n);

        for (i=1; i<=n; i++)
                for (j=1; j<=n; j++)
                        tempa[i][j] = a[i][j];
        ludcmp(tempa,n,indx,&d);
        for (j=1; j<=n; j++) d *= tempa[j][j];

        free_ivector(indx,1,n);
        free_dmatrix(tempa,1,n,1,n);
        }
        else d = a[1][1];

        return(d);
}

void solve_linear_equation(double **a, double *b, int n, double *x)
/*
        Solve linear algebra equation Ax = b
*/
{
        int i,j,*indx;
        double **tempa,d;

        indx = ivector(1,n);
        tempa = dmatrix(1,n,1,n);
        for (i=1; i<=n; i++) x[i] = b[i];
        for (i=1; i<=n; i++)
                for (j=1; j<=n; j++)
                        tempa[i][j] = a[i][j];

        ludcmp(tempa,n,indx,&d);
        lubksb(tempa,n,indx,x);

        free_ivector(indx,1,n);
        free_dmatrix(tempa,1,n,1,n);
}

// MATRIX ROUTINES
//      return X'
void dmatrix_tX(int m, int n, double **X, double **tX)
{
        int i, j;

        for (i=1; i<=m; i++)
                for (j=1; j<=n; j++)
                        tX[j][i] = X[i][j];
}
//      return A'A
void dmatrix_tXX(int m, int n, double **A, double **tAA)
{
        int i, j, k;

        for (j=1; j<=n; j++)
                for (k=1; k<=j; k++)
        {
                tAA[j][k] = 0;
                for (i=1; i<=m; i++)
                        tAA[j][k] += A[i][j]*A[i][k];
                tAA[k][j] = tAA[j][k];
        }
}
//      return X'Y, Y is a vector
void dmatrix_tXY(int m, int n, double **X, double *Y, double *tXY)
{
        int i, j, k;

        for (j=1; j<=n; j++)
        {
                tXY[j] = 0;
                for (i=1; i<=m; i++)
                        tXY[j] += X[i][j]*Y[i];
        }
}
//      return Z'X
void dmatrix_tZX(int m, int n, int k, double **Z, double **X, double **tZX)
{
        int i, j, l;

        for (j=1; j<=n; j++)
                for (l=1; l<=k; l++)
        {
                tZX[j][l] = 0;
                for (i=1; i<=m; i++)
                        tZX[j][l] += Z[i][j]*X[i][l];
        }
}
//      return Z'HZ
void dmatrix_tZHZ(int m, int n, double **Z, double **H, double **tZHZ)
{
        int i, j, k, l;

        for (k=1; k<=n; k++)
                for (l=1; l<=k; l++)
        {
                tZHZ[k][l] = 0;
                for (i=1; i<=m; i++)
                        for (j=1; j<=m; j++)
                                tZHZ[k][l] += Z[i][k]*H[i][j]*Z[j][l];
                tZHZ[l][k] = tZHZ[k][l];
        }
}
//      return ZHZ': Z(m X n), H(n X n)
void dmatrix_ZHtZ(int m, int n, double **Z, double **H, double **ZHtZ)
{
    int i, j, k, l;
    
    for (k=1; k<=m; k++)
        for (l=1; l<=k; l++)
        {
            ZHtZ[k][l] = 0;
            for (i=1; i<=n; i++)
                for (j=1; j<=n; j++)
                    ZHtZ[k][l] += Z[k][i]*H[i][j]*Z[l][j];
            ZHtZ[l][k] = ZHtZ[k][l];
        }
}
//      return AB
void dmatrix_AB(int m, int n, int k, double **A, double **B, double **AB)
{
        int i, j, l;

        for (i=1; i<=m; i++)
                for (l=1; l<=k; l++)
        {
                AB[i][l] = 0;
                for (j=1; j<=n; j++)
                        AB[i][l] += A[i][j]*B[j][l];
        }
}
//      return Ab
void dmatrix_Ab(int m, int n, double **A, double *b, double *Ab)
{
        int i, j;

        for (i=1; i<=m; i++)
        {
                Ab[i] = 0;
                for (j=1; j<=n; j++)
                        Ab[i] += A[i][j]*b[j];
        }
}

void meanvar(int n, double *x, double *meanx, double *varx)
{
	int i, l;
	double m, v;
// incorporate missing values
	m = 0; v = 0; l = 0;
	for (i=1; i<=n; i++)
		if (x[i] != -99)
		 { l++; m += x[i];} 
	m /= l;
	for (i=1; i<=n; i++) 
		if (x[i] != -99)
		v += pow(x[i]-m,2.0);
	v /= l-1;

	*meanx = m;
	*varx = v;
}

void print_matrix(int m, int n, double **X)
{
	int i, j;
    
	for (i=1; i<=m; i++)
	{
		for (j=1; j<=n; j++)
			printf("%f ", X[i][j]);
		printf("\n");
	}
}
/*
double rbinom(double p, long *idum)
{
	double u;
	
	u = ran2(idum);
	if (u <= p) return 1;
	else return 0;
}
*/
//generate Gamma(alpha, 1.0)
/*
double rgamma(double alpha, long *idum)
{
    int i, j;
    double n[3], b1[3], b2[3], c1[3], c2[3], v1, v2, w1[3], w2[3], y[3], x[3];
    double ran2();
    
    for (i=1; i<=1; i++)
    {
        if (alpha<= 0) printf("bad alpha in Gammdev");
        else if (alpha <= 0.4) n[i] = 1.0/alpha;
        else if (alpha <=4.0) n[i] =1.0/alpha *(1+(alpha-0.4)/3.6);
        else n[i] =1/sqrt(alpha);
        
        b1[i] = alpha - 1.0/n[i];
        b2[i] = alpha + 1.0/n[i];
        
        if (alpha <= 0.4) c1[i] = 0;
        else c1[i] = b1[i]*(log(b1[i])-1)/2.0;
        
        c2[i] = b2[i]*(log(b2[i])-1)/2.0;
    }
    
    do {
        do {
            v1 = ran2(idum);
            v2 = ran2(idum);
            
            for (i=1; i<=1; i++) {
                w1[i] = c1[i]+log(v1);
                w2[i] = c2[i]+log(v2);
                y[i] = n[i]*(b1[i]*w2[i]-b2[i]*w1[i]); }
        } while (y[1] <0);
        
        for (i=1; i<=1; i++) x[i] = n[i]*(w2[i]-w1[i]);
    } while (log(y[1])<x[1]);
    
    return exp(x[1]);
}
*/


double pythag(double a, double b)
{
        double absa,absb;
        absa=fabs(a);
        absb=fabs(b);
        if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
        else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void svdcmp(double **a, int m, int n, double w[], double **v)
{
        double pythag(double a, double b);
        int flag,i,its,j,jj,k,l,nm;
        double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

        rv1=dvector(1,n);
        g=scale=anorm=0.0;
        for (i=1;i<=n;i++) {
                l=i+1;
                rv1[i]=scale*g;
                g=s=scale=0.0;
                if (i <= m) {
                        for (k=i;k<=m;k++) scale += fabs(a[k][i]);
                        if (scale) {
                                for (k=i;k<=m;k++) {
                                        a[k][i] /= scale;
                                        s += a[k][i]*a[k][i];
                                }
                                f=a[i][i];
                                g = -SIGN(sqrt(s),f);
                                h=f*g-s;
                                a[i][i]=f-g;
                                for (j=l;j<=n;j++) {
                                        for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
                                        f=s/h;
                                        for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                                }
                                for (k=i;k<=m;k++) a[k][i] *= scale;
                        }
                }
                w[i]=scale *g;
                g=s=scale=0.0;
                if (i <= m && i != n) {
                        for (k=l;k<=n;k++) scale += fabs(a[i][k]);
                        if (scale) {
                                for (k=l;k<=n;k++) {
                                        a[i][k] /= scale;
                                        s += a[i][k]*a[i][k];
                                }
                                f=a[i][l];
                                g = -SIGN(sqrt(s),f);
                                h=f*g-s;
                                a[i][l]=f-g;
                                for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
                                for (j=l;j<=m;j++) {
                                        for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
                                        for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
                                }
                                for (k=l;k<=n;k++) a[i][k] *= scale;
                        }
                }
                anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
        }
        for (i=n;i>=1;i--) {
                if (i < n) {
                        if (g) {
                                for (j=l;j<=n;j++)
                                        v[j][i]=(a[i][j]/a[i][l])/g;
                                for (j=l;j<=n;j++) {
                                        for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
                                        for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
                                }
                        }
                        for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
                }
                v[i][i]=1.0;
                g=rv1[i];
                l=i;
        }
        for (i=IMIN(m,n);i>=1;i--) {
                l=i+1;
                g=w[i];
                for (j=l;j<=n;j++) a[i][j]=0.0;
                if (g) {
                        g=1.0/g;
                        for (j=l;j<=n;j++) {
                                for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
                                f=(s/a[i][i])*g;
                                for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                        }
                        for (j=i;j<=m;j++) a[j][i] *= g;
                } else for (j=i;j<=m;j++) a[j][i]=0.0;
                ++a[i][i];
        }
        for (k=n;k>=1;k--) {
                for (its=1;its<=30;its++) {
                        flag=1;
                        for (l=k;l>=1;l--) {
                                nm=l-1;
                                if ((double)(fabs(rv1[l])+anorm) == anorm) {
                                        flag=0;
                                        break;
                                }
                                if ((double)(fabs(w[nm])+anorm) == anorm) break;
                        }
                        if (flag) {
                                c=0.0;
                                s=1.0;
                                for (i=l;i<=k;i++) {
                                        f=s*rv1[i];
                                        rv1[i]=c*rv1[i];
                                        if ((double)(fabs(f)+anorm) == anorm) break;
                                        g=w[i];
                                        h=pythag(f,g);
                                        w[i]=h;
                                        h=1.0/h;
                                        c=g*h;
                                        s = -f*h;
                                        for (j=1;j<=m;j++) {
                                                y=a[j][nm];
                                                z=a[j][i];
                                                a[j][nm]=y*c+z*s;
                                                a[j][i]=z*c-y*s;
                                        }
                                }
                        }
                        z=w[k];
                        if (l == k) {
                                if (z < 0.0) {
                                        w[k] = -z;
                                        for (j=1;j<=n;j++) v[j][k] = -v[j][k];
                                }
                                break;
                        }
                        if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
                        x=w[l];
                        nm=k-1;
                        y=w[nm];
                        g=rv1[nm];
                        h=rv1[k];
                        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
                        g=pythag(f,1.0);
                        f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                        c=s=1.0;
                        for (j=l;j<=nm;j++) {
                                i=j+1;
                                g=rv1[i];
                                y=w[i];
                                h=s*g;
                                g=c*g;
                                z=pythag(f,h);
                                rv1[j]=z;
                                c=f/z;
                                s=h/z;
                                f=x*c+g*s;
                                g = g*c-x*s;
                                h=y*s;
                                y *= c;
                                for (jj=1;jj<=n;jj++) {
                                        x=v[jj][j];
                                        z=v[jj][i];
                                        v[jj][j]=x*c+z*s;
                                        v[jj][i]=z*c-x*s;
                                }
                                z=pythag(f,h);
                                w[j]=z;
                                if (z) {
                                        z=1.0/z;
                                        c=f*z;
                                        s=h*z;
                                }
                                f=c*g+s*y;
                                x=c*y-s*g;
                                for (jj=1;jj<=m;jj++) {
                                        y=a[jj][j];
                                        z=a[jj][i];
                                        a[jj][j]=y*c+z*s;
                                        a[jj][i]=z*c-y*s;
                                }
                        }
                        rv1[l]=0.0;
                        rv1[k]=f;
                        w[k]=x;
                }
        }
        free_dvector(rv1,1,n);
}

// dmaxtrix_zero
void dmatrix_zero(int m, int n, double **A)
{
	int i, j;
    
	for (i=1; i<=m; i++)
		for (j=1; j<=n; j++)
			A[i][j] = 0;
}
// dvector_zero
void dvector_zero(int n, double *a)
{
	int i;
	
	for (i=1; i<=n; i++)
		a[i] = 0;
}
void meancov(int n, int p, double **X, double *meanX, double **covX)
{
	int i, j, k;
    
	for (j=1; j<=p; j++)
	{
		meanX[j] = 0;
		for (k=1; k<=p; k++)
			covX[j][k] = 0;
	}
    
	for (j=1; j<=p; j++)
	{
		for (i=1; i<=n; i++)
			meanX[j] += X[i][j];
		meanX[j] /= (double)n;
	}
	for (j=1; j<=p; j++)
	{
		for (k=1; k<=j; k++)
			for (i=1; i<=n; i++)
                covX[j][k] += (X[i][j]-meanX[j])*(X[i][k]-meanX[k])/(n-1.0);
	}
	for (j=1; j<=p; j++)
		for (k=1; k<j; k++)
			covX[k][j] = covX[j][k];
}

