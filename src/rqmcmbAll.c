#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>

#define MAXN 200000
#define MAXP 100
#define XACCF .001
#define XACCD .001

#define MAXIT 100
#define BMAXIT 100

#define UNUSED (-1.11e30)
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

/*Numerical Recipes Random Number generator*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0(long *idum){
  long k;
  float ans;
  
  *idum ^= MASK;
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  ans=AM*(*idum);
  *idum ^= MASK;
  return ans;
}



int allZero;

/****************************************************************/
/*  Various Utilities off the Web (for use with sort2())        */ 
/****************************************************************/

unsigned long *lvector(long nl, long nh)
/* allocate a long int vector with subscript range v[nl..nh] */
{
    long *v = (long *)malloc((size_t)((nh-nl+2) * sizeof(long)));
    if (v == NULL) error("allocation failure in lvector()");
    return (v-nl+1);
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free a long int vector allocated with lvector() */
{
    free(v+nl-1);
}


/****************************************************************/
/*  Sort Function                                               */ 
/****************************************************************/

void sort2(unsigned long n, double arr[], double brr[]){
  unsigned long i,ir=n,j,k,l=1,*istack;
  int jstack=0;
  double a,b,temp;
  istack=lvector(1,NSTACK);

  for (;;) { 
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	b=brr[j];
	for (i=j-1;i>=l;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;
      }
      if (!jstack) {
	free_lvector(istack,1,NSTACK);
	return;
      }
      ir=istack[jstack]; 
      l=istack[jstack-1];
      jstack -= 2;
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
      SWAP(brr[k],brr[l+1])
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	SWAP(brr[l],brr[ir])
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	SWAP(brr[l+1],brr[ir])
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
	SWAP(brr[l],brr[l+1])
      }
      i=l+1; 
      j=ir;
      a=arr[l+1]; 
      b=brr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break; 
	SWAP(arr[i],arr[j]) 
	SWAP(brr[i],brr[j])
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      brr[l+1]=brr[j];
      brr[j]=b;
      jstack += 2;
      if (jstack > NSTACK) printf("NSTACK too small in sort2.\n");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}


/****************************************************************/
/*  Multiply two vectors x and c, each of length pp             */ 
/****************************************************************/

double mprodx(double *x, double *c, int pp){
  int i;
  double sum;
  
  sum=0;
  for(i=0;i<pp;i++){
    sum=sum+x[i]*c[i];}
  
  return sum;
}

/****************************************************************/
/*  Sign Function                                               */ 
/****************************************************************/

double sign(double x){
  double sign;

  if(x>0.0)
    sign=1;
  else if (x<0.0)
    sign=-1;
  else 
    sign=0;	
  
  return sign;
}


/****************************************************************/
/* This function finds a weighted taustar quantile, instead of  */
/* using bisection to solve the MCMB equations. See new paper.  */
/****************************************************************/

double func(double *x, double *y,  double tau, double *tTilda, double *A,  double 
sum_right, double sumxij, double sumabsxij, int j, int pp, int nn){
  int i,m;
  double xj[MAXN], yj[MAXN], z[MAXN], wt[MAXN],wtsum;
  unsigned long mm;
  double taustar, pwtsum, ans, large;
  
  for(i=0;i<nn;i++){
    yj[i]=y[i];
    xj[i]=x[i*pp+j];
  }
  xj[nn]=-sum_right/tau;
  yj[nn]=10e16;
  wtsum=sumabsxij + fabs(xj[nn]);
  
  /*update the first n elements*/

  /*sort function sorts arr[1...n], so set the first elements of z and wt to 0, so that the actual values in z and wt start from position z[1] and wt[1]*/
  m=1;
  z[0]=0; wt[0]=0; 
  for(i=0;i<nn;i++){
    if(fabs(xj[i]) > 10e-16){
      z[m]=(y[i]-mprodx(&x[i*pp],tTilda,pp)+tTilda[j]*xj[i])/xj[i];
      wt[m]=fabs(xj[i])/wtsum;
      m=m+1;
      /*printf("i=%d xj[i]=%f m=%d z[m]=%f wt[m]=%f\n",i,xj[i],m,z[m],wt[m]);*/}
    else{printf("fabs(xj[i])<10e-16\n");}
  }
  z[m]=10e16*sign(xj[nn]);
  large=z[m];
  wt[m]=fabs(xj[nn])/wtsum;
  
  /*calculate taustar*/
  taustar=(tau-0.5)*(sumxij+xj[nn])/(wtsum)+0.5;
  
  if(m==0){
    printf("Error: one design variable contains all 0s.\n");
    allZero=1;}
  
  mm=m;
  sort2(m, z, wt);

  pwtsum=0; i=1;
  ans=z[1];
  while(pwtsum<=taustar && i<nn+1){
    pwtsum=pwtsum+wt[i];
    ans=z[i];
    i++;
  }

  /*resample if pick the largest z*/ 
  if(fabs(ans) > 10e15){
    printf("Picked infinity; need to resample\n");
    return 1.0;
  }
  
  /*printf("q=%f sumabsxij=%f  fabs(xj[nn])=%f taustar=%f\n", ans, sumabsxij, fabs(xj[nn]), taustar);*/

  return ans;

}

/****************************************************************/
/* RQMCMB                                                       */
/****************************************************************/


void rqmcmb(double *x, double *y,  double *tau, double *theta_tilda, double *A, double *zstar, 
	    double *sumxij, double *sumabsxij, long *n, long *p, int *success, double *theta, long *MAXK, long *seed){
  
  int i, j, jj, k, nn, pp;
  double sum, s[MAXP],tau2, tTilda[MAXP];
  long t1;
  int rand_ind;
  extern int allZero;
  
  pp=(int) *p;
  nn=(int) *n;
  tau2=(double) *tau;
  allZero=0;
  
  Rprintf("RQMCMB is running.\n");
  
  for(i=0;i<pp;i++){
    theta[i]=theta_tilda[i];
    tTilda[i]=theta_tilda[i];
  }
  
  /*time(&t1);
    srand(seed[0]);*/
  success[0]=1;  
  
  for(k=0;k<*MAXK; k++){
    
    /*bootstrap from x[i][j] individually for each j, return sum */
    for(j=0;j<pp;j++){
      sum=0; 
      for(i=0;i<nn;i++){
        rand_ind=(int)(nn*j+nn*ran0(seed));
	sum=sum+zstar[rand_ind];
      }
      s[j]=sqrt(nn)/sqrt(nn-pp)* sum;
    }
        
    for(j=0;j<pp;j++){
      theta[(k+1)*pp+j]=func(x, y, tau2, tTilda, A, s[j], sumxij[j], sumabsxij[j], j, pp, nn);
      if(allZero==1){
	success[0]=0;
	return;}
      if(theta[(k+1)*pp+j]==1.0){
        for(jj=0;jj<pp;jj++){
          tTilda[jj]=theta[k*pp+jj];}
        k=k-1;
        break;
      }
      else { 
	tTilda[j]=theta[(k+1)*pp+j];
      }
    }
  }
  
}












