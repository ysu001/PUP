/*$Header: /home/usr/shimonyj/JSSutil/RCS/JSSstatistics.c,v 1.3 2008/12/31 01:42:25 shimonyj Exp $*/
/*$Log: JSSstatistics.c,v $
 * Revision 1.3  2008/12/31  01:42:25  shimonyj
 * Improve error output in gser
 *
 * Revision 1.2  2008/11/11  00:03:07  avi
 * correct rcs header
 *
 * Revision 1.1  2008/11/11  00:00:45  avi
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* used by betacf, gser, gcf */
#define MAXIT	100	
#define EPS 	3.0e-7
#define FPMIN	1.0e-30
#define ITMAX	100


/* fit to a straight line, NR page 665 added linear cor coeff */
/* linear cor coeff is not the same as the Pearson R coefficient */
void fit(float x[], float y[], int ndata, float sig[], int mwt, float *a,
	float *b, float *siga, float *sigb, float *chi2, float *q, float *lcc)
{
	float gammq(float a, float x);
	int i;
	float wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
	
	*b=0;
	if (mwt) {
		ss=0.0;
		for (i=1;i<=ndata;i++) {
			wt=1.0/(sig[i]*sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	}
	else {
		for (i=1;i<=ndata;i++) {
			sx += x[i];
			sy += y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt) {
		for (i=1;i<=ndata;i++) {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
		}
	}
	else {
		for (i=1;i<=ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
	}
	*b /= st2;
	*a = (sy-sx*(*b))/ss;
	*siga = sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb = sqrt(1.0/st2);
	*lcc = -sx/(ss*st2*(*siga)*(*sigb));
	*chi2 = 0.0;
	if (mwt==0) {
		for (i=1;i<=ndata;i++)
			*chi2 += (y[i]-(*a)-(*b)*x[i])*(y[i]-(*a)-(*b)*x[i]);
		*q=1.0;
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *= sigdat;
		*sigb *= sigdat;
	}
	else {
		for (i=1;i<=ndata;i++) 
			*chi2 += ((y[i]-(*a)-(*b)*x[i])/sig[i])*((y[i]-(*a)-(*b)*x[i])/sig[i]);
			*q = gammq(0.5*(ndata-2),0.5*(*chi2)); 
	}
}


/* f-test using NR p.619 */
void ftest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *f, float *prob)
{
	void avevar(float data[],unsigned long n,float *ave,float *var);
	float betai(float a,float b,float x);
	float var1,var2,ave1,ave2,df1,df2;
	
	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	if (var1 > var2) {
		*f=var1/var2;
		df1=n1-1;
		df2=n2-1;
	}
	else {
		*f=var2/var1;
		df1=n2-1;
		df2=n1-1;
	}
	*prob = 2.0*betai(0.5*df2,0.5*df1,df2/(df2+df1*(*f)));
	if (*prob > 1.0) *prob=2.0-*prob;
}

/* return ave and var of distribution NR p. 617 */
void avevar(float data[], unsigned long n, float *ave, float *var)
{
	unsigned long j;
	float s,ep;
	
	for (*ave=0.0, j=1; j<=n; j++) *ave += data[j];
	*ave /= n;
	*var=ep=0.0;
	for (j=1; j<=n; j++) {
		s=data[j]-(*ave);
		ep += s;
		*var += s*s;
	}
	*var = (*var - ep*ep/n)/(n-1);
}

/* returns log of the gamma function NR p. 214 */
float gammln(float xx)
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


/* returns the incomplete beta function NR p. 227 */
float betai(float a, float b, float x)
{
	float betacf(float a, float b, float x);
	float gammln(float xx);
	float bt;
	
	if (x<0.0 || x>1.0) {
		printf("Error(betai):bad input x\n");
		return(0);
	}
	if (x == 0.0 || x == 1.0) bt=0.0;
	else bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0)) return(bt*betacf(a,b,x)/a);
	else return(1.0-bt*betacf(b,a,1.0-x)/b);
}

/* evaluates continued fraction for betai NR p. 227 */
float betacf(float a, float b, float x)
{
	int m,m2;
	float aa,c,d,del,h,qab,qam,qap;
	
	qab=a+b; qap=a+1.0; qam=a-1.0;
	c=1.0; d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d; h=d;
	
	for (m=1;m<=MAXIT;m++) {
		m2 = 2*m;
		aa = m*(b-m)*x/((qam+m2)*(a+m2));
		d = 1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c = 1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d = 1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c = 1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m>MAXIT) {
		printf("Error(betacf):a or b too big, or MAXIT too small\n");
		return(0.0);
	}
	return(h);
}

/* paired t-test NR p. 618 small prob indicates significant difference in means */
void tptest(float data1[],float data2[],unsigned long n,float *t,float *prob)
{
	void avevar(float data[],unsigned long n,float *ave,float *var);
	float betai(float a,float b,float x);
	unsigned long j;
	float var1,var2,ave1,ave2,sd,df,cov=0.0;
	
	avevar(data1,n,&ave1,&var1);
	avevar(data2,n,&ave2,&var2);
	
	for (j=1;j<=n; j++)
		cov += (data1[j]-ave1)*(data2[j]-ave2);
	cov /= df = n-1;
	sd = sqrt((var1+var2-2.0*cov)/n);
	*t = (ave1-ave2)/sd;
	*prob = betai(0.5*df,0.5,df/(df+(*t)*(*t)));
}

/* chi square statistic NR p. 621 */
void chsone(float bins[],float ebins[],int nbins,int knstrn,float *df,float *chsq,float *prob)
{
	float gammq(float a, float x);
	int j;
	float temp;
	
	*df = nbins - knstrn;
	*chsq = 0.0;
	for (j=1; j<=nbins; j++) {
		if (ebins[j] <= 0.0) {
			printf("Error(chsone):bad expected number\n");
			return;
		}
		temp = bins[j] - ebins[j];
		*chsq += temp*temp/ebins[j];
	}
	*prob = gammq(0.5*(*df),0.5*(*chsq));
}

void chstwo(float bins1[],float bins2[],int nbins,int knstrn,float *df,float *chsq,float *prob)
{
	float gammq(float a, float x);
	int j;
	float temp;
	
	*df = nbins - knstrn;
	*chsq = 0.0;
	for (j=1; j<=nbins; j++) {
		if (bins1[j] == 0.0 && bins2[j] == 0.0) --(*df);
		else {
			temp = bins1[j] - bins2[j];
			*chsq += temp*temp/(bins1[j]+bins2[j]);
		}
	}
	*prob = gammq(0.5*(*df),0.5*(*chsq));
}

/* return the incomplete gamma function P and Q NR p. 218 */
float gammp(float a, float x)
{
	void gcf(float *gammcf,float a,float x,float *gln);
	void gser(float *gamser,float a, float x,float *gln);
	float gamser,gammcf,gln;
	
	if (x<0.0 || a <= 0.0) {
		printf("Error(gammp):invalid arguments a=%f x=%f\n",a,x);
		return(0);
	}
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return(gamser);
	}
	else {
		gcf(&gammcf,a,x,&gln);
		return(1.0-gammcf);
	}
}

float gammq(float a, float x)
{
	void gcf(float *gammcf,float a,float x,float *gln);
	void gser(float *gamser,float a, float x,float *gln);
	float gamser,gammcf,gln;
	
	if (x<0.0 || a <= 0.0) {
		printf("Error(gammq):invalid arguments\n");
		return(0);
	}
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return(1.0-gamser);
	}
	else {
		gcf(&gammcf,a,x,&gln);
		return(gammcf);
	}
}

void gser(float *gamser, float a, float x, float *gln)
{
	float gammln(float xx);
	int n;
	float sum,del,ap;
	
	*gln = gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) printf("Error(gser):x less than 0\n");
		*gamser = 0.0;
		return;
	}
	else {
		ap = a;
		del = sum = 1.0/a;
		for (n=1; n<=ITMAX; n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser = sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		printf("Error(gser):a=%f too large, x=%f, ITMAX too small\n",a,x);
		return;
	}
}

void gcf(float *gammcf, float a, float x, float *gln)
{
	float gammln(float xx);
	int i;
	float an,b,c,d,del,h;
	
	*gln = gammln(a);
	b = x + 1.0 - a;
	c = 1.0/FPMIN;
	d = 1.0/b;
	h = d;
	for (i=1; i<=ITMAX; i++) {
		an = -i*(i-a);
		b += 2.0;
		d = an*d + b;
		if (fabs(d) < FPMIN) d = FPMIN;
		c = b + an/c;
		if (fabs(c) < FPMIN) c = FPMIN;
		d = 1.0/d;
		del = d*c;
		h *= del;
		if (fabs(del - 1.0) < EPS) break;
	}
	if (i > ITMAX) printf("Error(gcf):a too large, ITMAX too small\n");
	*gammcf = exp(-x + a*log(x)-(*gln))*h;
}

/* error function routines, NR pg. 220 */
float erfunc(float x)
{
	float gammp(float a, float x);
	
	if (x<0.0) return(-gammp(0.5,x*x));
	else return(gammp(0.5,x*x));
}

float erffc(float x)
{
	float gammp(float a, float x);
	float gammq(float a, float x);

	if (x<0.0) return(1.0+gammp(0.5,x*x));
	else return(gammq(0.5,x*x));
}

float erfcc(float x)
{
	float t,z,ans;
	
	z = fabs(x);
	t = 1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));

	if (x>=0.0) return(ans);
	else return(2.0-ans);
}
