/*$Header: /data/petsun4/data1/src_solaris/sqrt_4dfp/RCS/bertulani_betai.c,v 1.3 2006/09/25 00:33:07 avi Exp $*/
/*$Log: bertulani_betai.c,v $
 * Revision 1.3  2006/09/25  00:33:07  avi
 * #include <stdio.h>
 *
 * Revision 1.2  2005/09/19  06:09:47  avi
 * correct typo in comment field
 *
 * Revision 1.1  2005/09/10  00:18:23  avi
 * Initial revision
 **/
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

double beta (double a, double b) {
	return exp (lgamma(a) + lgamma(b) - lgamma(a + b));
}

double lbeta (double a, double b) {
	return lgamma(a) + lgamma(b) - lgamma(a + b);
}

double dbetaidx (double a, double  b, double x) {
	if (x <= 0.0 || x >= 1.0) {
		fprintf (stderr, "dbetaidx argument out of range\n");
		exit (-1);
	}
	return pow(x, a - 1.) * pow(1. - x, b - 1.) / beta(a, b);
}

void betai_test () {
/*	main (int argc, char *argv[]) {	*/
	double betai (double, double, double);
	double x, a, b, q, p;
	int i;
	
	a = 10; b = 0.5;
	printf ("#1/beta(%f,%f)=%12.6e\n", a, b, 1./beta(a,b));
	for (i = 1; i < 100; i++) {
		x = (double) i / 100.;
		q = betai(a, b, x);
		printf ("%15.10f%15.6e%15.6e%10.4f\n", x, betai(a, b, x), dbetaidx(a, b, x), log(q));
	}
}

/****************************************************************/
/* Returns the incomplete beta function, I_x(a,b)		*/
/* I_x(a,b) = (int_0^x t^(a-1) * (1-t)^(b-1)  dt)/B(a,b)	*/
/*	where B(a,b) is the beta function			*/
/* B(a,b) = int_0^1 t^(a-1) * (1-t)^(b-1) dt			*/
/* 	= Gamma(a) * Gamma (b) / Gamma(a+b)			*/
/* C.A. Bertulani May/16/2000					*/
/****************************************************************/
double betai (double a, double b, double x) {
	double betacf (double a, double b, double x);
	double bt;

	if (x < 0.0 || x > 1.0) {
		fprintf (stderr, "betai argument out of range\n");
		exit (-1);
	}
	if (x == 0.0) return 0.0;
	if (x == 1.0) return 1.0;
/**********************************************/
/* factors in front of the continued fraction */
/**********************************************/
	bt = exp (lgamma(a + b) - lgamma(a) - lgamma(b) + a*log(x) + b*log(1.0 - x));
	if (x < (a+1.0)/(a+b+2.0)) {	
/***********************************/
/* use continued fraction directly */
/***********************************/
		return bt*betacf(a,b,x)/a;
	} else {
/******************************************************************/
/* use continued faction after making the symmetry transformation */
/******************************************************************/
		return 1.0-bt*betacf(b,a,1.0-x)/b;
	}
}

#define MAXIT 100
/***************************************************************************************/
/* evaluate continued fraction for incomplete beta function by modified Lentz's method */
/***************************************************************************************/
double betacf (double a, double b, double x) {
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;	
	qam=a-1.0;
/********************************/
/* first step of Lentz's method */
/********************************/
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < DBL_MIN) d=DBL_MIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
/************************/
/* even recurrence step */
/************************/
		d=1.0+aa*d;
		if (fabs(d) < DBL_MIN) d=DBL_MIN;
		c=1.0+aa/c;
		if (fabs(c) < DBL_MIN) c=DBL_MIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
/***********************/
/* odd recurrence step */
/***********************/
		d=1.0+aa*d;
		if (fabs(d) < DBL_MIN) d=DBL_MIN;
		c=1.0+aa/c;
		if (fabs(c) < DBL_MIN) c=DBL_MIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < DBL_EPSILON) break;
	}
	if (m > MAXIT) {
		fprintf (stderr, "a or b too big, or MAXIT too small in betacf");
	}
	if (0) printf ("betacf iters = %d\n", m);
	return h;
}
