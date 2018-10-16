/*$Header: /data/petsun4/data1/src_solaris/sqrt_4dfp/RCS/statsub.c,v 1.2 2008/03/14 03:48:22 avi Exp $*/
/*$Log: statsub.c,v $
 * Revision 1.2  2008/03/14  03:48:22  avi
 * linux compliant
 *
 * Revision 1.1  2005/09/10  00:18:33  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#ifndef _FEATURES_H		/* defined by math.h in linux */
	#include <ieeefp.h>	/* needed in Solaris but not linux */
#endif
#include <string.h>

#define LN_2SQRTPI	1.26551212348464539647
#define LN2PI		1.83787706640934548355

/**********************************/
/* externals in bertulani_betai.c */
/**********************************/
double betai	(double a, double b, double x);
double beta	(double a, double b);
double lbeta	(double a, double b);

double qerfc (double z) {
	double	oz2, a, b, q;
	int	i, j, k;

	oz2 = 1./(z*z);
	q = a = j = k = 1; if (oz2 > DBL_MIN) for (i = 1; i < 11; i++) {
		k *= -1;
		b = j*oz2;
		a *= b;
		j *= (2*i + 1);
		if (0) printf ("i=%d j=%d k=%d a=%.5e b=%.6e oz2=%.6e\n", i, j, k, a, b, oz2);
		q += k*a;
		if (a < DBL_EPSILON || b > 0.5) break;
	};
	return q;
}

double lnp_z (double z) {
	double	x;

	x = z/M_SQRT2;
	if (z < 38.) {
		return -M_LN2 + log(erfc(x));
	} else {
		return log(qerfc(z)) - (x*x + log(x) + LN_2SQRTPI);
	}
}

double p_z (double z) {
	double	x;

	x = z/M_SQRT2;
	if (z < 38.) {
		return 0.5*erfc(x);
	} else {
		return 0.5*exp(-x*x)*qerfc(z)/(sqrt(M_PI)*x);
	}
}

double dlnpdz (double z) {
	double	x;

	x = z/M_SQRT2;
	if (z < 38.) {
		return -sqrt(2.0/M_PI)*exp(-x*x)/erfc(x);
	} else {
		return -z/qerfc(z);
	}
}

double z_lnp (double lnp) {
	double	q, z, d;
	int	i;

	if (lnp > -M_LN2) {
		fprintf (stderr, "ln(p) > -M_LN2\n");
		exit (-1);
	}
	q = -2.0*lnp - LN2PI;
	z = (q > 0) ? sqrt(q) : 0.1;
	i = 0; do {
		d = (lnp_z(z) - lnp)/dlnpdz (z);
		z -= d;
		if (i++ > 5) {
			fprintf (stderr, "z_lnp iter limit (5) reached\n");
			break;
		}
	} while (fabs(d) > 1.e-6);
	return z;
}

/*********************************************************/
/* this algorithm has 3 digit accuracy for t >  5 (D=10) */
/*                and 6 digit accuracy for t > 20 (D=10) */
/* 0.5*betai(0.5*D, 0.5, D/(D + t*t)) is always accurate */
/* except as t->Inf and D >> 1, when infinitely small	 */
/* betai values are not different from 0		 */
/*********************************************************/
double p_t (double t, double D) {
	double	q, a, b, ot2;
	int	i, j, k;

	ot2 = 1./(t*t);
	q = a = j = k = 1; if (ot2 > DBL_MIN) for (i = 1; i < 11; i++) {
		k *= -1;
		b = j*(D/(D + 2*i))*ot2;
		a *= b;
		if (0) printf ("i=%d j=%d k=%d a=%f b=%f\n", i, j, k, a, b);
		j *= (2*i + 1);
		q += k*a;
		if (a < DBL_EPSILON || b > 0.5) break;
	}
	return q*pow(1.0 + t*t/D, -0.5*(D-1.0)) / (t*sqrt(D)*beta(0.5*D, 0.5));
}

double lnp_t (double t, double D) {
	double	q, a, b, ot2;
	int	i, j, k;

	q = -M_LN2 + log(betai(0.5*D, 0.5, D/(D + t*t)));
	if (finite(q)) return q;
	if (0) printf("using Woolrich\n");
	ot2 = 1./(t*t);
	q = a = j = k = 1; if (ot2 > DBL_MIN) for (i = 1; i < 11; i++) {
		k *= -1;
		b = j*(D/(D + 2*i))*ot2;
		a *= b;
		if (0) printf ("i=%d j=%d k=%d a=%f b=%f\n", i, j, k, a, b);
		j *= (2*i + 1);
		q += k*a;
		if (a < DBL_EPSILON || b > 0.5) break;
	}
	return log(q) - 0.5*(D-1.0)*log(1.0 + t*t/D) - log(t) - 0.5*log(D) - lbeta(0.5*D, 0.5);
}

/*********************************************************/
/* this algorithm is not accurate for low f		 */
/* betai(0.5*D2, 0.5*D1, D2/(D2 + D1*f)) is always	 */
/* accurate except as f->Inf and returned betai values	 */
/* are not different from 0				 */
/*********************************************************/
double p_f (double f, double D1, double D2) {
	double	q, a, b, r;
	int	i, k;

	r = D2/(D1*f);
	q = a = k = 1; for (i = 1; i < 32; i++) {
		k *= -1;
		b = r*(2*i - D1)/(2*i + D2);
		a *= b;
		q += k*a;
		if (0) printf ("i=%d k=%d a=%f b=%f q=%f\n", i, k, a, b, q);
		if (fabs(a) < DBL_EPSILON || (i > 20 && fabs(b) > 0.5)) break;
	}
	if (0) printf("p_f i=%d\n", i);
	return 2.0*q*pow(D1/D2, 0.5*D1) / (D1*beta(0.5*D2, 0.5*D1)
			*pow(1. + 1./r, 0.5*(D1 + D2) - 1.0)*pow(f, 1.0 - 0.5*D1));
}

double lnp_f (double f, double D1, double D2) {
	double	q, a, b, r;
	int	i, k;

	q = log(betai(0.5*D2, 0.5*D1, D2/(D2 + D1*f)));
	if (finite(q)) return q;
	if (0) printf("using Woolrich\n");
	r = D2/(D1*f);
	q = a = k = 1; for (i = 1; i < 32; i++) {
		k *= -1;
		b = r*(2*i - D1)/(2*i + D2);
		a *= b;
		q += k*a;
		if (0) printf ("i=%d k=%d a=%f b=%f q=%f\n", i, k, a, b, q);
		if (fabs(a) < DBL_EPSILON || (i > 20 && fabs(b) > 0.5)) break;
	}
	if (0) printf("p_f i=%d\n", i);
	return M_LN2 + log(q) + 0.5*D1*log(D1/D2) - log(D1) - lbeta(0.5*D2, 0.5*D1)
			- (0.5*(D1 + D2) - 1.0)*log(1. + 1./r) - (1.0 - 0.5*D1)*log(f);
}

/* int main (int argc, char *argv[]) { */
void	statsub_test () {
	double	p, z, x, D, D1, D2, t, f;
	double	a, b;
	int	i, j, k;

	D1 = 41.; D2 = 33.;
	t = 1.; for (i = 1; i <= 100; i++) {
		z = (double) i; f = t;
		x = z/M_SQRT2;
if (1)		printf ("%15.6f%15.6f%15.6f\n",	z, log(0.5*erfc(z/M_SQRT2)), lnp_z(z));
if (0)		printf ("%10d%15.6f%15.6f\n",	i, lnp_z(z), 0.5*sqrt(2.0/M_PI)*exp(-x*x)/p_z(z));
if (0)		printf ("%15.6f%15.6f%15.6f\n",	z, lnp_z(z), dlnpdz(z));
if (0)		printf ("%15.6f%15.6f%15.6f\n",	z, lnp_z(z), z_lnp(lnp_z(z)));
if (0)		printf ("%10d%15.6f%15.6e\n",	i, lnp_z(z), p_z(z));
if (0)		printf ("%15.6f%15.6e%15.6e\n",	t, p_t (t, D), 0.5*betai(0.5*D, 0.5, D/(D + t*t)));
if (0)		printf ("%15.6e%15.6e%15.6e\n",	f, p_f (f, D1, D2), 0.5*betai(0.5*D2, 0.5*D1, D2/(D2 + D1*f)));
if (0)		printf ("%15.6e%15.6e%15.6e\n",	f, lnp_f (f, D1, D2), log(betai(0.5*D2, 0.5*D1, D2/(D2 + D1*f))));
if (0)		printf ("%15.6e%15.6e%15.6e\n",	t, log(0.5*betai(0.5*D, 0.5, D/(D + t*t))), lnp_t(t, D));
		z *= 2.0;
		t *= 2.0;
	}
}
