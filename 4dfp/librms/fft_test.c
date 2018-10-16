/*$Header: /data/petsun4/data1/src_solaris/librms/RCS/fft_test.c,v 1.1 2006/02/19 01:28:22 avi Exp avi $*/
/*$Log: fft_test.c,v $
 * Revision 1.1  2006/02/19  01:28:22  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#define	NMAX	200

void errm (char* program) {
	fprintf (stderr, "%s: memory allocation error\n", program);
	exit (-1);
}

main (int argc, char *argv[]) {
/*void find_bad_n () {*/
	int		c, i, j, k;
/*******/
/* FFT */
/*******/
	float		a[NMAX], b[NMAX];
	int		one = 1, negone = -1, m, n, zero = 0;

	printf ("n="); fflush (stdout);
	for (n = 4; n < NMAX; n++) {
		m = n;
/*
		m = npad_ (&n, &zero);
		if (m >= NMAX) continue;
		if (n ==  99) continue;
		if (n == 195) continue;
		if (n == 291) continue;
		if (n == 297) continue;
		if (n == 325) continue;
		if (n == 363) continue;
		if (n == 387) continue;
		if (n == 390) continue;
		if (n == 392) continue;
		if (n == 483) continue;
		if (n == 485) continue;
		if (n == 495) continue;
testing stopped here */
		for (i = 0; i < m; i++) a[i] = b[i] = 0.;
		a[0] = 1.;
		printf (" %d %d %.2f\n", n, m, a[0]); fflush (stdout);
		fft_ (a, b, &one, &m, &one, &negone);
		/*realt_ (a, b, &one, &no2, &one, &negone);*/
	}
	printf ("\n");  fflush (stdout);
}

void transform_gaussian () {
#define N	64
	int	n = N;
	int	one = 1, negone = -1;
	int	i, ii, j, k;
	float	a[N], b[N], c[N], e[N], o[N];
	float	sigma = 5.;
	double	q, r, t, z;

	for (i = - n/2; i < n/2; i++) {
		ii = (i < 0) ? i + n : i;
		z = i - 0.5;
		c[ii] = exp(-0.5*z*z/(sigma*sigma));;
	}
	if (0) for (k = n, i = 1; i < n/2; i++) {
		c[--k] = -c[i];
		if (1) {c[i] *=  i; c[k] *= -i;}
	}
	for (i = 0; i < n; i++) {a[i] = c[i]; b[i] = 0.0;}
	for (i = 0; i < n; i++) printf ("%10d%10.4f%10.4f\n", i, a[i], b[i]);
printf ("OK\n");
	fft_ (a, b, &one, &n, &one, &negone);
	for (i = 0; i <= n/2; i++) {
		r = sqrt (a[i]*a[i] + b[i]*b[i]);
		t = atan2 (-b[i], a[i]);
		printf ("%10d%10.4f%10.4f%10.4f%10.4f\n", i, a[i], b[i], r, t*n/(2.*M_PI));
	}
}

void precision_test () {
	extern double dnormal ();
#define N	64
	int	n = N;
	int	one = 1, negone = -1;
	int	i, j, k;
	float	a[N], b[N], c[N], d[N];
	
	for (i = 0; i < n; i++) {
		a[i] = c[i] = dnormal ();
		b[i] = d[i] = dnormal ();
	}
	fft_ (a, b, &one, &n, &one, &negone);
	fft_ (a, b, &one, &n, &one, &one);
	for (i = 0; i < n; i++) {
		printf ("%10d%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n",
		i, a[i], b[i], c[i], d[i], c[i]-a[i], d[i]-b[i]);
	}

	exit (0);
}
