/*$Header*/
/*$Log: lud.c,v $
 * Revision 1.2  2007/02/13  07:12:22  avi
 * code cleaned
 **/
/* 2002.10.01 - LU decomposition R. Sachs */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define TINY	1.0e-20;

typedef float type;

void ludcmp (type **a, int n, int *indx, type *d) {
	int	i, imax, j, k;
	type	big, dum, sum, temp;
	type	*vv;

	vv = (type *) calloc (n, sizeof(type));
	if (!vv) {
		fprintf (stderr, "calloc fail ludcmp\n");
		exit (-1);
	}
	*d = 1.0;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n ; j++) if ((temp = fabs(a[i][j])) > big) big = temp;
		if (big == 0.0) { 
			fprintf (stderr, "ludcmp singular matrix i, j = %d %d\n", i, j);
			exit (-1);
                } 
		vv[i] = 1.0/big;
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i]*fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = TINY;
		if (j != n) {
			dum = 1.0/(a[j][j]);
			for (i = j + 1; i < n; i++) a[i][j] *= dum;
		}
	}
        free (vv);
}

void lubksb (type **a, int n, int *indx, type b[]) {
	int	i, ii = 0, ip, j;
	type	sum;

	for (i = 0; i < n; i++) {	/* forward substitution */
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		for (j = 0; j <= i - 1; j++) sum -= a[i][j]*b[j];
		b[i] = sum;
	}
	for (i = n - 1; i >= 0; i--) {	/* back substitution */
		sum = b[i];
		for (j = i + 1; j < n; j++) sum -= a[i][j]*b[j];
		b[i] = sum/a[i][i];
	}
}

type norm (type **a, int n) {
	int	i, j;
	type	sum, rmx;

	for (rmx = 0.0, i = 0; i < n; i++ ) {
		for (sum = 0.0, j = 0; j < n; j++) sum += fabs (a[i][j]);
		if (sum > rmx) rmx = sum;
	}
	return rmx;
}
