/*$Header: /data/petsun4/data1/src_solaris/librms/RCS/dsvdcmp0.c,v 1.5 2013/01/20 05:51:39 avi Exp $*/
/*$Log: dsvdcmp0.c,v $
 * Revision 1.5  2013/01/20  05:51:39  avi
 * reconstitution documentation
 *
 * Revision 1.4  2012/11/26  07:18:42  avi
 * void dsvdsqrtinv ()
 *
 * Revision 1.3  2012/09/05  22:07:43  avi
 * dsvdinv();
 *
 * Revision 1.2  2011/03/05  02:32:48  avi
 * ndsvdcmp0 () now sorts eigenvalues in reliably decreasing order
 *
 * Revision 1.1  2011/02/14  02:05:26  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define SIGN(a,b) ((b) > 0 ? fabs(a) : -fabs(a))
#define IMIN(a,b) ((a) < (b) ? (a) : (b))
#define FMAX(a,b) ((a) > (b) ? (a) : (b))

void svd_errm () {
	fprintf (stderr, "memory allocation error in module dsvdcmp0\n");
	exit (-1);
}

double **svd_calloc_double2 (int n1, int n2) {
	int	i;
	double	**a;

	if (!(a = (double **) malloc (n1 * sizeof (double *)))) svd_errm ();
	if (!(a[0] = (double *) calloc (n1 * n2, sizeof (double)))) svd_errm ();
	for (i = 1; i < n1; i++) a[i] = a[0] + i*n2;
	return a;
}

void svd_free_double2 (double **a) {
	free (a[0]);
	free (a);
}

double pythag (double a, double b) {
	double absa, absb, q;
	
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) {
		q = absb/absa;
		return absa*sqrt(1.0+q*q);
	} else {
		q = absa/absb;
		return (absb == 0.0) ? 0.0 : absb*sqrt(1.0+q*q);
	}
}

/********************************************************************************/
/* version of svdcmp (NR p. 68) in which arrays a, w, and v use C indexing	*/
/********************************************************************************/
/* If dsvdcmp0 is used with FORTRAN compatible arrays (first index rows) then	*/
/* [E], [F], [G] are ncol x npts and V is npts x npts and W is npts		*/
/* tol is least allowable ratio of smallest to greatest eigenvalue		*/
/*	[F] <- [E];								*/					
/*	     dsvdcmp0 (E, ncol, npts, W, V);					*/
/*	n = ndsvdcmp0 (E, ncol, npts, W, V, tol);				*/
/*	original [E] can be reconstituted as [Vt]*[W]*[E]			*/
/*	for (j = 0; j < ncol; j++) for (i = 0; i < npts; i++) {			*/
/*		for (k = 0; k < n; k++) G[j][i] += V[i][k]*W[k]*E[j][k];	*/
/*	}									*/
/* now [G] = [F] = original [E]							*/
/********************************************************************************/
void dsvdcmp0 (double **a, int m, int n, double *w, double **v) {
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	if (!(rv1 = (double *) malloc (n*sizeof(double)))) svd_errm ();
	g = scale = anorm = 0.0;
	for (i = 0; i < n; i++) {
		l=i+1;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		if (i < m) {
			for (k = i; k < m; k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k = i; k < m; k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				a[i][i] = f-g;
				for (j = l; j < n; j++) {
					for (s = 0.0, k = i; k < m; k++) s += a[k][i]*a[k][j];
					f = s/h;
					for (k = i; k < m; k++) a[k][j] += f*a[k][i];
				}
				for (k = i; k < m; k++) a[k][i] *= scale;
			}
		}
		w[i] = scale*g;
		g = s = scale = 0.0;
		if (i < m && i != n - 1) {
			for (k = l; k < n; k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k = l; k < n; k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				a[i][l] = f-g;
				for (k = l; k < n; k++) rv1[k] = a[i][k]/h;
				for (j = l; j < m; j++) {
					for (s = 0.0, k = l; k < n; k++) s += a[j][k]*a[i][k];
					for (         k = l; k < n; k++) a[j][k] += s*rv1[k];
				}
				for (k = l; k < n; k++) a[i][k] *= scale;
			}
		}
		anorm = FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i = n - 1; i >= 0; i--) {
		if (i < n - 1) {
			if (g) {
				for (j = l; j < n; j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j = l; j < n; j++) {
					for (s = 0.0, k = l; k < n; k++) s += a[i][k]*v[k][j];
					for (         k = l; k < n; k++) v[k][j] += s*v[k][i];
				}
			}
			for (j = l; j < n; j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g = rv1[i];
		l=i;
	}
	for (i = IMIN(m - 1, n - 1); i >= 0; i--) {
		l = i+1;
		g = w[i];
		for (j = l; j < n; j++) a[i][j]=0.0;
		if (g) {
			g = 1.0/g;
			for (j = l; j < n; j++) {
				for (s = 0.0, k = l; k < m; k++) s += a[k][i]*a[k][j];
				f = (s/a[i][i])*g;
				for (k = i; k < m; k++) a[k][j] += f*a[k][i];
			}
			for (j = i; j < m; j++) a[j][i] *= g;
		}
		else for (j = i; j < m; j++) a[j][i] = 0.0;
		a[i][i] += 1.0;
	}
	for (k = n - 1; k >= 0; k--) {
		for (its = 0; its < 30; its++) {
			flag=1;
			for (l = k; l > 0; l--) {	/* correct apparent error in original code */
				nm = l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0; s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g = w[i];
					h = pythag(f,g);
					w[i] = h;
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					for (j = 0; j < m; j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y*c + z*s;
						a[j][i]  = z*c - y*s;
					}
				}
			}
			z=w[k];
			if (l==k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j = 0; j < n; j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its >= 30) { /* correct error in original code */
				printf ("dsvdcmp0: no convergence in 30 iterations\n");
			}
			x = w[l];
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c = s = 1.0;
			for (j=l;j<=nm;j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c + g*s;
				g = g*c - x*s;
				h = y*s;
				y *= c;
				for (jj = 0; jj < n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x*c + z*s;
					v[jj][i] = z*c - x*s;
				}
				z = pythag(f,h);
				w[j] = z;
				if (z) {
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g + s*y;
				x = c*y - s*g;
				for (jj = 0; jj < m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y*c + z*s;
					a[jj][i] = z*c - y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free (rv1);
}

int mdsvdcmp0 (double **a, int m, int n, double *w, double **v, double tol) { /* obsolete */
	int	i, ii, j, *lw;
	double	wmax = 0.0;

	if (!(lw = (int *) calloc (n, sizeof(int)))) svd_errm ();
	for (i = 0; i < n; i++) if (w[i] > wmax) wmax = w[i];
	for (i = 0; i < n; i++) if (w[i] > wmax*tol) lw[i]++;
	
	for (ii = i = 0; i < n; i++ ) {
		if (lw[i]) {
			if (ii != i) {
				for (j = 0; j < n; j++) v[j][ii] = v[j][i];
				for (j = 0; j < m; j++) a[j][ii] = a[j][i];
							   w[ii] =    w[i];
			}
			ii++;
		}
	}

	free (lw);
	return ii;
}

int ndsvdcmp0 (double **a, int m, int n, double *w, double **v, double tol) {
	int	i, ii, j, k;
	double	t, wmax = 0.0;
	
	for (i = 0; i < n; i++) {
		for (k = ii = i; k < n; k++) if (w[k] > w[ii]) ii = k;
		if (ii != i) {
			for (j = 0; j < n; j++) {
				t = v[j][ii]; v[j][ii] = v[j][i]; v[j][i] = t;
			}
			for (j = 0; j < m; j++) {
				t = a[j][ii]; a[j][ii] = a[j][i]; a[j][i] = t;
			}
			t = w[ii]; w[ii] = w[i]; w[i] = t;
		}
	}
	for (i = 1; i < n; i++) if (w[i] < w[0]*tol) break;
	return i;
}

void dsvdinv (double **A, int n, double *det) {
	double	**V, **H, *W;
	int	i, j, k;

	V = svd_calloc_double2 (n, n);
	H = svd_calloc_double2 (n, n);
	if (!(W = (double *) malloc (n*sizeof (double)))) svd_errm ();

	dsvdcmp0 (A, n, n, W, V);
	*det = 1.; for (i = 0; i < n; i++) *det *= W[i];
	for (j = 0; j < n; j++) for (i = 0; i < n; i++) H[j][i] = A[j][i]/W[i];
	for (j = 0; j < n; j++) for (i = 0; i < n; i++) {
		for (A[i][j] = k = 0; k < n; k++) A[i][j] += V[i][k]*H[j][k];
	}

	svd_free_double2 (V); svd_free_double2 (H); free (W);
}

void dsvdsqrtinv (double **A, int n, double *det) {
	double	**V, **H, *W;
	int	i, j, k;

	V = svd_calloc_double2 (n, n);
	H = svd_calloc_double2 (n, n);
	if (!(W = (double *) malloc (n*sizeof (double)))) svd_errm ();

	dsvdcmp0 (A, n, n, W, V);
	*det = 1.; for (i = 0; i < n; i++) *det *= W[i];
	for (j = 0; j < n; j++) for (i = 0; i < n; i++) H[j][i] = A[j][i]/sqrt(W[i]);
	for (j = 0; j < n; j++) for (i = 0; i < n; i++) {
		for (A[i][j] = k = 0; k < n; k++) A[i][j] += V[i][k]*H[j][k];
	}

	svd_free_double2 (V); svd_free_double2 (H); free (W);
}
