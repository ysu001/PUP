/*$Header: /data/petsun4/data1/src_solaris/librms/RCS/dgeigen.c,v 1.4 2011/08/30 03:43:19 avi Exp $*/
/*$Log: dgeigen.c,v $
 * Revision 1.4  2011/08/30  03:43:19  avi
 * trap ill-conditioned sig2
 *
 * Revision 1.3  2011/08/21  02:47:25  avi
 * complete rewrite now with correct algorithm
 *
 * Revision 1.2  2009/01/20  03:57:17  avi
 * cosmetic changes
 *
 * Revision 1.1  2009/01/20  01:01:07  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <librms.h>

void errm (char* program);

void dgeigen (double *sig1, double *sig2, double *lambda, double *E, int *pn) {
	double	*C, *Ct, *M, *W1, *W2, dmin;
	int	n, nn, i, j, n1;
	int	debug = 0;
	int	status = 0;

	n = *pn;
	nn = n*n;
	n1 = n + 1;
	if (!(C  = (double *) malloc (nn*sizeof (double)))) errm ("dgeiegn");
	if (!(Ct = (double *) malloc (nn*sizeof (double)))) errm ("dgeiegn");
	if (!(M  = (double *) malloc (nn*sizeof (double)))) errm ("dgeiegn");
	if (!(W1 = (double *) malloc (nn*sizeof (double)))) errm ("dgeiegn");
	if (!(W2 = (double *) malloc (nn*sizeof (double)))) errm ("dgeiegn");

/*************************************/
/* find sqrt([sig2]) and its inverse */
/*************************************/
	dmatcop_ (sig2, M, pn);
	if (debug) printf ("before diagonalizing sig2\n");
	deigen_  (M, W2, pn);						/* M now has lambda2 */
	if (debug) printf ("after  diagonalizing sig2\n");
	if (debug) printf ("condition number before inversion = %.5g\n", M[0]/M[nn - 1]);
	for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
		if (j == i) {
			if (M[i + n*j] > 0.0) {
				M[i + n*j] = 1./sqrt(M[i + n*j]);	/* M now has 1/sqrt(lambda2) */
			} else {
				M[i + n*j] = 1.e6*M[0];
				status++;
			}
		} else {
			M[i + n*j] = 0.0;
		}
	}
	if (status) printf ("dgeigen warning: ill-conditioned sig2\n");
	dmatmul_   (W2, M, C, pn);				/* variables named as in Strang p. 344 */
	dtranspos_ (C,  Ct, pn);

/*********************************************/
/* compute [Ct][Sig1][C] and leave in lambda */
/*********************************************/
	dmatmul_ (Ct, sig1, M, pn); dmatmul_ (M, C, lambda, pn);

/***********************************************/
/* factor [Ct][Sig1][C] into [W1][lambda][W1t] */
/***********************************************/
	if (debug) printf ("before diagonalizing sig1\n");
	deigen_ (lambda, W1, pn);
	if (debug) printf ("after  diagonalizing sig1\n");

/***********************************/
/* compute eigenvectors as [C][W1] */
/***********************************/
	dmatmul_ (C, W1, E, pn);

	free (C); free (Ct); free (M); free (W1); free (W2);
}
