/*$Header: /data/petsun4/data1/src_solaris/normalize_4dfp/RCS/img2hist.c,v 1.3 2006/09/29 17:08:09 avi Exp $*/
/*$Log: img2hist.c,v $
 * Revision 1.3  2006/09/29  17:08:09  avi
 * Solaris 10
 *
 * Revision 1.2  1999/03/24  00:26:57  avi
 * range tayloring using marg[]
 *
 * Revision 1.1  1996/07/15  05:42:50  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <endianio.h>

#define CUM_D	0.01

static char rcsid[] = "$Id: img2hist.c,v 1.3 2006/09/29 17:08:09 avi Exp $";
void img2hist (float *imgt, int nvox, float *minval, float *maxval, int *hist, int nbin) {
	double		range;
	int		i, j, k;
	int		*marg;
	int		debug = 0;

	printf ("%s\n", rcsid);
	*minval = FLT_MAX; *maxval = -FLT_MAX;
	for (j = 0; j < nbin; j++) hist[j] = 0;
	if (!(marg = (int *) calloc (nbin, sizeof (int)))) errm ("img2hist");

	for (i = 0; i < nvox; i++) {
		if (imgt[i] > *maxval) *maxval = imgt[i];
		if (imgt[i] < *minval) *minval = imgt[i];
	}
	range = *maxval - *minval;
	if (debug) printf ("img2hist: min=%10.4f\tmax=%10.4f\trange=%10.4f\n", *minval, *maxval, range);

	for (i = 0; i < nvox; i++) {
		k = (nbin * (imgt[i] - *minval)) / range;
		if (k < nbin && k > 0) marg[k]++;
	}
	for (k = 2; k < nbin; k++) marg[k] += marg[k - 1];

	if (debug) for (k = 0; k < nbin; k++) {
		printf ("%6d%10.6f\n", k, (float) marg[k] / (float) marg[nbin - 1]);
	}
	for (i = 1; i < nbin; i++) if (((float) marg[i] / (float) marg[nbin - 1]) > CUM_D) break;
	for (j = 1; j < nbin; j++) if (((float) marg[j] / (float) marg[nbin - 1]) > 1.0-CUM_D) break;
	if (debug) printf ("img2hist: first_bin=%d\tlast_bin=%d\n", i, j);
	*minval +=	(float) (range * i) / nbin;
	range   *=	(float) (j - i) / nbin;
	*maxval  =	*minval + range;

	for (i = 0; i < nvox; i++) {
		j = (int) ((imgt[i] - *minval) * nbin / (*maxval - *minval));
		if (j > 0 && j < nbin) hist[j]++;
	}
	free (marg);
}

