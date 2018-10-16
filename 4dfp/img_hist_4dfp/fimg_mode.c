/**************************************************************************************/
/* Copyright  2000, 2001, 2002, 2003, 2004, 2005, 2006                                */
/* Washington University, Mallinckrodt Institute of Radiology.                        */
/* All Rights Reserved.                                                               */
/* This software may not be reproduced, copied, or distributed without written        */
/* permission of Washington University. For further information contact A. Z. Snyder. */
/**************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#define NBIN	128
#define CUM_D	0.001

float fimg_mode (float* fimg, int nval) {
	int		hist[NBIN], marg[NBIN];
	int		nbin = NBIN;
	int		i, j, k, imax;
	int		debug = 0;
	float		u, fmax, fmin, range, binwid;
	double		slope, curva;

	fmax = -FLT_MAX;
	fmin =  FLT_MAX;
	for (u = i = 0; i < nval; i++) {
		if (fimg[i] > fmax) fmax = fimg[i];
		if (fimg[i] < fmin) fmin = fimg[i];
		u += fimg[i];
	}
	u /= nval;
	range = fmax - fmin;
	if (debug) printf ("fmin=%f fmax=%f range=%f\n", fmin, fmax, range);

	for (k = 0; k < nbin; k++) marg[k] = hist[k] = 0;
	for (i = 0; i < nval; i++) {
		k = (nbin * (fimg[i] - fmin)) / range;
		if (k < nbin) marg[k]++;
	}
	for (k = 1; k < nbin; k++) marg[k] += marg[k - 1];
	if (debug) for (k = 0; k < nbin; k++) printf ("%-10d%10d\n", k, marg[k]);

	for (i = 1; i < nbin; i++) if ((float) marg[i] / (float) marg[nbin - 1] >       CUM_D) break;
	for (j = 1; j < nbin; j++) if ((float) marg[j] / (float) marg[nbin - 1] > 1.0 - CUM_D) break;
	if (debug) printf ("i=%d j=%d\n", i, j);
	fmin = range * i / nbin;
	fmax = range * j / nbin;
	range = fmax - fmin;
	binwid = range / nbin;
	if (debug) printf ("fmin=%f fmax=%f range=%f\n", fmin, fmax, range);

	for (k = 0; k < nval; k++) {
		i = (fimg[k] - fmin) / binwid;
		if (i < 0 || i >= nbin) continue;
		hist[i]++;
	}
	if (debug) for (k = 0; k < nbin; k++) printf ("%-10d%10d\n", k, hist[k]);

	imax = 0;
	for (i = 1; i < nbin - 1; i++) {
		if (hist[i] > hist[i - 1] && hist[i] > hist[i + 1]) {
			if (debug) printf ("peak at i=%d\n", i);
			if (!imax || hist[i] > hist[imax]) imax = i;
		}
	}
	if (debug) printf ("imax=%d\n", imax);
	if (imax == 0 || imax == nbin - 1) {
		fprintf (stdout, "fimg_mode: returning mean instead of mode\n");
		return u;
	}

/***************************/
/* parabolic interpolation */
/***************************/
	slope = 0.5 * (hist[imax + 1] - hist[imax - 1]);
	curva = hist[imax + 1] + hist[imax - 1] - 2*hist[imax];
	u = fmin + binwid*imax - binwid*slope/curva;
	if (debug) printf ("slope=%.2f curva=%.2f binwid=%f mode=%f\n", slope, curva, binwid, u);
	return u;
}

void hsmooth (int *hist, int nbin, int nsmooth) {
	int		i, j;
	float		fhist[NBIN], ghist[NBIN];

	for (i = 0; i < NBIN; i++) fhist[i] = hist[i];
	for (j = 0; j < nsmooth; j++) {
		for (i = 0; i < NBIN; i++) ghist[i] = fhist[i];
		for (i = 1; i < NBIN - 1; i++) fhist[i] = 0.25*(ghist[i-1] + 2.*ghist[i] + ghist[i+1]);
	}
	for (i = 0; i < NBIN; i++) hist[i] = (int) fhist[i] + 0.5;
}

/***********************************************************/
/* does not work correctly when compiled with optimizer -O */
/***********************************************************/
float fimg_modetn (float* fimg, int nval, float histmin, int nsmooth) {
	int		hist[NBIN], marg[NBIN];
	int		nbin = NBIN;
	int		i, j, k, imax;
	int		debug = 0;
	float		u, fmax, fmin, range, binwid;
	double		slope, curva;

	fmax = -FLT_MAX;
	fmin =  FLT_MAX;
	for (u = i = 0; i < nval; i++) {
		if (fimg[i] > fmax) fmax = fimg[i];
		if (fimg[i] < fmin) fmin = fimg[i];
		u += fimg[i];
	}
	u /= nval;
	range = fmax - fmin;
	if (debug) printf ("fmin=%f fmax=%f range=%f\n", fmin, fmax, range);

	for (k = 0; k < nbin; k++) marg[k] = hist[k] = 0;
	for (i = 0; i < nval; i++) {
		k = (nbin * (fimg[i] - fmin)) / range;
		if (k < nbin) marg[k]++;
	}
	for (k = 1; k < nbin; k++) marg[k] += marg[k - 1];
	if (debug) for (k = 0; k < nbin; k++) printf ("%-10d%10d\n", k, marg[k]);

	for (i = 1; i < nbin; i++) if ((float) marg[i] / (float) marg[nbin - 1] >       CUM_D) break;
	for (j = 1; j < nbin; j++) if ((float) marg[j] / (float) marg[nbin - 1] > 1.0 - CUM_D) break;
	if (debug) printf ("i=%d j=%d\n", i, j);
	fmin = range * i / nbin;
	fmax = range * j / nbin;
	range = fmax - fmin;
	binwid = range / nbin;
	if (debug) printf ("fmin=%f fmax=%f range=%f\n", fmin, fmax, range);

	for (k = 0; k < nval; k++) {
		i = (fimg[k] - fmin) / binwid;
		if (i < 0 || i >= nbin) continue;
		hist[i]++;
	}
	if (debug) for (k = 0; k < nbin; k++) printf ("%-10d%10d\n", k, hist[k]);

	hsmooth (hist, nbin, nsmooth);
	imax = 0;
	for (i = 1; i < nbin - 1; i++) {
		u = fmin + binwid*i;
		if (hist[i] > hist[i - 1] && hist[i] > hist[i + 1] && u > histmin) {
			if (debug) printf ("peak at bin %d\tordinate=%.2f\t%10d\n", i, u, hist[i]);
			if (!imax || hist[i] > hist[imax]) imax = i;
		}
	}
	if (debug) printf ("imax=%d\n", imax);
	if (imax == 0 || imax == nbin - 1) {
		fprintf (stdout, "fimg_mode: returning mean instead of mode\n");
		return u;
	}

/***************************/
/* parabolic interpolation */
/***************************/
	slope = 0.5 * (hist[imax + 1] - hist[imax - 1]);
	curva = hist[imax + 1] + hist[imax - 1] - 2*hist[imax];
	u = fmin + binwid*imax - binwid*slope/curva;
	if (debug) printf ("slope=%.2f curva=%.2f binwid=%f mode=%f\n", slope, curva, binwid, u);
	return u;
}
