/**************************************************************************************/
/* Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007			      */
/* Washington University, Mallinckrodt Institute of Radiology.			      */
/* All Rights Reserved.                                                               */
/* This software may not be reproduced, copied, or distributed without written        */
/* permission of Washington University. For further information contact A. Z. Snyder. */
/**************************************************************************************/
/*$Header: /home/usr/shimonyj/diff4dfp/RCS/get_dti_params.c,v 1.4 2012/06/15 01:30:37 avi Exp $*/
/*$Log: get_dti_params.c,v $
 * Revision 1.4  2012/06/15  01:30:37  avi
 * tolerate (meaningless) vector components on b = 0 input lines
 *
 * Revision 1.3  2007/08/30  05:09:08  avi
 * JSSutil.h compliant
 *
 * Revision 1.2  2000/12/19  01:46:58  avi
 * copyright
 *
 * Revision 1.1  2000/10/05  05:46:35  avi
 * Initial revision
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <JSSutil.h>

void get_dti_nrc_test () {
	extern int	dti_dimen (char *file);
	extern int	get_dti_params_nrc (char *file, int n, float *b_vals, float **q_vals);
	int		i, j, k, n;
	float		**q_vals, *b_vals;
	char		file[] = "/home/usr/shimonyj/diff4dfp/tp7_params.dat";

	if ((n = dti_dimen (file)) <= 0) exit (-1);
	b_vals = vector (1, n);
	q_vals = matrix (1, n, 1, 3);

	get_dti_params_nrc (file, n, b_vals, q_vals);

	free_matrix (q_vals, 1, n, 1, 3);
	free_vector (b_vals, 1, n);
	exit (0);
}

int dti_dimen (char *file) {
	FILE		*fp;
	char		*ptr, string[256];
	int		i, k, n;

	if (access (file, R_OK)	|| !(fp = fopen (file, "r"))) return -1;
	while (1) {
		if (!(fgets (string, 256, fp))) {
			fclose (fp); return -1;
		}
		if (ptr = strchr (string, '#')) *ptr = '\0';
		if (k = sscanf (string, "%d", &n) != 1) continue;
		fclose (fp);
		printf ("dti_dimen=%d\n", n);
		return n;
	}
}

int get_dti_params_nrc (char *file, int n, float *b_vals, float **q_vals) {
	FILE		*fp;
	char		*ptr, string[256];
	float		v[4];
	double		q;
	int		i, j, k, l;

	if (access (file, R_OK)	|| !(fp = fopen (file, "r"))) return -1;
	i = 0; while (i <= n) {
		if (!(fgets (string, 256, fp))) goto ERRF;
		if (ptr = strchr (string, '#')) *ptr = '\0';
		if (!i) {
			if (k = sscanf (string, "%d", &j) != 1) continue;
			if (j != n) goto ERRF;
		} else {
			k = sscanf (string, "%f%f%f%f", v + 0, v + 1, v + 2, v + 3);
			if (!k) continue;
			b_vals[i] = v[0] / 1000.;
			if (k == 1 || v[0] == 0.) {
				for (j = 1; j <= 3; j++) q_vals[i][j] = 0.0;
			} else if (k == 4) {
				q = 0;
				for (j = 1; j <= 3; j++) q += v[j] * v[j];
				for (j = 1; j <= 3; j++) q_vals[i][j] = v[j] / sqrt (q);
			} else {
				fprintf (stderr, "input line field count not 1 or 4\n");
				return -1;
			}
			printf ("%10.4f%10.4f%10.4f%10.4f\n", b_vals[i], q_vals[i][1], q_vals[i][2], q_vals[i][3]);
		}
		i++;
	}
	fclose (fp);
	return 0;
ERRF:	fprintf (stderr, "%s parse error\n", file);
	fclose (fp);
	return -1;
}
