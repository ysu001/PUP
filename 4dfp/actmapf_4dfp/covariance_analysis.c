/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/covariance_analysis.c,v 1.4 2006/11/25 23:49:43 avi Exp $*/
/*$Log: covariance_analysis.c,v $
 * Revision 1.4  2006/11/25  23:49:43  avi
 * Solaris 10
 *
 * Revision 1.3  2005/12/12  06:47:59  avi
 * obtain values nsub, mlag, nreg from input
 * also output ACR averaged over regions
 *
 * Revision 1.2  2005/05/25  23:43:56  justinv
 * nsub, nreg, and max lag options
 *
 * Revision 1.1  2005/02/01  21:58:47  avi
 * Initial revision
 **/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>			/* R_OK */
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#define MAXL		1024		/* maximum input line length */

void errm (char* program) {
	fprintf (stderr, "%s: memory allocation error\n", program);
	exit (-1);
}

void errr (char* program, char* filespc) {
	fprintf (stderr, "%s: %s read error\n", program, filespc);
	exit (-1);
}

void errw (char* program, char* filespc) {
	fprintf (stderr, "%s: %s write error\n", program, filespc);
	exit (-1);
}

void errf (char* program, char* filespc) {
	fprintf (stderr, "%s: %s format error\n", program, filespc);
	exit (-1);
}

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;

	i = m = 0;
	while (m < maxp) {
		while (!isgraph ((int) string[i]) && string[i]) i++;
		if (!string[i]) break;
		srgv[m++] = string + i;
		while (isgraph ((int) string[i])) i++;
		if (!string[i]) break;
		string[i++] = '\0';
	}
	return m;
}

void getroot (char *filespc, char *filroot) {
	char	*str;
	strcpy (filroot, filespc);
	while (str = strrchr (filroot, '.')) {
			if (!strcmp (str, ".rec"))	*str = '\0';
		else	if (!strcmp (str, ".dat"))	*str = '\0';
		else	if (!strcmp (str, ".txt"))	*str = '\0';
		else	if (!strcmp (str, ".lst"))	*str = '\0';
		else	break;
	}
}

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: covariance_analysis.c,v 1.4 2006/11/25 23:49:43 avi Exp $";
static char	program[MAXL];
static int	debug = 0;

float ***calloc_float3 (int n1, int n2, int n3) {
	int	i, j;
	float	***a;

	if (!(a = (float ***) malloc (n1 * sizeof (float **)))) errm (program);
	for (i = 0; i < n1; i++) {
		if (!(a[i] = (float **) malloc (n2 * sizeof (float *)))) errm (program);
		for (j = 0; j < n2; j++) {
			if (!(a[i][j] = (float *) calloc (n3, sizeof (float)))) errm (program);
		}
	}
	return a;
}

float **calloc_float2 (int n1, int n2) {
	int	i;
	float	**a;

	if (!(a = (float **) malloc (n1 * sizeof (float *)))) errm (program);
	for (i = 0; i < n1; i++) {
		if (!(a[i] = (float *) calloc (n2, sizeof (float)))) errm (program);
	}
	return a;
}

char **calloc_char2 (int n1, int n2) {
	int	i;
	char	**a;

	if (!(a = (char **) malloc (n1 * sizeof (char *)))) errm (program);
	for (i = 0; i < n1; i++) {
		if (!(a[i] = (char *) calloc (n2, sizeof (char)))) errm (program);
	}
	return a;
}

void free3 (void ***a, int n1, int n2) {
	int	i, j;

	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) free (a[i][j]);
		free (a[i]);
	}
	free (a);
}

void free2 (void **a, int n1) {
	int	i;

	for (i = 0; i < n1; i++) free (a[i]);
	free (a);
}

void usage (char *program) {
	fprintf (stderr, "Usage:\t%s <lstfile>\n", program);
	fprintf (stderr, "\toption\n");
	exit (1);
}

int main (int argc, char *argv[]) {
/**********************/
/* filename variables */
/**********************/
	FILE		*tmpfp, *ACRfp, *outfp, *lstfp;
	char		outroot[MAXL], outfile[MAXL], lstfile[MAXL], ACRfile[MAXL];
	char		**ROI;			/* ROI[nreg][32];		32 char ROI names */

/*************************/
/* timeseries processing */
/*************************/
	int		mlag;
	int		isub, nsub;
	int		ireg, jreg, nreg;
	float		***ACR;			/* ACR[nsub][nreg][2*mlag + 1];	input autocorrelation */
	float		**ACRbar;		/* ACRbar[nreg][2*m + 1];	autocorrelation average over subjects */
	float		*gammabar, sum;		/* gammabar[2*mlag + 1];	Bartlett's gamma */
	float		*ACRbar1;		/* autocorrelation averaged over subjects and ROI */

/***********/
/* utility */
/***********/
	char		command[MAXL], string[MAXL], *ptr, *srgv[MAXL];
	int		c, i, j, k, l, m;

/*********/
/* flags */
/*********/
	int		status = 0;

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'd': debug++;		break;
			}
		}
		else switch (k) {
			case 0:	strcpy (lstfile, argv[i]);	k++; break;
		}	
	}
	if (k < 1) usage (program);

/*****************************/
/* read lstfile and ACRfiles */
/*****************************/
	printf ("Reading: %s\n", lstfile);
	if (!(lstfp = fopen (lstfile, "r"))) errr (program, lstfile);
/***********************************************/
/* count number of subjects = listed ACR files */
/***********************************************/
	nsub = 0;
	while (fgets (ACRfile, MAXL, lstfp)) {
		ptr = strrchr (ACRfile, '\n'); *ptr = '\0'; /* strip off trailing '\n' */
		if (access (ACRfile, R_OK)) errr (program, ACRfile);
		nsub++;
	}
	rewind (lstfp);

/**********************************************/
/* determine nreg and mlag from first ACRfile */
/**********************************************/
	fgets (ACRfile, MAXL, lstfp);
	ptr = strrchr (ACRfile, '\n'); *ptr = '\0'; /* strip off trailing '\n' */
	if (!(ACRfp = fopen (ACRfile, "r"))) errr (program, ACRfile);
	while (fgets (string, MAXL, ACRfp) && !mlag) {
		if (string[0] == '#') {
			ptr = string; while (ptr = strpbrk (ptr, "#")) *ptr = '\x20'; /* convert '#' to space */
			nreg = split (string, srgv, MAXL) - 1;
		} else {
			m = split (string, srgv, MAXL);
			if (m != nreg + 1) errf (program, ACRfile);
			mlag = -atoi (srgv[0]);
			if (mlag <= 0) errf (program, ACRfile);
		}
	}
	fclose (ACRfp);
	rewind (lstfp);

	printf ("nsub=%d nreg%d mlag=%d\n", nsub, nreg, mlag);
/************/
/* allocate */
/************/
	ROI	= calloc_char2  (nreg, 32);
	ACR	= calloc_float3 (nsub, nreg, 2*mlag + 1);
	ACRbar	= calloc_float2 (nreg, 2*mlag + 1);
	if (!(gammabar = (float *) calloc (2*mlag + 1, sizeof (float)))) errm (program);

/************************/
/* read ACR[isub][j][l] */
/************************/
	for (isub = 0; isub < nsub; isub++) {
		fgets (ACRfile, MAXL, lstfp);
		ptr = strrchr (ACRfile, '\n'); *ptr = '\0'; /* strip off trailing '\n' */
		if (!(ACRfp = fopen (ACRfile, "r"))) errr (program, ACRfile);
		printf ("Reading: %s\n", ACRfile);
		l = -mlag; while (fgets (string, MAXL, ACRfp)) {
			if (string[0] == '#') {
				printf ("%s", string);
				ptr = string; while (ptr = strpbrk (ptr, "#")) *ptr = '\x20'; /* convert '#' to space */
				k = split (string, srgv, MAXL) - 1;
				if (k == nreg) for (j = 0; j < nreg; j++) strncpy (ROI[j], srgv[j + 1], 32);
			} else {
				m = split (string, srgv, MAXL);
				if (m != nreg + 1) errf (program, ACRfile);
				if (atoi (srgv[0]) != l) errf (program, ACRfile);
				for (j = 0; j < nreg; j++) {
					ACR[isub][j][l + mlag] = atof (srgv[j + 1]);
				}
				l++;
			}
		}
		fclose (ACRfp);
	}
	fclose (lstfp);

/*****************************/
/* average ACR over subjects */
/*****************************/
	for (ireg = 0; ireg < nreg; ireg++) for (l = 0; l < 2*mlag + 1; l++) {
		for (isub = 0; isub < nsub; isub++) ACRbar[ireg][l] += ACR[isub][ireg][l];
		ACRbar[ireg][l] /= nsub;
	}
	getroot (lstfile, outroot);
	sprintf (outfile, "%s_avg.dat", outroot);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	fprintf (stderr, "Writing: %s\n", outfile);
	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]); fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
	fprintf (outfp, "#%9s", "lag");
	for (i = 0; i < nreg; i++) fprintf (outfp, "%10s", ROI[i]); fprintf (outfp, "\n");
	for (l = 0; l < 2*mlag + 1; l++) {
		fprintf (outfp, "%10d", -mlag + l);
		for (ireg = 0; ireg < nreg; ireg++) fprintf (outfp, "%10.4f", ACRbar[ireg][l]);
		fprintf (outfp, "\n");
	}
	if (fclose (outfp)) errw (program, outfile);

/*******************************/
/* average ACRbar over regions */
/*******************************/
	if (!(ACRbar1 = (float *) calloc (2*mlag + 1, sizeof (float)))) errm (program);
	for (l = 0; l < 2*mlag + 1; l++) {
		for (ireg = 0; ireg < nreg; ireg++) ACRbar1[l] += ACRbar[ireg][l];
		ACRbar1[l] /= nreg;
	}
	sprintf (outfile, "%s_ROIavg.dat", outroot);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	fprintf (stderr, "Writing: %s\n", outfile);
	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]); fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
	fprintf (outfp, "#%9s%10s\n", "lag", "ROIavg");
	for (l = 0; l < 2*mlag + 1; l++) {
		fprintf (outfp, "%10d%10.4f\n", -mlag + l, ACRbar1[l]);
	}
	if (fclose (outfp)) errw (program, outfile);
	free (ACRbar1);

/****************************/
/* compute Bartlett's gamma */
/****************************/
	for (sum = l = 0; l < 2*mlag + 1; l++) {
		for (ireg = 0; ireg < nreg; ireg++) for (jreg = 0; jreg < nreg; jreg++) {
			if (jreg == ireg) continue;
			gammabar[l] += ACRbar[ireg][l]*ACRbar[jreg][l];
		}
		gammabar[l] /= nreg*(nreg - 1);
		sum += gammabar[l];
	}
	getroot (lstfile, outroot);
	sprintf (outfile, "%s_gammabar.dat", outroot);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	fprintf (stderr, "Writing: %s\n", outfile);
	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]); fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
	fprintf (outfp, "#sum of gammabar = %.4f\n", sum);
	for (l = 0; l < 2*mlag + 1; l++) {
		fprintf (outfp, "%10d", -mlag + l);
		fprintf (outfp, "%10.4f\n", gammabar[l]);
	}
	if (fclose (outfp)) errw (program, outfile);
	printf ("Bartlett correction factor (sum of gammabar) = %.4f\n", sum);

	free (gammabar);
	free2 ((void *) ROI, nreg);
	free3 ((void *) ACR, nsub, nreg);
	free2 ((void *) ACRbar, nreg);
	exit (status);
}
