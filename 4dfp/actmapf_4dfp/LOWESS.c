/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/LOWESS.c,v 1.3 2011/01/16 02:51:42 avi Exp $*/
/*$Log: LOWESS.c,v $
 * Revision 1.3  2011/01/16  02:51:42  avi
 * in output order y_fit, y_orig -> y_orig, y_fit
 *
 * Revision 1.2  2010/07/23  21:11:11  avi
 * correct W()
 *
 * Revision 1.1  2010/07/23  04:37:34  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>			/* getpid () */

#define MAXL		256
#define MAX(a,b)	(a>b? a:b)
#define MIN(a,b)	(a<b? a:b)

/*************/
/* externals */
/*************/
void	dmatinv_ (double *a, int *n, double *det);				/* dmatinv.f */

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: LOWESS.c,v 1.3 2011/01/16 02:51:42 avi Exp $";
static char	program[MAXL];
static int	debug = 0;

typedef struct {
	double	xorig, x, y, yhat, h, err;
	double	beta[2];
	int	fit;
} POINT;

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

double W (double x) {
	double q;

	q = 1. - x*x;
	return (q > 0.) ? q*q : 0.;
}

static int sortbyx (const void *pi, const void *pj) {
	POINT	*i, *j;
	double	q;

	i = (POINT *) pi;
	j = (POINT *) pj;
	q = i->xorig - j->xorig;
	return (q < 0.) ? -1 : 1;
}

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n#%s\n", rcsid);
}

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;

	if (ptr = strchr (string, '#')) *ptr = '\0';
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

void usage (char *program) {
	printf ("Usage:\t%s <profile>\n", program);
	printf ("\toption\n");
	printf ("\t-r<int>	specify LOWESS fitting radius (default = npts/2)\n");
	exit (1);
}

int main (int argc, char *argv[]) {
/**********************/
/* filename variables */
/**********************/
	FILE		*profp, *outfp;
	char		profile[MAXL], outfile[MAXL];
	float		fmin, fmax, fmode;
	double		A[4], C[2], det, w;
	int		npts, ncol, r = 0, kmax, imax;
	POINT		*points, point;

/***********/
/* utility */
/***********/
	char		command[MAXL], string[MAXL], *ptr, *srgv[MAXL];
	int		c, i, j, k, m;
	int		two = 2;
	double		q;

/*********/
/* flags */
/*********/
	int		status = 0;
	int		threecol = 0;

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
				case 'd': debug++;			break;
				case '3': threecol++;			break;
				case 'r': r = atoi (ptr);		*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	strcpy  (profile, argv[i]);		k++; break;
		}	
	}
	if (k < 1) usage (program);

/****************/
/* read profile */
/****************/
	printf ("Reading: %s\n", profile);
	if (!(profp = fopen (profile, "r"))) errr (program, profile);
	npts = 0; while (fgets (string, MAXL, profp)) {
		if (!(m = split (string, srgv, MAXL))) continue;
		if (!npts) {
			ncol = m;
			printf ("ncol=%d\n", ncol);
		} else {
			if (m != ncol) {
				fprintf (stderr, "%s: %s format error\n", program, profile);
				exit (-1);
			}
		}
		npts++;
	}
	printf ("npts=%d\n", npts);

	if (!(points =	(POINT *) calloc (npts, sizeof (POINT)))) errm (program);
	rewind (profp);
	for (i = 0; i < npts; i++) {
		fgets (string, MAXL, profp);
		if (!(m = split (string, srgv, MAXL))) continue;
		points[i].xorig 	= atof (srgv[0]);
		points[i].y		= atof (srgv[1]);
	}
	fclose (profp);

	qsort ((void *) points, npts, sizeof (POINT), sortbyx);
	if (debug) for (i = 0; i < npts; i++) printf ("%10.4f%10.4f\n", points[i].xorig, points[i].y);
	for (i = 0; i < npts; i++) {
		points[i].x = -1. + 2.*(points[i].xorig - points[0].xorig)/(points[npts-1].xorig - points[0].xorig);
	}

	if (r == 0 || r < 0 || r > npts/2) {
		imax = npts/2;
		r = imax - 1;
	} else {
		imax = npts - r;
	}

	for (i = r; i < imax; i++) {
		for (j = 0; j < 4; j++) A[j] =0.;
		for (j = 0; j < 2; j++) C[j] =0.;
		kmax = (r == imax - 1) ? npts - 1 : i + r;
		points[i].h = MAX(fabs(points[i].x - points[kmax].x), fabs(points[i].x - points[i - r].x));
		for (k = i - r; k <= kmax; k++) {
			if (r == imax - 1) {
				w = 1.;
			} else {
				w = W ((points[k].x - points[i].x)/points[i].h);
			}
			A[0] += w;
			A[1] += w*points[k].x;
			A[2] =  A[1];
			A[3] += w*points[k].x*points[k].x;
			C[0] += w*points[k].y;
			C[1] += w*points[k].x*points[k].y;
		}
		dmatinv_ (A, &two, &det);
		for (j = 0; j < 2; j++) {
			points[i].beta[j] = A[j]*C[0] + A[j + 2]*C[1];
		}
		points[i].fit = 1;
	}

	sprintf (outfile, "%s.fit", profile);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);
	write_command_line (outfp, argc, argv);
	fprintf (outfp, "#%9s%10s%10s\n", "x", "y_orig", "y_fit");
	for (i = 0; i < npts; i++) {
		if (points[i].fit) {
				points[i].yhat = points[i].beta[0]	+ points[i].beta[1]*points[i].x;
		} else {
			if (i < r) {
				points[i].yhat = points[r].beta[0]	+ points[r].beta[1]*points[i].x;
			} else {
				points[i].yhat = points[imax-1].beta[0] + points[imax-1].beta[1]*points[i].x;
			}
		}
		fprintf (outfp, "%10.4f%10.4f%10.4f\n", points[i].xorig, points[i].y, points[i].yhat);
	}
	fclose (outfp);

	free (points);
	exit (status);
}
