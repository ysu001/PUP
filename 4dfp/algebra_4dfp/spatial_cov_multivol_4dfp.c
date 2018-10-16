/*$Header: /data/petsun4/data1/src_solaris/algebra_4dfp/RCS/spatial_cov_multivol_4dfp.c,v 1.2 2010/07/08 04:07:30 avi Exp $*/
/*$Log*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>		/* R_OK */
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#ifndef _FEATURES_H		/* defined by math.h in Linux */
	#include <ieeefp.h>
	#include <sunmath.h>	/* isnormal() */
#endif
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MAXL		256	/* maximum string length */

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: spatial_cov_multivol_4dfp.c,v 1.2 2010/07/08 04:07:30 avi Exp $";
static char	program[MAXL];
static int	debug = 0;

float **calloc_float2 (int n1, int n2) {
	int	i;
	float	**a;

	if (!(a = (float **) malloc (n1 * sizeof (float *)))) errm (program);
	if (!(a[0] = (float *) calloc (n1 * n2, sizeof (float)))) errm (program);
	for (i = 1; i < n1; i++) a[i] = a[0] + i*n2;
	return a;
}

void free_float2 (float **a) {
	free (a[0]);
	free (a);
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
}

void usage () {
	fprintf (stderr, "Usage:	%s <(4dfp) image> <(4dfp) mask>\n", program);
	fprintf (stderr, "	option\n");
	fprintf (stderr, "	-Z	compute covariance with respect to zero (default wrt image mean)\n");
	fprintf (stderr, "	-c<flt>	scale output covariance matrix values by specified factor\n");
	fprintf (stderr, "N.B.:	%s counts only defined (not NaN or 1.e-37 or 0.0) voxels\n", program);
	exit (1);
}

int main (int argc, char *argv[]) {
	FILE		*fp;
	IFH		ifh;
	char		imgroot[2][MAXL], imgfile[MAXL];
	char		outfile[MAXL];
	float		**img1, *mean, **cov, scale = 1.0;
	int		*imgm;
	float		voxdim[3], voxdim1[3];
	int		orient, orient1, imgdim[2][4], nvox, vdim, isbig[2];

/***********/
/* utility */
/***********/
	int		c, i, j, k, l;
	char		*ptr, command[MAXL];
	double		q;

/*********/
/* flags */
/*********/
	int		defined;
	int		debug = 0;
	int		Z_flag = 0;
	int		status = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'd': debug++;			break;
				case 'Z': Z_flag++;			break;
				case 'c': scale = atof (ptr);		*ptr = '\0';	break;
			}
		} else switch (k) {
			case 0:
			case 1:	getroot (argv[i], imgroot[k++]);	break;
		}	
	}
	if (k < 2) usage ();

/*********************************/
/* check dimensional consistency */
/*********************************/
	for (i = 0; i < 2; i++) {
		printf ("input %3d: %s\n", i + 1, imgroot[i]);
		if (!i) {
			if (get_4dfp_dimoe (imgroot[i], imgdim[i],  voxdim, &orient,  isbig + i)) exit (-1);
			if (Getifh (imgroot[i], &ifh)) errr (program, imgroot[i]);
		} else {
			if (get_4dfp_dimoe (imgroot[i], imgdim[i], voxdim1, &orient1, isbig + i)) exit (-1);
			status = (orient1 != orient);
			for (k = 0; k < 3; k++) status |= (imgdim[i][k] != imgdim[0][k]);
			for (k = 0; k < 3; k++) status |= (fabs (voxdim1[k] - voxdim[k]) > 1.e-5);
		}
		if (status) {
			fprintf (stderr, "%s: %s %s nvox mismatch\n", program, imgroot[0], imgroot[i]);
			exit (-1);
		}
	}
	vdim = imgdim[0][0]*imgdim[0][1]*imgdim[0][2];
	img1 = calloc_float2 (2, vdim);
	cov  = calloc_float2 (imgdim[0][3], imgdim[0][3]);
	if (!(imgm = (int *) calloc (vdim, sizeof (int)))) errm (program);
	if (!(mean = (float *) calloc (imgdim[0][3], sizeof (float)))) errm (program);

/*************/
/* read mask */
/*************/
	sprintf (imgfile, "%s.4dfp.img", imgroot[1]);
	printf ("Reading: %s\n", imgfile);
	if (!(fp = fopen (imgfile, "rb")) || eread (img1[1], vdim, isbig[1], fp)
	|| fclose (fp)) errr (program, imgfile);
	for (j = 0; j < vdim; j++) {
		defined = isnormal (img1[1][j]) && img1[1][j] != 0.0 && img1[1][j] != (float) 1.e-37;
		imgm[j] = (defined) ? 1 : 0;
	}

/***************************/
/* open multi-volume image */
/***************************/
	sprintf (imgfile, "%s.4dfp.img", imgroot[0]);
	printf ("Reading: %s\n", imgfile);
	if (!(fp = fopen (imgfile, "rb"))) errr (program, imgfile);

/**************************/
/* compute defined voxels */
/**************************/
	for (k = 0; k < imgdim[0][3]; k++) {
		if (eread (img1[0], vdim, isbig[0], fp)) errr (program, imgfile);
		for (j = 0; j < vdim; j++) if (imgm[j]) {
			defined = isnormal (img1[1][j]) && img1[1][j] != 0.0 && img1[1][j] != (float) 1.e-37;
			if (!defined) imgm[j] = 0;
		}
	}
	rewind (fp);
	for (nvox = j = 0; j < vdim; j++) if (imgm[j]) nvox++;

/************************/
/* compute volume means */
/************************/
	if (!Z_flag) for (k = 0; k < imgdim[0][3]; k++) {
		if (eread (img1[0], vdim, isbig[0], fp)) errr (program, imgfile);
		for (q = j = 0; j < vdim; j++) if (imgm[j]) q += img1[0][j];
		mean[k] = q/nvox;
	}
	rewind (fp);

/******************************/
/* compute spatial covariance */
/******************************/
	printf ("computing covariance matrix row");
	for (k = 0; k < imgdim[0][3]; k++) {printf (" %d", k + 1); fflush (stdout);
		if (fseek (fp, (long) k*vdim*sizeof (float), SEEK_SET)
		||  eread (img1[0], vdim, isbig[0], fp)) errr (program, imgfile);
	for (l = k; l < imgdim[0][3]; l++) {
		if (fseek (fp, (long) l*vdim*sizeof (float), SEEK_SET)
		||  eread (img1[1], vdim, isbig[0], fp)) errr (program, imgfile);
		for (q = j = 0; j < vdim; j++) if (imgm[j]) q += (img1[0][j] - mean[k])*(img1[1][j] - mean[l]);
		cov[k][l] = cov[l][k] = q/nvox;
	}}
	printf ("\n");
	if (fclose (fp)) errr (program, imgfile);

/****************/
/* write output */
/****************/
	sprintf (outfile, "%s_covariance_matrix.dat", imgroot[0]);
	if (!(fp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);
	write_command_line (fp, argc, argv);
	for (k = 0; k < imgdim[0][3]; k++) {
		for (l = 0; l < imgdim[0][3]; l++) fprintf (fp, " %9.6f", scale*cov[k][l]);
		fprintf (fp, "\n");
	}
	if (fclose (fp)) errw (program, outfile);

	free_float2 (img1);
	free_float2 (cov);
	free (imgm);
	free (mean);
	exit (status);
}
