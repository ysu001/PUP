/*$Header: /data/petsun4/data1/src_solaris/algebra_4dfp/RCS/spatial_corr_4dfp.c,v 1.3 2012/06/07 23:01:13 avi Exp $*/
/*$Log: spatial_corr_4dfp.c,v $
 * Revision 1.3  2012/06/07  23:01:13  avi
 * option -M
 *
 * Revision 1.2  2010/05/03  21:08:27  avi
 * correct usage; option -c
 *
 * Revision 1.1  2010/03/17  20:31:44  avi
 * Initial revision
 **/

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
static char	rcsid[] = "$Id: spatial_corr_4dfp.c,v 1.3 2012/06/07 23:01:13 avi Exp $";
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

void usage () {
	fprintf (stderr, "Usage:	%s <image_x> <mask_x> <image_y> <mask_y> [output_text_file]\n", program);
	fprintf (stderr, "	option\n");
	fprintf (stderr, "	-C	compute covariance (default correlation)\n");
	fprintf (stderr, "	-M	suppress removal of patial means\n");
	fprintf (stderr, "	-c<flt>	scale output covariance matrix values by specified factor\n");
	fprintf (stderr, "N.B.:	image dimensions must match\n");
	fprintf (stderr, "N.B.:	%s counts only defined (not NaN or 1.e-37 or 0.0) voxels\n", program);
	exit (1);
}

int main (int argc, char *argv[]) {
	FILE		*fp;
	IFH		ifh;
	char		imgroot[4][MAXL], imgfile[MAXL];
	char		outfile[MAXL] = "";
	float		**img1;
	int		*imgm;
	float		voxdim[3], voxdim1[3], scale = 1.0;
	int		orient, orient1, imgdim[4][4], nvox, vdim, isbig[4];
	double		meanx, meany, top, botx, boty, r, cov;

/***********/
/* utility */
/***********/
	int		c, i, j, k, m;
	char		*ptr, command[MAXL];

/*********/
/* flags */
/*********/
	int		defined;
	int		debug = 0;
	int		status = 0;
	int		corr = 1;		/* governs covariance or correlation computation */
	int		mean_out = 1;		/* flag removal of spatial means */

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'C': corr = 0;			break;
				case 'M': mean_out = 0;			break;
				case 'd': debug++;			break;
				case 'c': scale = atof (ptr);		*ptr = '\0';	break;
			}
		} else switch (k) {
			case 0:
			case 1:
			case 2:
			case 3:	getroot (argv[i], imgroot[k++]);	break;
			case 4:	strcpy (outfile, argv[i]); k++;		break;
		}	
	}
	if (k < 4) usage ();

/*********************************/
/* check dimensional consistency */
/*********************************/
	for (i = 0; i < 4; i++) {
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
	img1 = calloc_float2 (4, vdim);
	if (!(imgm = (int *) calloc (vdim, sizeof (int)))) errm (program);

/*************/
/* read data */
/*************/
	for (i = 0; i < 4; i++) {
		sprintf (imgfile, "%s.4dfp.img", imgroot[i]);
		printf ("Reading: %s\n", imgfile);
		if (!(fp = fopen (imgfile, "rb")) || eread (img1[i], vdim, isbig[i], fp)
		|| fclose (fp)) errr (program, imgfile);
	}

/****************************************************/
/* compute intersection of masks and defined voxels */
/****************************************************/
	for (j = 0; j < vdim; j++) {
		imgm[j] = 1;
		for (i = 0; i < 4; i++) {
			defined = isnormal (img1[i][j]) && img1[i][j] != 0.0 && img1[i][j] != (float) 1.e-37;
			if (!defined) imgm[j] = 0;
		}
	}

/****************/
/* remove means */
/****************/
	for (meanx = meany = nvox = j = 0; j < vdim; j++) if (imgm[j]) {
		meanx += img1[0][j];
		meany += img1[2][j];
		nvox++;
	}
	meanx /= nvox;
	meany /= nvox;
	if (debug) printf ("nvox=%d meanx=%f meany=%f\n", nvox, meanx, meany);
	if (mean_out) for (j = 0; j < vdim; j++) if (imgm[j]) {
		img1[0][j] -= meanx;
		img1[2][j] -= meany;
	}

/**************************/
/* open named output file */
/**************************/
	if (strlen (outfile)) {
		if (!(fp = fopen (outfile, "a"))) errw (program, outfile);
		printf ("Writing %s\n", outfile);
	}

/********************************/
/* compute products over voxels */
/********************************/
	for (top = j = 0; j < vdim; j++) if (imgm[j]) {
		top += img1[0][j]*img1[2][j];
	}
	top /= nvox;

	if (!corr) {
		cov = scale*top;
			     printf ("%s\t%s\tcov=%.6f\n", imgroot[0], imgroot[2], cov);
		if (strlen (outfile)) {
			fprintf (fp, "%s\t%s\tcov=%.6f\n", imgroot[0], imgroot[2], cov);
			if (fclose (fp)) errw (program, outfile);
		}
		goto DONE;
	}

	for (botx = boty = j = 0; j < vdim; j++) if (imgm[j]) {
		botx += img1[0][j]*img1[0][j];
		boty += img1[2][j]*img1[2][j];
	}
	botx /= nvox;
	boty /= nvox;
	r = top/sqrt(botx*boty);
		         printf ("%s\t%s\tr=%.6f\n", imgroot[0], imgroot[2], r);
		if (strlen (outfile)) {
			fprintf (fp, "%s\t%s\tr=%.6f\n", imgroot[0], imgroot[2], r);
			if (fclose (fp)) errw (program, outfile);
		}

DONE:	free_float2 (img1);
	free (imgm);
	exit (status);
}
