/*$Header: /data/petsun4/data1/src_solaris/algebra_4dfp/RCS/RFX_4dfp.c,v 1.2 2011/01/01 00:52:46 avi Exp $*/
/*$Log: RFX_4dfp.c,v $
 * Revision 1.2  2011/01/01  00:52:46  avi
 * better usage and cosmetics
 *
 * Revision 1.1  2010/12/31  07:26:32  avi
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
#define MAXR		4096	/* maximum number of input images */

float rnan (void) {
	union {
		float		r;
		unsigned int	j;
        } word;
	word.j = 0x7fffffff;
	return word.r;
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

/********************/
/* global variables */
/********************/
char		program[MAXL];
static char	rcsid[] = "$Id: RFX_4dfp.c,v 1.2 2011/01/01 00:52:46 avi Exp $";

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

double **calloc_double2 (int n1, int n2) {
	int	i;
	double	**a;

	if (!(a = (double **) malloc (n1 * sizeof (double *)))) errm (program);
	if (!(a[0] = (double *) calloc (n1 * n2, sizeof (double)))) errm (program);
	for (i = 1; i < n1; i++) a[i] = a[0] + i*n2;
	return a;
}

void free_double2 (double **a) {
	free (a[0]);
	free (a);
}

int **calloc_int2 (int n1, int n2) {
	int	i;
	int	**a;

	if (!(a = (int **) malloc (n1 * sizeof (int *)))) errm (program);
	if (!(a[0] = (int *) calloc (n1 * n2, sizeof (int)))) errm (program);
	for (i = 0; i < n1; i++) a[i] = a[0] + i*n2;
	return a;
}

void free_int2 (int **a) {
	free (a[0]);
	free (a);
}

void usage () {
	fprintf (stderr, "Usage:\t%s <list_group1> <(4dfp) Nimage_group1> [<list_group2> <(4dfp) Nimage_group2>]\n", program);
	fprintf (stderr, " e.g.:\t%s grp1_zfrm.lst grp1_N grp2_zfrm.lst grp2_N\n", program);
	fprintf (stderr, "   or:\t%s grp1_diff_zfrm.lst grp1_N\n", program);
	fprintf (stderr, "	option\n");
	fprintf (stderr, "	-o<str>	specify non-default output tmap 4dfp filename root\n");
	fprintf (stderr, "	-R	abbreviate output tmap rec file\n");
	fprintf (stderr, "	-N\toutput undefined voxels as NaN\n");
	fprintf (stderr, "	-Z\toutput undefined voxels as 0\n");
	fprintf (stderr, "	-E\toutput undefined voxels as 1.E-37 (default)\n");
	fprintf (stderr, "	-@<b|l>\toutput big or little endian (default first named <list_group1> endian)\n");
	fprintf (stderr, "	N.B.:	default output tmap fileroot is, e.g., grp1_zfrm_vs_grp2_zfrm_tmap\n");
	exit (1);
}

int main (int argc, char *argv[]) {
	FILE		*fp;
	IFH		ifh;
	char		imgroot[2][MAXR][MAXL], Nmgroot[2][MAXL], lstfile[2][MAXL];
	int		nimg[2], isbig[2][MAXR], isbigN[2];
	int		orient, orient1, imgdim[4], imgdim1[4], dimension, vdim;
	char		outroot[MAXL] = "", imgfile[MAXL], outfile[MAXL];
	float		voxdim[3], voxdim1[3];
	float		amin = FLT_MAX, amax = -FLT_MAX;
	float		*img1,  **imgN;
	double		**imgs, **imgv;
	int		**imgm, ngroup;
	char		control = '\0';

/***********/
/* utility */
/***********/
	double		u;
	int		c, i, j, k, m;
	char		*ptr, *str, command[MAXL];
	char		*srgv[MAXL];		/* string field pointers */

/*********/
/* flags */
/*********/
	int		defined;
	int		debug = 0;
	int		NaN_flag = 'E';		/* 'E' 1.e-37; 'Z' 0.0; 'N' NaN; */
	int		status = 0;
	int		saverec = 1;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	k = 0; for (i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'R': saverec = 0;			break;
				case 'd': debug++;			break;
				case 'o': getroot (ptr, outroot);	*ptr = '\0'; break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0:	strcpy (lstfile[0], argv[i]);		k++; break;
			case 1:	getroot (argv[i], Nmgroot[0]);		k++; break;
			case 2:	strcpy (lstfile[1], argv[i]);		k++; break;
			case 3:	getroot (argv[i], Nmgroot[1]);		k++; break;
		}	
	}
	if (k != 2 && k != 4) usage ();
	ngroup = k/2;

/********************/
/* parse list files */
/********************/
	for (k = 0; k < ngroup; k++) {
		nimg[k] = 0;
		if (!(fp = fopen (lstfile[k], "r"))) errr (program, lstfile[k]);
		printf ("Reading: %s\n", lstfile[k]);
		while (fgets (command, MAXL, fp)) {
			if (nimg[k] >= MAXR) {
				fprintf (stderr, "%s: maximum number of input images (%d) exceeded\n", program, MAXR);
				exit (-1);
			}
			if (ptr = strchr (command, '#'))  *ptr = '\0';
			if (!strlen (command)) continue;		/* skip blank lines */
			if (ptr = strchr (command, '\n')) *ptr = '\0';	/* strip terminal nl */
			i = m = 0; while (m < MAXL && i < MAXL) {
				while (!isgraph ((int) command[i]) && command[i]) i++;
				if (!command[i]) break;
				srgv[m++] = command + i;
				while (isgraph ((int) command[i])) i++;
				if (!command[i]) break;
				command[i++] = '\0';
			}
			getroot (srgv[0], imgroot[k][nimg[k]++]);
		}
		fclose (fp);
		if (!nimg[k]) usage ();
	
/*********************************/
/* check dimensional consistency */
/*********************************/
		for (i = 0; i < nimg[k]; i++) {
			printf ("image %3d: %s\n", i + 1, imgroot[k][i]);
			if (!i && !k) {
				if (get_4dfp_dimoe (imgroot[k][i], imgdim,  voxdim,  &orient,  isbig[k] + i)) exit (-1);
				if (Getifh (imgroot[k][i], &ifh)) errr (program, imgroot[k][i]);
			} else {
				if (get_4dfp_dimoe (imgroot[k][i], imgdim1, voxdim1, &orient1, isbig[k] + i)) exit (-1);
				status = (orient1 != orient);
				for (j = 0; j < 4; j++) status |= (imgdim1[j]      != imgdim[j]);
				for (j = 0; j < 3; j++) status |= (fabs (voxdim1[j] - voxdim[j]) > 1.e-4);
			}
			if (status) {
				fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot[0][0], imgroot[k][i]);
				exit (-1);
			}
		}
		if (get_4dfp_dimoe (Nmgroot[k], imgdim1, voxdim1, &orient1, isbigN + k)) exit (-1);
		status = (orient1 != orient);
		for (j = 0; j < 3; j++) status |= (imgdim1[j]      != imgdim[j]);
		for (j = 0; j < 3; j++) status |= (fabs (voxdim1[j] - voxdim[j]) > 1.e-4);
		if (status) {
			fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot[0][0], Nmgroot[k]);
			exit (-1);
		}
	}
	if (!control) control = (isbig[0][0]) ? 'b' : 'l';

	vdim = imgdim[0]*imgdim[1]*imgdim[2];
	dimension = vdim*imgdim[3];
	if (!(img1 = (float *) calloc (dimension, sizeof (float)))) errm (program);
	imgv = calloc_double2 (ngroup, dimension);	/* variance accumulator */
	imgs = calloc_double2 (ngroup, dimension);	/* sum      accumulator */
	imgm = calloc_int2    (ngroup, dimension);	/* voxel    counter */
	imgN = calloc_float2  (ngroup, vdim);		/* input N image */

	for (k = 0; k < ngroup; k++) {
/****************/
/* read N image */
/****************/
		printf ("Reading: %s\n", Nmgroot[k]);
		sprintf (imgfile, "%s.4dfp.img", Nmgroot[k]);
		if (!(fp = fopen (imgfile, "rb")) || eread (imgN[k], vdim, isbigN[k], fp)
		|| fclose (fp)) errr (program, imgfile);

/**********************************/
/* read and accumulate image data */
/**********************************/
		for (j = 0; j < nimg[k]; j++) {
			printf ("Reading: %s\n", imgroot[k][j]);
			sprintf (imgfile, "%s.4dfp.img", imgroot[k][j]);
			if (!(fp = fopen (imgfile, "rb")) || eread (img1, dimension, isbig[k][j], fp)
			|| fclose (fp)) errr (program, imgfile);
			for (i = 0; i < dimension; i++) {
				defined = isnormal (img1[i]) && img1[i] != (float) 1.e-37;
				if (!defined) continue;
				imgm[k][i]++;
				imgs[k][i] += img1[i];
				imgv[k][i] += img1[i]*img1[i];
			}
		}
		for (i = 0; i < dimension; i++) {
			u = (imgs[k][i] /= imgm[k][i]);
			imgv[k][i] -= u*u*imgm[k][i]; if (imgv[k][i] < 0.0) imgv[k][i] = 0.0;
			imgv[k][i] /= (imgm[k][i] - 1);
		}
/***************************************/
/* divide variance estimate by N image */
/***************************************/
		for (i = 0; i < vdim; i++) if (imgN[k][i]) for (j = 0; j < imgdim[3]; j++) {
			imgv[k][j*vdim + i] /= imgN[k][i];
		}
	}

/***************************************************************************/
/* compute voxelwize t statistic (Welch's approximate t' when ngroup == 2) */
/***************************************************************************/
	for (i = 0; i < dimension; i++) {
		if (ngroup == 2) {
			imgs[0][i] -= imgs[1][i];
			imgv[0][i] += imgv[1][i];
		}
		if (isnormal (imgs[0][i]) && isnormal (imgv[0][i])) {
			img1[i] = imgs[0][i] / sqrt (imgv[0][i]);
			if (img1[i] > amax) amax = img1[i];
			if (img1[i] < amin) amin = img1[i];
		} else switch (NaN_flag) {
			case 'Z': img1[i] = 0.0;		break;
			case 'E': img1[i] = (float) 1.e-37;	break;
			case 'N': img1[i] = rnan ();		break;
		}
	}

	if (!strlen (outroot)) {
/*******************/
/* compute outroot */
/*******************/
		strcpy (command, lstfile[0]);
		if (!(ptr = strrchr (command, '/'))) ptr = command; else ptr++;
		if (str = strrchr (ptr, '.')) *str = '\0';
		sprintf (outroot, "RFX_%s", ptr);
		if (ngroup == 2) {
			strcat (outroot, "_vs_");
			strcpy (command, lstfile[1]);
			if (!(ptr = strrchr (command, '/'))) ptr = command; else ptr++;
			if (str = strrchr (ptr, '.')) *str = '\0';
			strcat (outroot, ptr);
		}
		strcat (outroot, "_tmap");
	}
	sprintf (outfile, "%s.4dfp.img", outroot);

/****************/
/* write result */
/****************/
	printf ("Writing: %s\n", outfile);
	if (!(fp = fopen (outfile, "wb"))
	|| ewrite (img1, dimension, control, fp)
	|| fclose (fp)) errw (program, outfile);

/*******/
/* ifh */
/*******/
	Writeifh (program, outfile, &ifh, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s -r%.0fto%.0f", outroot, amin, amax);
	printf ("%s\n", command);
	status |= system (command);

/*******/
/* rec */
/*******/
	startrecle (outfile, argc, argv, rcsid, control);
	for (k = 0; k < ngroup; k++) {
		sprintf (command, "imglist %d\n", k + 1);
		printrec (command); catrec (lstfile[k]); printrec ("endimglist\n");
		if (saverec) {
			sprintf (imgfile, "%s.4dfp.img", Nmgroot[k]); catrec (imgfile);
			for (j = 0; j < nimg[k]; j++) {
				sprintf (imgfile, "%s.4dfp.img", imgroot[k][j]); catrec (imgfile);
			}
		}
	}
	endrec ();
	
	free (img1); free_double2 (imgs); free_double2 (imgv); free_int2 (imgm); free_float2 (imgN);
	exit (status);
}
