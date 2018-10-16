/*$Header: /data/petsun4/data1/src_solaris/algebra_4dfp/RCS/imgopr_4dfp.c,v 1.14 2010/12/23 02:43:21 avi Exp $*/
/*$Log: imgopr_4dfp.c,v $
 * Revision 1.14  2010/12/23  02:43:21  avi
 * float -> double imgo
 * UCHAR -> UINT imgc imgn
 *
 * Revision 1.13  2010/05/18  04:52:17  avi
 * operation -P (unsplit ROIs)
 *
 * Revision 1.12  2010/01/31  07:33:26  avi
 * extra space in command line buffer
 *
 * Revision 1.11  2009/07/17  06:09:52  avi
 * option -R
 *
 * Revision 1.10  2008/06/05  06:15:01  avi
 * enable binary operations in which single volume second input is paired with multi-volume first input
 * clarify operation specification in usage
 *
 * Revision 1.9  2007/07/05  04:37:00  avi
 * correct error in count_defined code (1.e37 -> (float) 1.e-37)
 * Linux compliant
 *
 * Revision 1.8  2006/09/24  02:55:52  avi
 * Solaris 10
 *
 * Revision 1.7  2006/09/13  02:09:48  avi
 * MAXR -> 4096
 *
 * Revision 1.6  2005/09/10  05:50:00  avi
 * fix usage and several bugs
 *
 * Revision 1.5  2005/09/10  04:21:47  avi
 * count defined (options -d and -u)
 * output defined (options -N -Z -E)
 *
 * Revision 1.4  2005/01/07  22:42:16  avi
 * -x and -y options
 * remove references to IFH and ifh.h
 *
 * Revision 1.3  2004/09/21  20:32:26  rsachs
 * Installed 'setprog'.
 *
 * Revision 1.2  2004/01/04  23:50:49  avi
 * prevent overwriting input by output
 *
 * Revision 1.1  2004/01/02  05:18:16  avi
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
#define UCHAR		unsigned char
#define UINT		unsigned int

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
static char	rcsid[] = "$Id: imgopr_4dfp.c,v 1.14 2010/12/23 02:43:21 avi Exp $";

void usage () {
	fprintf (stderr, "Usage:	%s -<operation><(4dfp) outroot> <(4dfp) image1> <(4dfp) image2> ...\n", program);
	fprintf (stderr, "	operation\n");
	fprintf (stderr, "	a	add\n");
	fprintf (stderr, "	s	subtract (image1 - image2)\n");
	fprintf (stderr, "	p	product\n");
	fprintf (stderr, "	r	ratio (image1 / image2)\n");
	fprintf (stderr, "	e	mean (expectation)\n");
	fprintf (stderr, "	v	variance\n");
	fprintf (stderr, "	g	geometric mean\n");
	fprintf (stderr, "	n	count defined (see -u option) voxels\n");
	fprintf (stderr, "	x	voxelwize maximum\n");
	fprintf (stderr, "	y	voxelwize minimum\n");
	fprintf (stderr, "	G	report serial number (counting from 1) of image with greatest value\n");
	fprintf (stderr, "	P	unsplit multiple ROIs into fidl compatible ROI file\n");
	fprintf (stderr, "	option\n");
	fprintf (stderr, "	-u	count only defined (not NaN or 1.e-37 or 0.0) voxels\n");
	fprintf (stderr, "	-R	suppress creation of rec file\n");
	fprintf (stderr, "	-N\toutput undefined voxels as NaN\n");
	fprintf (stderr, "	-Z\toutput undefined voxels as 0\n");
	fprintf (stderr, "	-E\toutput undefined voxels as 1.E-37 (default)\n");
	fprintf (stderr, "	-c<flt>	multiply output by specified scaling factor\n");
	fprintf (stderr, "	-l<lst>	read input file names from specified list file\n");
	fprintf (stderr, "	-@<b|l>\toutput big or little endian (default first input endian)\n");
	fprintf (stderr, "N.B.:	image dimensions must match except for binary operations {aspr} in which\n");
	fprintf (stderr, "	a 1 volume second image may be paired with a multi-volume first image\n");
	exit (1);
}

int main (int argc, char *argv[]) {
	FILE		*fp;
	IFH		ifh;
	char		imgroot[MAXR][MAXL], imgfile[MAXL], lstfile[MAXL] = "";
	char		outroot[MAXL] = "", outfile[MAXL];
	float		*img1;
	double		*imgs, *imgo;
	UINT		*imgn;		/* count defined voxels */
	UINT		*imgc;
	float		voxdim[3], voxdim1[3];
	float		sfactor = 1.0, amax = -FLT_MAX, amin = FLT_MAX;
	int		*nvox, nrun = 0, opr = 0;			/* operation */
	int		orient, orient1, imgdim[MAXR][4], dimension, vdim, isbig[MAXR];
	int		ndefined;
	char		control = '\0';

/***********/
/* utility */
/***********/
	double		q, u;
	int		c, i, j, k, m;
	char		*ptr, command[4*MAXL];			/* extra space for use in printrec */
	char		*srgv[MAXL];				/* list file string field pointers */

/*********/
/* flags */
/*********/
	int		debug = 0;
	int		count_defined = 0;
	int		NaN_flag = 'E';		/* 'E' 1.e-37; 'Z' 0.0; 'N' NaN; */
	int		bmode_flag = 0;	/* when true, first frame of imput 1 is applied to all frames of input 0 */
	int		status = 0;
	int		saverec = 1;

	printf ("%s\n", rcsid);
	setprog (program, argv);
	for (k = 0; k < MAXR; k++) imgroot[k][0] = '\0';
/************************/
/* process command line */
/************************/
	for (i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'R': saverec = 0;				break;
				case 'd': debug++;				break;
				case 'u': count_defined++;			break;
				case 'N': case 'Z': case 'E': NaN_flag = c;	break;
				case 'l': strcpy (lstfile, ptr);		*ptr =  '\0'; break;
				case 'P':
				case 'n': count_defined++;
				case 'a':
				case 's':
				case 'p':
				case 'r':
				case 'e':
				case 'v':
				case 'x':
				case 'y':
				case 'G':
				case 'g': opr = c; getroot (ptr, outroot);	*ptr =  '\0'; break;
				case 'c': sfactor = atof (ptr); 		*ptr =  '\0'; break;
				case '@': control = *ptr++;			*ptr = '\0'; break;
			}
		} else {
			getroot (argv[i], imgroot[nrun++]);
		}	
	}

/*******************/
/* parse list file */
/*******************/
	if (strlen (lstfile)) {
		if (!(fp = fopen (lstfile, "r"))) errr (program, lstfile);
		while (fgets (command, MAXL, fp)) {
			if (nrun >= MAXR) {
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
			getroot (srgv[0], imgroot[nrun++]);
		}
		fclose (fp);
	}
	if (!nrun || !strlen (outroot)) usage ();

/*********************************/
/* check dimensional consistency */
/*********************************/
	for (i = 0; i < nrun; i++) {
		printf ("input %3d: %s\n", i + 1, imgroot[i]);
		if (!strcmp (imgroot[i], outroot)) {
			fprintf (stderr, "%s: output %s matches one or more inputs\n", program, outroot);
			exit (-1);
		}
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
			fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot[0], imgroot[i]);
			exit (-1);
		}
	}
	bmode_flag = !status && strchr ("aspr", opr) && nrun == 2 && imgdim[0][3] > 1 && imgdim[1][3] == 1;
	if (!bmode_flag) for (i = 0; i < nrun; i++) {
		status |= (imgdim[i][3] != imgdim[0][3]);
		if (status) {
			fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot[0], imgroot[i]);
			exit (-1);
		}
	}
	if (!control) control = (isbig[0]) ? 'b' : 'l';

	vdim = imgdim[0][0]*imgdim[0][1]*imgdim[0][2];
	dimension = vdim*imgdim[0][3];
	if (!(img1 = (float *)  calloc (dimension, sizeof (float))))	errm (program);
	if (!(imgo = (double *) calloc (dimension, sizeof (double))))	errm (program);
	if (!(imgs = (double *) calloc (dimension, sizeof (double))))	errm (program);
	if (!(imgc = (UINT *)   calloc (dimension, sizeof (UINT))))	errm (program);
	if (!(imgn = (UINT *)   calloc (dimension, sizeof (UINT))))	errm (program);
	if (!(nvox = (int *)    calloc (nrun,      sizeof (int))))	errm (program);

/******************/
/* prepare arrays */
/******************/
	for (j = 0; j < dimension; j++) {
		switch (opr) {
			case 'p':			/* product */
			case 'g':			/* geometric mean */
				imgo[j] = 1.0;
				break;
			case 'G':			/* Tony Jack operation */
			case 'x':			/* voxelwize maximum */
				imgo[j] = -FLT_MAX;
				break;
			case 'y':			/* voxelwize minimum */
				imgo[j] =  FLT_MAX;
				break;
			default:
				break;
		}
	}

/******************/
/* start rec file */
/******************/
if (saverec) {
	sprintf (outfile, "%s.4dfp.img", outroot);
	startrecle (outfile, argc, argv, rcsid, control);
	if (strlen (lstfile)) {
		printrec ("imglist\n"); catrec (lstfile); printrec ("endimglist\n");
	}
	if (bmode_flag) {
		sprintf (command, "first frame of %s applied to all frames of %s\n", imgroot[1], imgroot[0]);
		printrec (command);
	}
}

/*******************/
/* execute algebra */
/*******************/
	for (i = 0; i < nrun; i++) {
		sprintf (imgfile, "%s.4dfp.img", imgroot[i]);
		printf ("Reading: %s\n", imgfile);
		j = imgdim[i][0]*imgdim[i][1]*imgdim[i][2]*imgdim[i][3];
		if (!(fp = fopen (imgfile, "rb")) || eread (img1, j, isbig[i], fp)
		|| fclose (fp)) errr (program, imgfile);
		if (bmode_flag && i == 1) for (k = 1; k < imgdim[0][3]; k++) for (j = 0; j < vdim; j++) {
			(img1 + k*vdim)[j] = img1[j];	/* duplicate first volume */
		}
		if (saverec) catrec (imgfile);

		for (j = 0; j < dimension; j++) {
			q = img1[j];
			if (count_defined && (img1[j] == 0.0 || img1[j] == (float) 1.e-37 || !isnormal (q))) continue;
			imgn[j]++;
			switch (opr) {
				case 'a':
				case 'e':
					imgo[j] += img1[j];				break;
				case 's':
					switch (i) {
						case 0: imgo[j]  = img1[j]; break;
						case 1: imgo[j] -= img1[j]; break;
					}						break;
				case 'p':
				case 'g':
					imgo[j] *= img1[j];				break;
				case 'n':
					imgo[j] =  imgn[j];				break;
				case 'r':
					switch (i) {
						case 0: imgo[j]  = img1[j]; break;
						case 1: imgo[j] /= img1[j]; break;
					}						break;
				case 'v':
					imgo[j] += img1[j];
					imgs[j]	+= img1[j]*img1[j];			break;
				case 'x':
					if (img1[j] > imgo[j]) imgo[j] = img1[j];	break;
				case 'y':
					if (img1[j] < imgo[j]) imgo[j] = img1[j];	break;
				case 'G':
					if (img1[j] > imgo[j]) {
						imgo[j] = img1[j];
						imgc[j] = i + 1;
					}						break;
				case 'P':
					imgo[j] = i + 2; nvox[i]++;			break;
				default:						break;
			}
		}
	}

	for (ndefined = j = 0; j < dimension; j++) {
		q  = 1.0 / imgn[j];
		switch (opr) {
			case 'e':
				imgo[j] *= q;
				break;
			case 'g':
				imgo[j] = pow (imgo[j], q);
				break;
			case 'v':
				u = imgo[j] * q;
				imgo[j]	= (imgs[j] - u*u/q) / (imgn[j] - 1);
				if (imgn[j] > 1 && imgo[j] < 0.0) imgo[j] = 0.0;
				break;
			case 'G':
				imgo[j] = imgc[j];
				break;
			default:
				break;
		}
		u = imgo[j];
		if (isnormal (u) || u == 0.0) {
			img1[j] = imgo[j]*sfactor;
			if (img1[j] > amax) amax = img1[j];
			if (img1[j] < amin) amin = img1[j];
			ndefined++;
		} else switch (NaN_flag) {
			case 'Z': img1[j] = 0.0;		break;
			case 'E': img1[j] = (float) 1.e-37;	break;
			case 'N': img1[j] = rnan ();		break;
		}
	}

/****************/
/* write result */
/****************/
	sprintf (outfile, "%s.4dfp.img", outroot);
	printf ("Writing: %s\n", outfile);
	if (!(fp = fopen (outfile, "wb"))
	|| ewrite (img1, dimension, control, fp)
	|| fclose (fp)) errw (program, outfile);

/*******/
/* ifh */
/*******/
	Writeifh (program, outfile, &ifh, control);
	if (opr == 'P') {
		sprintf (outfile, "%s.4dfp.ifh", outroot);
		if (!(fp = fopen (outfile, "a"))) errw (program, outfile);
		for (i = 0; i < nrun; i++) {
			fprintf (fp, "region names\t:= %-5d\t%s\t%10d\n", i, imgroot[i], nvox[i]);
		}
		if (fclose (fp)) errw (program, outfile);
	}

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s -r%.0fto%.0f", outroot, amin, amax);
	printf ("%s\n", command);
	status |= system (command);

/*******/
/* rec */
/*******/
if (saverec) {
	sprintf (command, "defined or zero output voxel count = %d out of %d total\n", ndefined, dimension);
	printrec (command);
	endrec ();
}

	free (img1); free (imgo); free (imgs); free (imgc); free (imgn); free (nvox);
	exit (status);
}
