/*$Header: /data/petsun4/data1/src_solaris/algebra_4dfp/RCS/ROI_resolve_4dfp.c,v 1.5 2012/04/18 02:05:51 avi Exp $*/
/*$Log: ROI_resolve_4dfp.c,v $
 * Revision 1.5  2012/04/18  02:05:51  avi
 * MAXR -> 1000
 *
 * Revision 1.4  2009/10/12  04:14:49  avi
 * correct rec code
 *
 * Revision 1.3  2009/09/08  23:50:08  avi
 * Solaris 10 and linux compliant
 * new (better) algorithm based on imgw and imgr
 *
 * Revision 1.2  2005/08/05  05:48:11  avi
 * compute virtual coordinates correctly using ifh center and mmppix
 *
 * Revision 1.1  2005/08/04  07:02:22  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>		/* R_OK */
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MAXL		256	/* maximum string length */
#define MAXR		1000	/* maximum number of ROIs */

void vrtflip  (int *iori, int *imgdim, float *centeri, float *mmppixi, float *centerr, float *mmppixr);	/* cvrtflip.c */

/********************/
/* global variables */
/********************/
char		program[MAXL];

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

typedef struct {
	char	imgroot[MAXL];
	float	roicom[3];
	int	nvox, nvox1, isbig;
} ROI_INFO;

void usage () {
	fprintf (stderr, "Usage:	%s <(4dfp) ROI1> <(4dfp) ROI2> <(4dfp) ROI3> ...\n", program);
	fprintf (stderr, "	option\n");
	fprintf (stderr, "	-l<lst>	read input file names from specified list file\n");
	fprintf (stderr, "	-@<b|l>\toutput big or little endian (default CPU endian)\n");
	fprintf (stderr, "N.B.:	output 4dfp fileroots are same as inputs with appended \"z\"\n");
	exit (1);
}

static char rcsid[] = "$Id: ROI_resolve_4dfp.c,v 1.5 2012/04/18 02:05:51 avi Exp $";
int main (int argc, char *argv[]) {
	FILE		*fp;
	char		imgfile[MAXL], ifhfile[MAXL], lstfile[MAXL] = "";
	char		outroot[MAXL], outfile[MAXL], tmpfile[MAXL];
	IFH		ifh, ifh1;
	float		mmppixr[3], centerr[3];
	float		**imag, *imgr;
	int		ix, iy, iz, jndex, dimension;
	int		*imgw;		/* tracks winning ROI at each voxel */
	char		control = '\0';

/**********************/
/* roi center of mass */
/**********************/
	int		nroi = 0;
	ROI_INFO	roi_info[MAXR];
	double		fndex[3], r1[3], rsq1;

/***********/
/* utility */
/***********/
	int		c, i, j, k, m;
	char		*ptr, command[MAXL];
	char		*srgv[MAXL];				/* list file string field pointers */

/*********/
/* flags */
/*********/
	int		debug = 0;
	int		status = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);
	for (k = 0; k < MAXR; k++) roi_info[k].imgroot[0] = '\0';

/************************/
/* process command line */
/************************/
	for (i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'd': debug++;				break;
				case 'l': strcpy (lstfile, ptr);		*ptr = '\0'; break;
				case '@': control = *ptr++;			*ptr = '\0'; break;
			}
		} else {
			getroot (argv[i], roi_info[nroi++].imgroot);
		}	
	}

/*******************/
/* parse list file */
/*******************/
	if (strlen (lstfile)) {
		if (!(fp = fopen (lstfile, "r"))) errr (program, lstfile);
		while (fgets (command, MAXL, fp)) {
			if (nroi >= MAXR) {
				fprintf (stderr, "%s: maximum roi count (%d) exceeded\n", program, MAXR);
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
			getroot (srgv[0], roi_info[nroi++].imgroot);
		}
		fclose (fp);
	}
	if (!nroi) usage ();

/*********************************/
/* check dimensional consistency */
/*********************************/
	for (i = 0; i < nroi; i++) {
		if (!strcmp (roi_info[i].imgroot, outroot)) {
			fprintf (stderr, "%s: output %s matches one or more inputs\n", program, outroot);
			exit (-1);
		}
		sprintf (ifhfile, "%s.4dfp.ifh", roi_info[i].imgroot);
		printf ("Reading: %s\n", ifhfile);
		if (!i) if (Getifh (ifhfile, &ifh))  errr (program, ifhfile);
		if (Getifh (ifhfile, &ifh1)) errr (program, ifhfile);
		roi_info[i].isbig = strcmp (ifh1.imagedata_byte_order, "littleendian");
		status = ifh1.orientation - ifh.orientation;
		for (k = 0; k < 3; k++) {
			status |= ifh1.matrix_size[k] - ifh.matrix_size[k];
			status |= (fabs ((double) (ifh1.mmppix[k] - ifh.mmppix[k])) > 1.e-5);
		}
		if (status) {
			fprintf (stderr, "%s: %s %s dimension mismatch\n",
						program, roi_info[0].imgroot, roi_info[i].imgroot);
			exit (-1);
		}
	}

/*****************************************/
/* virtual flip instead of x4dfp2ecat () */
/*****************************************/
	vrtflip (&ifh.orientation, ifh.matrix_size, ifh.center, ifh.mmppix, centerr, mmppixr);
	if (debug) {
		printf ("orient=%d\n", ifh.orientation);
		printf ("image dimensions         %10d%10d%10d\n",
						ifh.matrix_size[0], ifh.matrix_size[1], ifh.matrix_size[2]);
		printf ("virtually flipped mmppix %10.6f%10.6f%10.6f\n", mmppixr[0], mmppixr[1], mmppixr[2]);
		printf ("virtually flipped center %10.4f%10.4f%10.4f\n", centerr[0], centerr[1], centerr[2]);
	}
	dimension = 1; for (k = 0; k < 3; k++) dimension *= ifh.matrix_size[k];
	imag = calloc_float2 (nroi, dimension);
	if (!(imgw = calloc (dimension, sizeof (int))))		errm (program);
	if (!(imgr = calloc (dimension, sizeof (float))))	errm (program);

/**************************************/
/* compute center of mass for all ROI */
/**************************************/
	for (i = 0; i < nroi; i++) {
		sprintf (imgfile, "%s.4dfp.img", roi_info[i].imgroot);
		printf ("Reading: %s\n", imgfile);
		if (!(fp = fopen (imgfile, "rb")) || eread (imag[i], dimension, roi_info[i].isbig, fp)
		|| fclose (fp)) errr (program, imgfile);

		roi_info[i].roicom[0] = roi_info[i].roicom[1] = roi_info[i].roicom[2] = roi_info[i].nvox = jndex = 0;
		for (iz = 0; iz < ifh.matrix_size[2]; iz++) {fndex[2] = (double) (iz + 1);
		for (iy = 0; iy < ifh.matrix_size[1]; iy++) {fndex[1] = (double) (iy + 1);
		for (ix = 0; ix < ifh.matrix_size[0]; ix++) {fndex[0] = (double) (ix + 1);
			if (imag[i][jndex++] == 0.0) continue;
			for (k = 0; k < 3; k++) roi_info[i].roicom[k] += fndex[k]*mmppixr[k] - centerr[k];
			roi_info[i].nvox++;
		}}}
		for (k = 0; k < 3; k++) roi_info[i].roicom[k] /= roi_info[i].nvox;
	}

/*****************************************************/
/* compute winning ROI at each voxel and zero losers */
/*****************************************************/
	jndex = 0;
	for (iz = 0; iz < ifh.matrix_size[2]; iz++) {fndex[2] = (double) (iz + 1);
	for (iy = 0; iy < ifh.matrix_size[1]; iy++) {fndex[1] = (double) (iy + 1);
	for (ix = 0; ix < ifh.matrix_size[0]; ix++) {fndex[0] = (double) (ix + 1);
		for (i = 0; i < nroi; i++) if (imag[i][jndex] > 0.) {
			for (rsq1 = k = 0; k < 3; k++) {
				r1[k] = fndex[k]*mmppixr[k] - centerr[k] - roi_info[i].roicom[k];
				rsq1 += r1[k]*r1[k];
			}
			if (!imgw[jndex] || rsq1 < imgr[jndex]) {
				imgr[jndex] = rsq1;
				imgw[jndex] = i + 1;
			}
		}
		for (i = 0; i < nroi; i++) if (imgw[jndex] != i + 1) imag[i][jndex] = 0.;
		jndex++;
	}}}

/*****************/
/* write results */
/*****************/
	for (i = 0; i < nroi; i++) {
		roi_info[i].nvox1 = jndex = 0;
		for (iz = 0; iz < ifh.matrix_size[2]; iz++) {
		for (iy = 0; iy < ifh.matrix_size[1]; iy++) {
		for (ix = 0; ix < ifh.matrix_size[0]; ix++) {
			if (imag[i][jndex++] == 0.0) continue;
			roi_info[i].nvox1++;
		}}}
		sprintf (outroot, "%sz", roi_info[i].imgroot);
		sprintf (outfile, "%s.4dfp.img", outroot);
		printf ("Writing: %s\n", outfile);
		if (!(fp = fopen (outfile, "wb")) || ewrite (imag[i], dimension, control, fp)
		|| fclose (fp)) errr (program, outfile);
/*******/
/* ifh */
/*******/
		Writeifh (program, outfile, &ifh, control);

/*******/
/* hdr */
/*******/
		sprintf (command, "ifh2hdr %s", outroot);
		printf ("%s\n", command);
		status |= system (command);
/*******/
/* rec */
/*******/
		startrece (outfile, argc, argv, rcsid, control);
		sprintf (command, "center of mass (given ifh center and mmppix)%10.4f%10.4f%10.4f\n",
			roi_info[i].roicom[0], roi_info[i].roicom[1], roi_info[i].roicom[2]);
				printrec (command);
		sprintf (command, "voxel count %d -> %d\n", roi_info[i].nvox, roi_info[i].nvox1);
				printrec (command);
		sprintf (tmpfile, "%s.4dfp.img.rec", roi_info[i].imgroot);
		catrec (tmpfile);
		endrec ();
	}

	free (imgw); free (imgr);
	free_float2 (imag);
	exit (status);
}
