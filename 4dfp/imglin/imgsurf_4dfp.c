/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/imgsurf_4dfp.c,v 1.4 2010/08/11 07:20:39 avi Exp $*/
/*$Log: imgsurf_4dfp.c,v $
 * Revision 1.4  2010/08/11  07:20:39  avi
 * eliminate NPTMAX limit on number of input points
 *
 * Revision 1.3  2009/08/09  05:46:14  avi
 * make output more useful for ECoG analysis
 *
 * Revision 1.2  2007/09/06  03:39:01  avi
 * strip off path and extensions in constructing outfile
 *
 * Revision 1.1  2007/09/05  04:05:27  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>	/* R_OK */
#include <endianio.h>
#include <Getifh.h>

#define MAXL	256

/*************/
/* externals */
/*************/
extern void	vrtflip_ (int *orientation, int *imgdim, float *center, float *mmppix, float *centerr, float *mmppixr);
extern void	imgvalx_ (float *imgt, int *nx, int *ny, int *nz, float *center, float *mmppix, float *x, float *v, int *lslice);
extern void	fimgsurf_ (float *imgt, int *nx, int *ny, int *nz, float *center, float *mmppix, float *x0, float *grad, float *xopt);

/***********/
/* globals */
/***********/
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

void errx (char* program, char* filespc) {
	fprintf (stderr, "%s: %s format error\n", program, filespc);
	exit (-1);
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

static char rcsid [] = "$Id: imgsurf_4dfp.c,v 1.4 2010/08/11 07:20:39 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* image i/o */
/*************/
	FILE		*outfp, *pntfp, *imgfp;
	IFH		ifh;
	char		imgroot[MAXL], imgfile[MAXL], pntfile[MAXL]; 
	char		outroot[MAXL], outfile[MAXL];
	int		isbig;

/***********/
/* utility */
/***********/
	char		*ptr, *ptr1, string[MAXL], program[MAXL], command[MAXL], *srgv[MAXL];
	int		c, i, j, k, m;
	double		q, t;

/**************/
/* image data */
/**************/
	float		*imgt;
	float		mmppixr[3], centerr[3];		/* as in peak_4dfp.c */
	float		**x0, **xopt, **grad;
	int		imgdim[4], vdim, nx, ny, nz, ipt, npt;

/*********/
/* flags */
/*********/
	int		status = 0;

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);
/******************************/
/* get command line arguments */
/******************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
			}
		} else switch (k) {
		 	case 0: getroot (argv[i], imgroot);	k++; break;
		 	case 1: strcpy (pntfile, argv[i]);	k++; break;
		}	
	}
	if (k < 2) {
		printf ("Usage:	%s <(4dfp) image> <point_list>\n", program);
		printf ("\toption\n");
		printf ("N.B.:\t<point_list> lists loci in atlas coordinates (X Y Z) in mm\n");
		exit (1);
	}

/**************/
/* point list */
/**************/
	printf ("Reading: %s\n", pntfile);
	if (!(pntfp = fopen (pntfile, "r"))) errr (program, pntfile);
	npt = 0; while (fgets (string, MAXL, pntfp)) {
		m = split (string, srgv, MAXL);
		if (m && m != 3) errx (program, pntfile);
		npt++;
	}
	x0	= calloc_float2 (npt, 3);
	grad	= calloc_float2 (npt, 3);
	xopt	= calloc_float2 (npt, 3);
	rewind (pntfp);
	for (ipt = 0; ipt < npt; ipt++) {
		fgets (string, MAXL, pntfp);
		m = split (string, srgv, MAXL);
		if (m) for (k = 0; k < 3; k++) x0[ipt][k] = atof (srgv[k]);
	}
	fclose (pntfp);

/**************/
/* read image */
/**************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (Getifh (imgroot, &ifh)) errr (program, imgfile);
	for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];
	isbig = strcmp (ifh.imagedata_byte_order, "littleendian");
	nx = imgdim[0];
	ny = imgdim[1];
	nz = imgdim[2];
	vdim = imgdim[0]*imgdim[1]*imgdim[2];
	if (!(imgt = (float *) malloc (vdim * sizeof (float)))) errm (program);
	if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
	printf ("Reading: %s\n", imgfile);
	if (eread (imgt, vdim, isbig, imgfp)) errr (program, imgfile);
	fclose (imgfp);

/*****************************************/
/* virtual flip instead of x4dfp2ecat () */
/*****************************************/
	vrtflip_ (&ifh.orientation, ifh.matrix_size, ifh.center, ifh.mmppix, centerr, mmppixr);
	printf ("atlas mmppix     \t%10.6f%10.6f%10.6f\n", mmppixr[0], mmppixr[1], mmppixr[2]); 
	printf ("atlas center     \t%10.4f%10.4f%10.4f\n", centerr[0], centerr[1], centerr[2]);

	printf ("Calling fimgsurf on %d points\n", npt);
	for (ipt = 0; ipt < npt; ipt++) {
		fimgsurf_ (imgt, &nx, &ny, &nz, centerr, mmppixr, x0[ipt], grad[ipt], xopt[ipt]);
	}

/*****************************/
/* construct output filename */
/*****************************/
	if (!(ptr = strrchr (imgroot, '/'))) ptr = imgroot; else ptr++;
	strcpy (string, pntfile); if ((ptr1 = strrchr (string, '.'))) *ptr1 = '\0';
	if (!(ptr1 = strrchr (string, '/'))) ptr1 = string; else ptr1++;
	sprintf (outfile, "%s_%s_%s.dat", ptr, ptr1, program);

/*************************/
/* write output dat file */
/*************************/
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);
	fprintf (outfp, "%10s%10s%10s%10s%10s%10s%10s%10s\n", "point", "x0", "y0", "z0", "xsurf", "ysurf", "zsurf", "distance");
	for (i = 0; i < npt; i++) {
		for (q = k = 0; k < 3; k++) {
			t = x0[i][k] - xopt[i][k];
			q += t*t;
		}
		fprintf (outfp, "%10d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n", i + 1,
			x0[i][0],   x0[i][1],   x0[i][2],
			xopt[i][0], xopt[i][1], xopt[i][2],
			sqrt (q));
	}
	if (fclose (outfp)) errw (program, outfile);

/***********************************************/
/* write output modified input coordinate list */
/***********************************************/
	if (!(ptr = strrchr (pntfile, '/'))) ptr = pntfile; else ptr++;
	sprintf (outfile, "%s.surf", ptr);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);
	for (i = 0; i < npt; i++) {
		fprintf (outfp, "%10.4f%10.4f%10.4f\n", xopt[i][0], xopt[i][1], xopt[i][2]);
	}
	if (fclose (outfp)) errw (program, outfile);

	free (imgt);
	free_float2 (x0); free_float2 (xopt); free_float2 (grad); 
        exit (status);
}
