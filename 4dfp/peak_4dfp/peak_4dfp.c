/*$Header: /data/petsun4/data1/src_solaris/peak_4dfp/RCS/peak_4dfp.c,v 1.25 2012/03/30 23:11:05 avi Exp $*/
/*$Log: peak_4dfp.c,v $
 * Revision 1.25  2012/03/30  23:11:05  avi
 * correct bug in computation of del2v
 *
 * Revision 1.24  2011/04/08  05:05:01  avi
 * option -N
 *
 * Revision 1.23  2010/08/07  23:33:13  avi
 * option -F
 *
 * Revision 1.22  2010/03/24  03:42:31  avi
 * option -a
 *
 * Revision 1.21  2007/09/23  04:02:14  avi
 * eliminate typo
 *
 * Revision 1.20  2007/09/23  03:55:56  avi
 * linux compliant (eliminate #include <sunmath.h>)
 *
 * Revision 1.19  2007/06/04  04:15:25  avi
 * correct failure to set isbigm
 *
 * Revision 1.18  2006/09/26  01:46:55  avi
 * #include <sunmath.h> and prototype imgvalx_()
 *
 * Revision 1.17  2006/09/26  01:17:05  avi
 * scrupulously correct qsort usage
 *
 * Revision 1.16  2006/09/25  21:35:56  avi
 * nominally Solaris 10
 *
 * Revision 1.15  2006/09/07  06:47:01  avi
 * correctly output rec after ifh and hdr
 *
 * Revision 1.14  2006/09/07  06:35:18  avi
 * -q (suppress rec listing)
 *
 * Revision 1.13  2004/01/14  05:51:31  avi
 * NTOP increased to 1000
 * correct logic for deciding to call catrec () on mskfile
 *
 * Revision 1.12  2002/05/02  02:39:04  avi
 * fix bug in voxels per region accumulator
 *
 * Revision 1.11  2002/03/25  23:33:11  avi
 * make mask values biolar
 *
 * Revision 1.10  2002/02/08  07:36:31  avi
 * -m option (mask ROI output)
 * restore -O compilation as likely memory management error corrected
 *
 * Revision 1.9  2002/02/08  05:45:16  avi
 * make preblur an option
 * remove -f option
 *
 * Revision 1.8  2002/02/08  00:25:55  avi
 * -v (value thresholds in constructing initial list)
 * keep track of nvox in each ROI
 *
 * Revision 1.7  2002/01/31  06:33:34  avi
 * -o (output ROI file)
 *
 * Revision 1.6  2002/01/31  02:52:49  avi
 * consolidate peaks closer than dthresh
 *
 * Revision 1.5  2002/01/29  03:21:04  avi
 * close pair (pdist < dthresh) listing
 * compile without -O as otherwise core dumps keep happening
 *
 * Revision 1.4  2002/01/29  01:32:26  avi
 * initial malloc extremum memory to prevent unpredictable BusErrors
 * return necessary 0 in pcompare if (p1->v == p2->v)
 * call imgvalx only after extremum identified
 * correct computation of dvdx[] by including factor of 0.5
 *
 * Revision 1.3  2002/01/27  23:34:01  avi
 * -n option
 * -c option
 * srad == 0.0 causes smoothing logic to be skipped
 *
 * Revision 1.2  2001/07/11  06:11:45  avi
 * use hsphere_4dfp for convolution
 *
 * Revision 1.1  2001/07/09  21:26:52  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>	/* R_OK */
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MSIZE		2048
#define NTOP		1000
#define MAXL		256

typedef struct {
	float	x[3];
	float	v;
	float	del2v;
	float	weight;
	int	nvox;
	int	killed;
} EXTREMUM;

/*************/
/* externals */
/*************/
void vrtflip_ 	(int *orientation, int *imgdim, float *center, float *mmppix, float *centerr, float *mmppixr);
void imgvalx_	(float *imgt, int *nx, int *ny, int *nz, float *center, float *mmppix, float *x, float *v, int *lslice);
void index_flip (int orientation, int *imgdim, float *fndex);	/* below */

static int pcompare (const void *ptr1, const void *ptr2) {
	EXTREMUM	*p1, *p2;
	
	p1 = (EXTREMUM *) ptr1; p2 = (EXTREMUM *) ptr2;
	if (p1->v == p2->v) return 0;	/* necessary for correct function */
	return (fabs (p1->v) > fabs (p2->v)) ? -1 : 1;
}

static double pdist (EXTREMUM *p1, EXTREMUM *p2) {
	int		k;
	double		q;

	for (q = k = 0; k < 3; k++) q += (p2->x[k] - p1->x[k]) * (p2->x[k] - p1->x[k]);
	return sqrt (q);
}

void getrange (char *string, float *minval, float *maxval) {
        char	*str;

	str = strstr (string, "to");
	if (str) {
		*str = '\0';
		*minval = atof (string);
		*maxval = atof (str + 2);
	} else {
		*minval = 0.0;
		*maxval = atof (string);
	}
}

/*************************************************************/
/* replacement for FORTRAN-style intrinsic in Sun libsunmath */
/*************************************************************/
int nint (float x) {
	int	i;
	float	y;

	y = x + 0.5;
	i = (int) y;
	if (y < 0.) i--;
	return i;
}

/********************/
/* global variables */
/********************/
static	IFH		imgifh, mskifh;
static	float		mmppixr[3], centerr[3];
static	float		dthresh = 0.0;
static	int		ntop = NTOP;
static	int		ntot = 0;

static void consolidate (EXTREMUM *plst, int nlst) {
	float		x[3], t, dmin;
	int		i, j, imin, jmin, k;
	int		npair;

	do {	printf ("peak pairs closer than %.4f mm\n", dthresh);
		dmin = 1.e6;
		npair = 0;
		for (i = 0; i < ntop && i < nlst; i++)		{if (plst[i].killed) continue;
		for (j = i + 1; j < ntop && j < nlst; j++)	{if (plst[j].killed) continue;
			t = (float) pdist (plst + i, plst + j);
			if (t < dthresh) {
				printf ("%5d%5d%10.4f\n", i + 1, j + 1, t);
				if (t < dmin) {
					imin = i;
					jmin = j;
					dmin = t;
				}
				npair++;
			}
		}}
		printf ("npair = %d\n", npair);
		if (npair) {
			for (k = 0; k < 3; k++) {
				x[k] = plst[imin].weight*plst[imin].x[k] + plst[jmin].weight*plst[jmin].x[k];
			}
			plst[imin].weight += plst[jmin].weight;
			for (k = 0; k < 3; k++) plst[imin].x[k] = x[k] / plst[imin].weight;
			plst[jmin].killed++;
			npair--;
		}
	} while (npair);
}

/********************************************/
/* peaklist increments global variable ntot */
/********************************************/
static void peaklist (FILE *fp, EXTREMUM *plst, int nlst, int nvox_flag) {
	int		i, k;
	float		fndex[3];

	fprintf (fp, "%-5s%10s%10s%10s%10s%10s%10s%10s%10s%",
		"ROI", "index_x", "index_y", "index_z", "atlas_x", "atlas_y", "atlas_z", "value", "curvature");
	if (nvox_flag) fprintf (fp, "%10s", "nvox");
	fprintf (fp, "\n");
	for (i = 0; i < ntop && i < nlst; i++) {
		if (plst[i].killed) continue;
		for (k = 0; k < 3; k++) fndex[k] = (plst[i].x[k] + centerr[k]) / mmppixr[k];
		index_flip (imgifh.orientation, imgifh.matrix_size, fndex);
		fprintf (fp, "%-5d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.6f", ntot++ + 1,
			fndex[0], fndex[1], fndex[2],
			plst[i].x[0], plst[i].x[1], plst[i].x[2], plst[i].v, -plst[i].del2v);
		if (nvox_flag) fprintf (fp, "%10d", plst[i].nvox);
		fprintf (fp, "\n");
	}
}

static char rcsid[] = "$Id: peak_4dfp.c,v 1.25 2012/03/30 23:11:05 avi Exp $";
int main (int argc, char *argv[]) {
	FILE		*imgfp;			/* input file */
	FILE		*outfp;			/* output file */

/*************/
/* image i/o */
/*************/
	char		*ptr, string[MAXL], program[MAXL], command[MAXL];
	char		imgroot[MAXL], imgfile[MAXL], ifhfile[MAXL]; 
	char		mskroot[MAXL], mskfile[MAXL]; 
	char		tmproot[MAXL], tmpfile[MAXL]; 
	char		trailer[MAXL] = "", outroot[MAXL], outfile[MAXL], recfile[MAXL];

/***************/
/* computation */
/***************/
	float		*imgv, *imgr, *imgm;
	float		v, frac, t, del2v, dmin;
	float		fndex[3], x[3], w[3], dvdx[3], d2vdx2[3];
	float		ctneg = 0.0, ctpos = 0.0;
	float		vtneg = 0.0, vtpos = 0.0;
	float		srad = 0.0, orad = 0;
	int		ix, iy, iz, dim, nx, ny, nz;
	int		c, i, j, k, l, jmin;
	int		isbig, isbigm;
	char		control = '\0';

/**************/
/* peak lists */
/**************/
	EXTREMUM	*ppos, *pneg, *pall, *pallN, ptmp[1];
	int		mpos = 0, mneg = 0, npos = 0, nneg = 0, nall, nallN, Nvoxcrit = 1;
	int		npos0, npos1, npos2;
	int		nneg0, nneg1, nneg2;

/*********/
/* flags */
/*********/
	int		mask = 0;
	int		status = 0;
	int		debug = 0;
	int		quiet = 0;
	int		forceblur = 0;
	
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
				case 'q': quiet++;				break;
				case 'F': forceblur++;				break;
				case 's': srad = atof (ptr);			*ptr = '\0'; break;
				case 'd': dthresh = atof (ptr);			*ptr = '\0'; break;
				case 'n': ntop = atoi (ptr);			*ptr = '\0'; break;
				case 'N': Nvoxcrit = atoi (ptr);		*ptr = '\0'; break;
				case 'o': orad = atof (ptr);			*ptr = '\0'; break;
				case '@': control = *ptr++;			*ptr = '\0'; break;
				case 'c': getrange (ptr, &ctneg, &ctpos);	*ptr = '\0'; break;
				case 'v': getrange (ptr, &vtneg, &vtpos);	*ptr = '\0'; break;
				case 'm': getroot (ptr, mskroot); mask++;	*ptr = '\0'; break;
				case 'a': strncpy (trailer, ptr, MAXL - 1);	*ptr = '\0'; break;
			}
		} else switch (k) {
		 	case 0: getroot (argv[i], imgroot);	k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage: %s <file_4dfp>\n", program);
		printf (" e.g., %s grand_average_222[.4dfp.img] -s10\n", program);
		printf ("\toption\n");
		printf ("\t-s<flt>\tpreblur with hard sphere kernel of specified radius\n");
		printf ("\t-n<int>\tlimit initial pos and neg peak list lengths (default=%d)\n", NTOP);
		printf ("\t-c<flt>[to<flt>] specify sign inverted curvature thresholds (default none)\n");
		printf ("\t-v<flt>[to<flt>] specify peak value thresholds (default none)\n");
		printf ("\t-d<flt>\tconsolidate extremum pairs closer than specified distance\n");
		printf ("\t-o<flt>\toutput a fidl compatible 4dfp format ROI file with regions of specified radius\n");
		printf ("\t-m<str>\tapply named mask file to output ROIs\n");
		printf ("\t-N<int>\tspecify output ROI minimum voxel count (default = 1)\n");
		printf ("\t-a<str>\tappend specified string to ROI output filename\n");
		printf ("\t-q\tquiet mode (suppress rec file listing)\n");
		printf ("\t-F\tforce preblur image creation even if hsphere_4dfp result exists\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:\toperations controlled by options -s, -n, -c, -v, -d, -o, -m, -N are applied serially in listed order\n");
		printf ("N.B.:\tall distances are in mm\n");
		printf ("N.B.:\toption -s<flt> creates a blurred image by invoking hsphere_4dfp\n");
		printf ("N.B.:\toption -F is intended for use in iterative scripts and has no effect absent -s<flt>\n");
		exit (1);
	}

/************/
/* read ifh */
/************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (srad > 0.0) {
		sprintf (tmproot, "%s_%.0fmm", imgroot, srad);
	} else {
		strcpy  (tmproot, imgroot);
	}
	sprintf (tmpfile, "%s.4dfp.img", tmproot);
	sprintf (ifhfile, "%s.4dfp.ifh", imgroot);

	fprintf (stdout, "Reading: %s\n", ifhfile);
	if (Getifh (ifhfile, &imgifh)) errr (program, ifhfile);
	isbig = strcmp (imgifh.imagedata_byte_order, "littleendian");
	if (!control) control = (isbig) ? 'b' : 'l';
	printf ("image dimensions \t%10d%10d%10d%10d\n",
		imgifh.matrix_size[0], imgifh.matrix_size[1], imgifh.matrix_size[2], imgifh.matrix_size[3]);
	printf ("image mmppix     \t%10.6f%10.6f%10.6f\n",
			imgifh.mmppix[0], imgifh.mmppix[1], imgifh.mmppix[2]);
	printf ("image center     \t%10.4f%10.4f%10.4f\n",
			imgifh.center[0], imgifh.center[1], imgifh.center[2]);
	printf ("image orientation\t%10d\n", imgifh.orientation);

/********************************************************************/
/* check that optionally specified mask is dimensionally consistent */
/********************************************************************/
	if (mask) {
		sprintf (mskfile, "%s.4dfp.img", mskroot);
		sprintf (ifhfile, "%s.4dfp.ifh", mskroot);
		fprintf (stdout, "Reading: %s\n", ifhfile);
		if (Getifh (ifhfile, &mskifh)) errr (program, ifhfile);
		isbigm = strcmp (mskifh.imagedata_byte_order, "littleendian");
		status = imgifh.orientation - mskifh.orientation;
		for (k = 0; k < 3; k++) {
			status |= imgifh.matrix_size[k] - mskifh.matrix_size[k];
			status |= (fabs ((double) (imgifh.mmppix[k] - mskifh.mmppix[k])) > 0.0001);
		}
		if (status) {
			fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot, mskroot);
			exit (-1);
		}
	}

	dim = 1; for (k = 0; k < 3; k++) dim *= imgifh.matrix_size[k];
	if (!(imgv = (float *) malloc (dim * sizeof (float)))) errm (program);
	if (!(imgr = (float *) malloc (dim * sizeof (float)))) errm (program);
	if (!(imgm = (float *) malloc (dim * sizeof (float)))) errm (program);
	nx = imgifh.matrix_size[0];
	ny = imgifh.matrix_size[1];
	nz = imgifh.matrix_size[2]; 

/*****************************************/
/* virtual flip instead of x4dfp2ecat () */
/*****************************************/
	vrtflip_ (&imgifh.orientation, imgifh.matrix_size, imgifh.center, imgifh.mmppix, centerr, mmppixr);
	printf ("atlas mmppix     \t%10.6f%10.6f%10.6f\n", mmppixr[0], mmppixr[1], mmppixr[2]); 
	printf ("atlas center     \t%10.4f%10.4f%10.4f\n", centerr[0], centerr[1], centerr[2]);

/****************************************************/
/* read filtered image or create using hpshere_4dfp */
/****************************************************/
	if (forceblur || access (tmpfile, R_OK)) {
		sprintf (command, "hsphere_4dfp %s %.4f", imgfile, srad);
		printf ("%s\n", command);
		if (system (command)) errw (program, tmpfile);
	}
	printf ("Reading: %s\n", tmpfile);
	if (!(imgfp = fopen (tmpfile, "rb")) || eread (imgv, dim, isbig, imgfp)
	|| fclose (imgfp)) errr (program, tmpfile);

/************************************************************************************************/
/* intial extremum memory allocation (realloc on unassigned pointer => unpredictable core dump) */
/************************************************************************************************/
	if (!(ppos = (EXTREMUM *) malloc ((mpos = MSIZE) * sizeof (EXTREMUM)))) errm (program);
	if (!(pneg = (EXTREMUM *) malloc ((mneg = MSIZE) * sizeof (EXTREMUM)))) errm (program);
/*******************/
/* compile extrema */
/*******************/
	printf ("peak value     thresholds %10.4f to %10.4f\n", vtneg, vtpos);
	printf ("peak curvature thresholds %10.4f to %10.4f\n", ctneg, ctpos);
	printf ("compiling extrema slice");
	for (iz = 1; iz < imgifh.matrix_size[2] - 1; iz++) {printf (" %d", iz + 1); fflush (stdout);
	for (iy = 1; iy < imgifh.matrix_size[1] - 1; iy++) {
	for (ix = 1; ix < imgifh.matrix_size[0] - 1; ix++) {
		i = ix + nx*(iy + ny*iz);
		dvdx[0] = 0.5*(-imgv[i - 1]	+ imgv[i + 1]);
		dvdx[1] = 0.5*(-imgv[i - nx]	+ imgv[i + nx]);
		dvdx[2] = 0.5*(-imgv[i - nx*ny]	+ imgv[i + nx*ny]);
		d2vdx2[0] = -2.0*imgv[i] + imgv[i - 1]		+ imgv[i + 1];
		d2vdx2[1] = -2.0*imgv[i] + imgv[i - nx]		+ imgv[i + nx];
		d2vdx2[2] = -2.0*imgv[i] + imgv[i - nx*ny]	+ imgv[i + nx*ny];
		for (k = 0; k < 3; k++) w[k] = -dvdx[k] / d2vdx2[k];
		for (del2v = k = 0; k < 3; k++) del2v += d2vdx2[k] / (mmppixr[k] * mmppixr[k]);
		fndex[0] = (float) (ix + 1) + w[0];
		fndex[1] = (float) (iy + 1) + w[1];
		fndex[2] = (float) (iz + 1) + w[2];
		for (k = 0; k < 3; k++) x[k] = fndex[k] * mmppixr[k] - centerr[k];

		if (imgv[i] > 0.0
		&& del2v <= -ctpos
		&& imgv[i] > imgv[i - 1]	&& imgv[i] > imgv[i + 1]
		&& imgv[i] > imgv[i - nx]	&& imgv[i] > imgv[i + nx]
		&& imgv[i] > imgv[i - nx*ny]	&& imgv[i] > imgv[i + nx*ny]) {
			imgvalx_ (imgv, &nx, &ny, &nz, centerr, mmppixr, x, &t, &l);
			if (l < 1) {printf ("undefined imgvalx point\n"); continue;}
			if (t < vtpos) continue;
			if (mpos <= npos) {
				mpos += MSIZE;
				if (!(ppos = (EXTREMUM *) realloc (ppos, mpos * sizeof (EXTREMUM)))) errm (program);
			}
			for (k = 0; k < 3; k++) ppos[npos].x[k] = x[k];
			ppos[npos].v = t;
			ppos[npos].del2v = del2v;
			ppos[npos].weight = 1.0;
			ppos[npos].nvox = ppos[npos].killed = 0;
			npos++;
		}

		if (imgv[i] < 0.0
		&& del2v >= -ctneg
		&& imgv[i] < imgv[i - 1]	&& imgv[i] < imgv[i + 1]
		&& imgv[i] < imgv[i - nx]	&& imgv[i] < imgv[i + nx]
		&& imgv[i] < imgv[i - nx*ny]	&& imgv[i] < imgv[i + nx*ny]) {
			imgvalx_ (imgv, &nx, &ny, &nz, centerr, mmppixr, x, &t, &l);
			if (l < 1) {printf ("undefined imgvalx point\n"); continue;}
			if (t > vtneg) continue;
			if (mneg <= nneg) {
				mneg += MSIZE;
				if (!(pneg = (EXTREMUM *) realloc (pneg, mneg * sizeof (EXTREMUM)))) errm (program);
			}
			for (k = 0; k < 3; k++) pneg[nneg].x[k] = x[k];
			pneg[nneg].v = t;
			pneg[nneg].del2v = del2v;
			pneg[nneg].weight = 1.0;
			pneg[nneg].nvox = pneg[nneg].killed = 0;
			nneg++;
		}
	}}}
	printf ("\n"); fflush (stdout);
	printf ("before sorting npos=%d nneg=%d\n", npos, nneg); npos0 = npos; nneg0 = nneg;

	qsort ((void *) ppos, npos, sizeof (EXTREMUM), pcompare);
	qsort ((void *) pneg, nneg, sizeof (EXTREMUM), pcompare);

	printf ("positive peaks\n");
	ntot = 0;
	peaklist (stdout, ppos, npos, 0); npos1 = ntot;
	if (dthresh > 0.0) consolidate (ppos, npos);

	printf ("negative peaks\n");
	ntot = 0;
	peaklist (stdout, pneg, nneg, 0); nneg1 = ntot;
	if (dthresh > 0.0) consolidate (pneg, nneg);

	printf ("final peak list\n");
	ntot = 0;
	peaklist (stdout, ppos, npos, 0);
	peaklist (stdout, pneg, nneg, 0);

/********************************************/
/* combine positive and negative peak lists */
/********************************************/
	if (!(pall = (EXTREMUM *) malloc (ntot * sizeof (EXTREMUM)))) errm (program);
	nall = 0;
	for (npos2 = i = 0; i < ntop && i < npos; i++) if (!ppos[i].killed) {pall[nall++] = ppos[i]; npos2++;}
	for (nneg2 = i = 0; i < ntop && i < nneg; i++) if (!pneg[i].killed) {pall[nall++] = pneg[i]; nneg2++;}
	assert (ntot == nall);
	free (ppos); free (pneg);

	if (orad == 0.0) goto DONE;
/***************************/
/* create output ROI image */
/***************************/
		printf ("Reading: %s\n", imgfile);
		if (!(imgfp = fopen (imgfile, "rb")) || eread (imgv, dim, isbig,  imgfp)
		|| fclose (imgfp)) errr (program, imgfile);
	if (mask) {
		printf ("Reading: %s\n", mskfile);
		if (!(imgfp = fopen (mskfile, "rb")) || eread (imgm, dim, isbigm, imgfp)
		|| fclose (imgfp)) errr (program, mskfile);
	}

/****************************/
/* count voxels in each ROI */
/****************************/
	for (iz = 0; iz < imgifh.matrix_size[2]; iz++) {
	for (iy = 0; iy < imgifh.matrix_size[1]; iy++) {
	for (ix = 0; ix < imgifh.matrix_size[0]; ix++) {
		i = ix + nx*(iy + ny*iz);
		fndex[0] = (float) (ix + 1);
		fndex[1] = (float) (iy + 1);
		fndex[2] = (float) (iz + 1);
		for (k = 0; k < 3; k++) ptmp[0].x[k] =  fndex[k]*mmppixr[k] - centerr[k];
		dmin = 1.e6;
		for (j = 0; j < nall; j++) {
			t = (float) pdist (ptmp, pall + j);
			if (t < dmin) {
				jmin = j;
				dmin = t;
			}
		}
		k = !mask || (fabs (imgm[i]) > 1.e-37);
		if (k && (dmin < orad)) pall[jmin].nvox++;
	}}}

/*******************************/
/* apply voxel count criterion */
/*******************************/
	for (i = 0; i < ntot; i++) if (pall[i].nvox < Nvoxcrit) pall[i].killed = 1;
	if (!(pallN = (EXTREMUM *) malloc (ntot * sizeof (EXTREMUM)))) errm (program);
	for (nallN = i = 0; i < ntot; i++) if (!pall[i].killed) pallN[nallN++] = pall[i];
	free (pall);

	printf ("creating ROI slice");
	for (iz = 0; iz < imgifh.matrix_size[2]; iz++) {printf (" %d", iz + 1); fflush (stdout);
	for (iy = 0; iy < imgifh.matrix_size[1]; iy++) {
	for (ix = 0; ix < imgifh.matrix_size[0]; ix++) {
		i = ix + nx*(iy + ny*iz);
		imgr[i] = 0.0;
		fndex[0] = (float) (ix + 1);
		fndex[1] = (float) (iy + 1);
		fndex[2] = (float) (iz + 1);
		for (k = 0; k < 3; k++) ptmp[0].x[k] =  fndex[k]*mmppixr[k] - centerr[k];
		dmin = 1.e6;
		for (j = 0; j < nallN; j++) {
			t = (float) pdist (ptmp, pallN + j);
			if (t < dmin) {
				jmin = j;
				dmin = t;
			}
		}
		k = !mask || (fabs (imgm[i]) > 1.e-37);
		if (k && (dmin < orad)) imgr[i] = jmin + 2;
	}}}
	printf ("\n"); fflush (stdout);

	if (!(ptr = strrchr (imgroot, '/'))) ptr = imgroot; else ptr++;
	if (strlen (trailer)) {
		sprintf (outroot, "%s_ROI_%s", ptr, trailer);
	} else {
		sprintf (outroot, "%s_ROI", ptr);
	}
	sprintf (outfile, "%s.4dfp.img", outroot);
	printf ("Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "wb")) || ewrite (imgr, dim, control, outfp)
	|| fclose (outfp)) errw (program, outfile);

/********************************************/
/* output ROI ifh and make analyze hdr file */
/********************************************/
	if (Writeifh (program, outfile, &imgifh, control)) errw (program, outroot);
	sprintf (outfile, "%s.4dfp.ifh", outroot);
	if (!(outfp = fopen (outfile, "a"))) errw (program, outfile);
	for (i = 0; i < nallN; i++) {
		fprintf (outfp, "region names\t:= %-5droi_%+03d_%+03d_%+03d%10d\n", i,
			nint (pallN[i].x[0]), nint (pallN[i].x[1]), nint (pallN[i].x[2]), pallN[i].nvox);
	}
	fclose (outfp);

/********************/
/* make analyze hdr */
/********************/
	sprintf (command, "ifh2hdr %s -r%d", outroot, nallN + 1);
	printf ("%s\n", command); status |= system (command);

/******************/
/* output ROI rec */
/******************/
	sprintf   (outfile, "%s.4dfp.img", outroot);
	startrece (outfile, argc, argv, rcsid, control);
	sprintf   (recfile, "%s.rec", outfile);
	if (!(outfp = fopen (recfile, "a"))) errw (program, recfile);
	fprintf (outfp, "Peak value     thresholds %10.4f to %10.4f\n", vtneg, vtpos);
	fprintf (outfp, "Peak curvature thresholds %10.4f to %10.4f\n", ctneg, ctpos);
	fprintf (outfp, "Peak counts before sorting pos=%d neg=%d\n", npos0, nneg0);
	fprintf (outfp, "Peak counts after  sorting pos=%d neg=%d\n", npos1, nneg1);
	fprintf (outfp, "Peak counts after %.4f mm threshold consolidation pos=%d neg=%d\n", dthresh, npos2, nneg2);
	ntot = 0; peaklist (outfp, pallN, nallN, 1);
	free (pallN);
	fprintf (outfp, "N.B.: indices count from 1 and include orientation specific 4dfptoanalyze flips\n");
	fclose (outfp);
	catrec (tmpfile);
	if (mask) catrec (mskfile);
	endrec ();
	sprintf (command, "brec %s -2", recfile); if (!quiet) system (command);

DONE:	free (imgv); free (imgr); free (imgm);
	exit (status);
}

/************************************************/
/* convert (flip) array indices 4dfp to analyze */
/************************************************/
void index_flip (int orientation, int* imgdim, float* fndex) {
	switch (orientation) {
		case 4:	fndex[0] = (float) imgdim[0] + 1. - fndex[0];
		case 3:	fndex[2] = (float) imgdim[2] + 1. - fndex[2];
		case 2:	fndex[1] = (float) imgdim[1] + 1. - fndex[1];
		return;
		break;
	}
}
