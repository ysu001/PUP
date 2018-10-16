/*$Header: /data/petsun4/data1/src_solaris/mprtir/RCS/partitiond_gfc_4dfp.c,v 1.9 2011/11/09 08:03:31 avi Exp $*/
/*$Log: partitiond_gfc_4dfp.c,v $
 * Revision 1.9  2011/11/09  08:03:31  avi
 * option -M (maximum intensity scaling limit)
 * correct bug in final output gfc image computation loop
 *
 * Revision 1.8  2011/11/09  06:00:50  avi
 * correct use of preblur
 *
 * Revision 1.7  2011/04/06  01:32:07  avi
 * correct alpha parameter in call to fimgblur3d_() (1.0 -> 0.1)
 *
 * Revision 1.6  2010/11/29  06:36:05  avi
 * correct appearance of usage and datfile
 *
 * Revision 1.5  2008/05/08  04:55:56  avi
 * several flow-control adjustments to improve stability and trap algorithmic failure
 *
 * Revision 1.4  2008/01/02  01:24:40  avi
 * correct memory access error in getxyzn0()
 *
 * Revision 1.3  2007/12/31  22:51:29  avi
 * correct rec file generation
 * Solaris10 and unix compliant
 *
 * Revision 1.2  2007/12/29  01:41:46  avi
 * cosmetics
 *
 * Revision 1.1  2007/12/28  21:58:55  avi
 * Initial revision
**/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <rec.h>
#include <ifh.h>
#include <Getifh.h>
#include <endianio.h>

#define ALPHA		0.1		/* relative height of imgblur3d_() smoothing kernel skirt */
#define UCHAR		unsigned char
#define MAXR		4096		/* regions */
#define MR		24		/* gfc regions */
#define MAXL		256		/* string length */
#define MAXD		65536		/* maximum image dimension */
#define MEMB		4096		/* region buffer memory increment */
#define LIMITER		8
#define BANDWIDTH	200.
#define DIFFCRIT	0.0002
#define ZTHRESH		180.
#define GFC_THRESH	0.
#define GFC_CEIL	10000.
#define SIGMA		1.0
#define SCONST		4.0		/* space constant in units of mm */
#define NORDER		2		/* gain field correction polynomial order */
#define NTERM		(((NORDER+3)*(NORDER+2)*(NORDER+1))/6)	/* number of polynomial coefficients */

typedef struct {
	unsigned short	ix;
	unsigned short	iy;
	unsigned short	iz;
	unsigned short	nn;
} XYZN;

typedef struct {
	float	u1, u2;
	int	count;
	int	ne;
	int	nnmax;
	XYZN	xyzn0;
	int	nxyzn;
	XYZN	*xyzn;
} REGION;

/*******************************/
/* partitiond global variables */
/*******************************/
static float		g_bandwidth = BANDWIDTH;
static float		g_sigma = SIGMA;
static float		g_zthresh = ZTHRESH, g_gfc_thresh = GFC_THRESH, g_gfc_ceil = GFC_CEIL;
static float		g_sconst = SCONST;
static REGION		g_region[MAXR];
static int		g_nr = 0, g_mr = MR;
static short		g_list[MAXR];
static int		g_nl;

static void enter (int ir) {
	int		i;

	if (--ir < 0) return;
	for (i = 0; i < g_nl; i++) if (ir == g_list[i]) return;
	g_list[g_nl++] = ir;
}

/***********************/
/* global image memory */
/***********************/
static int	g_dimension, g_xdim, g_ydim, g_zdim;
static float	g_voxdim[3];
static float	*g_img0,		/* original image */
		*g_img1,		/* corrected image */
		*g_imgg;		/* gain field image */
static short	*g_imgo;
static UCHAR	*g_imgb;		/* neighbors in same partition count */

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

static void	assign (int ir, XYZN xyzn);						/* below */
static void	getmean (void);								/* below */
static void	getxyzn0 (void);							/* below */
static void	listreg (FILE *fp);							/* below */
static int	rcompare (const void *p1, const void *p2);				/* below */
void	imgblur3d_ (float* fwhm, float* alpha, float* voxsiz, float* img0, int* imgdim, float* img1);	/* fimgblur.f */
void	fitgain3d_ (float *, int *nx, int *ny, int *nz, short *,
			float fv[], int *ncat, int *norder, float a[]);			/* fitgain3d.f */
void	evalgain3d_ (float *, int *nx, int *ny, int *nz, int *norder, float a[]);	/* fitgain3d.f */
void	fitgain3dd_ (float *, int *nx, int *ny, int *nz, short *,
			float fv[], int *ncat, int *norder, float a[], float *drms);	/* fitgain3dd.f */
void	powstring_ (int *norder, int *lenfield, char *string);				/* powsub.f */
void	fnegdef_ (float a[]);								/* fnegdef.f */

static char g_rcsid[] = "$Id: partitiond_gfc_4dfp.c,v 1.9 2011/11/09 08:03:31 avi Exp $";
int main (int argc, char *argv[]) {
	FILE 		*imgfp, *datfp;
	char		imgroot[MAXL], imgfile[MAXL], gfcfile[MAXL], ifhfile[MAXL];
	char		outroot[MAXL], outfile[MAXL], datfile[MAXL], volfile[MAXL], recfile[MAXL];
	char		program[MAXL], string[MAXL],  command[MAXL], *ptr;
	int		isbig;
	char		control = '\0';

	IFH		ifh;
	UCHAR		ix, iy, iz;
	XYZN		xyzn;
	int		iv, iv0, jv, imgdim[3];
	int		ir, jr, nn;
	int		c, i, j, k, m, n;
	double		dela, delg, u1;
	float		vmin, vmax, fwhm = 0., alpha = ALPHA, *img_orig=NULL;
	float		maxfac = 0.0;

/**************/
/* gain field */
/**************/
	int		norder = NORDER, nterm = NTERM;
	int		niter = 0, limiter = LIMITER;
	int		eight = 8;		/* field length for call to powstring_ () */
	float		a[NTERM],		/* polynomial coefficients */
			a0[NTERM];
	float		drms, diffcrit = DIFFCRIT;
	float		u1gfc[MAXR];
	double		sum0, sum, q, t;

/*********/
/* flags */
/*********/
	int		test = 0;
	int		status = 0;
	int		debug = 0;
	int		gfreeze = 0;
	int		verbose = 0;
	int		negdef = 0;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) if (*argv[i] == '-') {
		strcpy (string, argv[i]); ptr = string;
		while (c = *ptr++) switch (c) {
			case 'v': verbose++;		break;
			case 'n': negdef++;		break;
			case 'g': gfreeze++;		break;
			case '@': control = *ptr++;		*ptr = '\0'; break;
			case 'b': g_bandwidth	= atof (ptr);	*ptr = '\0'; break;
			case 'p': fwhm		= atof (ptr);	*ptr = '\0'; break;
			case 'e': diffcrit	= atof (ptr);	*ptr = '\0'; break;
			case 'i': g_sigma	= atof (ptr);	*ptr = '\0'; break;
			case 'l': limiter	= atoi (ptr);	*ptr = '\0'; break;
			case 'm': g_mr		= atoi (ptr);	*ptr = '\0'; break;
			case 's': g_sconst	= atof (ptr);	*ptr = '\0'; break;
			case 'z': g_zthresh	= atof (ptr);	*ptr = '\0'; break;
			case 'M': maxfac	= atof (ptr);	*ptr = '\0'; break;
			case 'r': getrange (ptr, &g_gfc_thresh, &g_gfc_ceil); *ptr = '\0'; break;
		}
	} else switch (k) {
		case 0: getroot (argv[i], imgroot); k++;  break;
	}
	if (k < 1) {
		fprintf (stderr, "Usage:\t%s <imgroot>\n", program);
		fprintf (stderr, " e.g.,\t%s vc1440_mpr_n4_111_t88.4dfp\n", program); /*??*/
		fprintf (stderr, "\toption\n");
		fprintf (stderr, "\t-g	freeze initial gain field\n");
		fprintf (stderr, "\t-n	force negative definite quadratic gain field\n");
		fprintf (stderr, "\t-v	verbose mode\n");
		fprintf (stderr, "\t-p<flt> pre-blur by specified FWHM in mm\n");
		fprintf (stderr, "\t-b<flt>\tspecify bandwidth in intensity units (default=%.1f)\n", BANDWIDTH);
		fprintf (stderr, "\t-e<flt>\tspecify drms convergence criterion (default=%f)\n", DIFFCRIT);
		fprintf (stderr, "\t-i<flt>\tspecify sigma (default=%f)\n", SIGMA);
		fprintf (stderr, "\t-l<int>\tspecify iteration limit (default=%d)\n", LIMITER);
		fprintf (stderr, "\t-m<flt>\tspecify gfc computation region count (default=%d)\n", MR);
		fprintf (stderr, "\t-s<flt>\tspecify space constant in mm (default=%f)\n", SCONST);
		fprintf (stderr, "\t-z<flt>\tspecify background threshold (default=%.1f)\n", ZTHRESH);
		fprintf (stderr, "\t-M<flt>\tspecify maximum correction factor\n");
		fprintf (stderr, "\t-r<flt>[to<flt>]\tspecify gfc range (default=%.1fto%.1f)\n", GFC_THRESH, GFC_CEIL);
		exit (1);
	}

/*******************************************/
/* read/initialize gain field coefficients */
/*******************************************/
	sprintf (gfcfile, "%s.4dfp.gfc", imgroot);
	if (imgfp = fopen (gfcfile, "r")) {
		printf ("Reading: %s\n", gfcfile);
		for (k = 0; k < nterm; k++) fscanf (imgfp, "%f", a + k);
		fclose (imgfp);
	} else {
		a[0] = 1.0;
		for (k = 1; k < nterm; k++) a[k] = 0.0;
	}
	for (k = 0; k < nterm; k++) fprintf (stdout, "%8.4f", a[k]); fprintf (stdout, "\n");

/******************************/
/* get input image dimensions */
/******************************/
	if (Getifh (imgroot, &ifh) != 0) errr (program, imgroot);
	isbig = strcmp (ifh.imagedata_byte_order, "littleendian");
	if (!control) control = (isbig) ? 'b' : 'l';

	for (k = 0; k < 3; k++) g_voxdim[k] = ifh.scaling_factor[k];
	if (debug) printf ("%10.6f%10.6f%10.6f\n", g_voxdim[0], g_voxdim[1], g_voxdim[2]);
	g_dimension = ((g_xdim = ifh.matrix_size[0])
		      *(g_ydim = ifh.matrix_size[1])
		      *(g_zdim = ifh.matrix_size[2]));
/************************************************/
/* check image index limits against XYZN typedef */
/*************************************************/
	if (g_xdim > MAXD || g_ydim > MAXD || g_zdim > MAXD) {
		fprintf (stderr, "%s: %s dimension exceeds %d\n", program, imgroot, MAXD);
		exit (-1);
	}
	if (ifh.number_of_bytes_per_pixel != 4) {
		printf ("%s: cannot process %d bytes per pixel", program, ifh.number_of_bytes_per_pixel);
		exit (-1);
	}

/************/
/* allocate */
/************/
	if (!(g_img1 = (float *) malloc (g_dimension * sizeof (float)))
	||  !(g_img0 = (float *) malloc (g_dimension * sizeof (float)))
	||  !(g_imgg = (float *) malloc (g_dimension * sizeof (float)))
	||  !(g_imgo = (short *) calloc (g_dimension,  sizeof (short)))
	||  !(g_imgb = (UCHAR *) calloc (g_dimension,  sizeof (UCHAR)))) errm (program);

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	printf ("Reading: %s\n", imgfile);
	if (!(imgfp = fopen (imgfile, "rb")) || eread (g_img0, g_dimension, isbig, imgfp)
	|| fclose (imgfp)) errr (program, imgfile);

/******************/
/* pre-blur image */
/******************/
	if (!gfreeze && fwhm > 0) {
		if (!(img_orig = (float *) malloc (g_dimension * sizeof (float)))) errm (program);
/*************************************/
/* img_orig will hold original image */
/*************************************/
		memcpy (img_orig, g_img0, g_dimension*sizeof(float));
		memcpy (g_img1,   g_img0, g_dimension*sizeof(float));
		printf ("blur FWHM = %.2f mm\n", fwhm);
		imgdim[0] = g_xdim; imgdim[1] = g_ydim; imgdim[2] = g_zdim;
		imgblur3d_ (&fwhm, &alpha, g_voxdim, g_img1, imgdim, g_img0);	/* g_img1 returned unusable */
/************************************************************************/
/* g_img0 now points to blurred image to be used in further computation */
/************************************************************************/
	}

/*******************/
/* compute outroot */
/*******************/
	if (ptr = strrchr (imgroot, '/')) ptr++; else ptr = imgroot;
	getroot (ptr, outroot);

/*********************/
/* apply initial gfc */
/*********************/
	evalgain3d_ (g_imgg, &g_xdim, &g_ydim, &g_zdim, &norder, a);
	for (sum = sum0 = i = 0; i < g_dimension; i++) {
		sum0 += g_img0[i];
		sum  += (g_img1[i] = g_img0[i]/g_imgg[i]);
	}
	for (i = 0; i < g_dimension; i++) g_img1[i] *= (sum0/sum);

	if (gfreeze) goto WRITE;

/****************/
/* open datfile */
/****************/
	sprintf (datfile, "%s_gfc.dat", outroot);
	if (!(datfp = fopen (datfile, "a"))) errw (program, datfile);
	powstring_ (&norder, &eight, string);
	fprintf (datfp, "#\t%s\n", string);
	fprintf (datfp, "%d\t", niter);
	for (k = 0; k < nterm; k++) fprintf (datfp, "%8.4f", a[k]);
	fprintf (datfp, "%8.4f\n", 0.0); fflush (datfp);

/***********************/
/* clear region buffer */
/***********************/
	memset (&g_region, '\0', MAXR * sizeof (REGION));
	iv0 = ir = 0;

/***************************/
/* find contiguous regions */
/***************************/
CONTIG:	while (ir < MAXR) {
		for (iv = iv0; iv < g_dimension; iv++) if (!g_imgo[iv]) break;
		if (iv == g_dimension) break;
		iv0 = iv;

		g_imgo[iv] = ir + 1;
		g_region[ir].u1 = (floor (g_img1[iv] / g_bandwidth) + 0.5) * g_bandwidth;
		k = iv;
		j = g_dimension;
		j /= g_zdim; iz = k / j; k -= iz * j;
		j /= g_ydim; iy = k / j; k -= iy * j;
		j /= g_xdim; ix = k / j; k -= ix * j;
		xyzn.ix = ix; xyzn.iy = iy; xyzn.iz = iz; xyzn.nn = 0;
		assign (ir, xyzn);

		do {
			m = g_region[ir].ne++;
			iv = g_region[ir].xyzn[m].ix + g_xdim*(g_region[ir].xyzn[m].iy + g_ydim*g_region[ir].xyzn[m].iz);
			for (j = 0; j < 6; j++) {
				xyzn = g_region[ir].xyzn[m];
				switch (j) {
				case 0:	k =  -1;	if (xyzn.ix < 1)		continue; xyzn.ix--; break;
				case 1:	k *= -1;	if (xyzn.ix > g_xdim - 2)	continue; xyzn.ix++; break;
				case 2:	k *= -g_xdim;	if (xyzn.iy < 1)		continue; xyzn.iy--; break;
				case 3:	k *= -1;	if (xyzn.iy > g_ydim - 2)	continue; xyzn.iy++; break;
				case 4:	k *= -g_ydim;	if (xyzn.iz < 1)		continue; xyzn.iz--; break;
				case 5:	k *= -1;	if (xyzn.iz > g_zdim - 2)	continue; xyzn.iz++; break;
				}
				jv = iv + k;
				if (g_imgo[jv] || g_img0[jv] < g_zthresh) continue;
				dela = g_img1[jv] - g_region[ir].u1;
				if (fabs (dela) < 0.5 * g_bandwidth) {
					g_imgo[jv] = ir + 1;
					assign (ir, xyzn);
				}
			}
			if (verbose && !(g_region[ir].count % 1000)) {
				xyzn = g_region[ir].xyzn[m];
				printf ("checkr0: region=%4d unexa=%8d coords=%4d%4d%4d u1=%8.2f\n",
					ir, g_region[ir].count - m, xyzn.ix, xyzn.iy, xyzn.iz, g_region[ir].u1);
			}
		} while (g_region[ir].count > g_region[ir].ne);

		if (g_region[ir].count == 1) {
			iv = g_region[ir].xyzn[0].ix + g_xdim*(g_region[ir].xyzn[0].iy + g_ydim*g_region[ir].xyzn[0].iz);
			g_imgo[iv] = -1;
			g_region[ir].count = g_region[ir].ne = 0;
		} else {
			ir++;
		}
	}

/************************/
/* sort regions by size */
/************************/
	qsort ((void *) g_region, ir, sizeof (REGION), rcompare);
	if (ir == MAXR) {
		if (verbose) printf ("ir %d -> %d iv0=%10d ivmax=%d\n", ir, MAXR/2, iv0, g_dimension);
		while (--ir > MAXR/2) {
			for (m = 0; m < g_region[ir].count; m++) {
				xyzn = g_region[ir].xyzn[m]; iv = xyzn.ix + g_xdim*(xyzn.iy + g_ydim*xyzn.iz);
				g_imgo[iv] = 0;
			}
			g_region[ir].count = g_region[ir].ne = 0;
		}
		goto CONTIG;
	}
	g_nr = ir;
	for (ir = g_nr; ir < MAXR; ir++) if (g_region[ir].xyzn) free (g_region[ir].xyzn);

	test = 1;
ITER:	getmean ();	/* analyze g_img1 */
	getxyzn0 ();	/* analyze g_img1 */
/********************/
/* reassign regions */
/********************/
	for (i = 0; i < g_dimension; i++) g_imgo[i] = 0;
	printf ("total number of regions %d\n", g_nr);
	for (ir = 0; ir < g_nr; ir++) {
		g_region[ir].ne = g_region[ir].count = 0;
		if (g_region[ir].u1 < g_gfc_thresh || g_region[ir].u1 > g_gfc_ceil) continue;
		assign (ir, g_region[ir].xyzn0);
		do {
			m = g_region[ir].ne++;
			xyzn = g_region[ir].xyzn[m];
			iv = xyzn.ix + g_xdim*(xyzn.iy + g_ydim*xyzn.iz);
			for (j = 0; j < 6; j++) {
				xyzn = g_region[ir].xyzn[m];
				switch (j) {
				case 0:	k =  -1;	if (xyzn.ix < 1)		continue; xyzn.ix--; break;
				case 1:	k *= -1;	if (xyzn.ix > g_xdim - 2)	continue; xyzn.ix++; break;
				case 2:	k *= -g_xdim;	if (xyzn.iy < 1)		continue; xyzn.iy--; break;
				case 3:	k *= -1;	if (xyzn.iy > g_ydim - 2)	continue; xyzn.iy++; break;
				case 4:	k *= -g_ydim;	if (xyzn.iz < 1)		continue; xyzn.iz--; break;
				case 5:	k *= -1;	if (xyzn.iz > g_zdim - 2)	continue; xyzn.iz++; break;
				}
				jv = iv + k;
				if (g_imgo[jv] || g_img0[jv] < g_zthresh) continue;
				dela = g_img1[jv] - g_region[ir].u1;
				delg = g_sconst * (g_img1[jv] - g_img1[iv]) / g_voxdim[j/2];
				if (sqrt (dela*dela + delg*delg) < g_sigma*g_bandwidth) {
					g_imgo[jv] = ir + 1;
					assign (ir, xyzn);
				}
			}
			if (verbose && !(g_region[ir].count % 1000)) {
				xyzn = g_region[ir].xyzn[m];
				printf ("checkr1: region=%4d unexa=%8d coords=%4d%4d%4d u1=%8.2f\n",
				ir, g_region[ir].count - m, xyzn.ix, xyzn.iy, xyzn.iz, g_region[ir].u1);
			}
		} while (g_region[ir].count > g_region[ir].ne);
	}
	getmean (); getxyzn0 ();	/* analyze g_img1 */
	listreg (stdout);
	if (g_region[0].u1 < g_bandwidth) {
		fprintf (stderr, "%s: algorithmic failure\n", program);
		exit (-1);
	}	

	if (++niter > limiter || !test) goto COUNT;
	fprintf (stdout, "%s: iteration %d\n", program, niter);
/***************/
/* compute gfc */
/***************/
	for (ir = 0; ir < g_nr; ir++) u1gfc[ir] = g_region[ir].u1;
	for (k = 0; k < nterm; k++) a0[k] = a[k];
	fitgain3dd_ (g_img0, &g_xdim, &g_ydim, &g_zdim, g_imgo, u1gfc, &g_mr, &norder, a, &drms);
	if (negdef) fnegdef_ (a);
	printf ("niter=%d drms=%.6f diffcrit=%.6f\n", niter, drms, diffcrit);
	fprintf (datfp, "%d\t", niter);
	for (k = 0; k < nterm; k++) fprintf (datfp, "%8.4f", a[k]);
	fprintf (datfp, "%8.4f\n", 100*drms); fflush (datfp);
	evalgain3d_ (g_imgg, &g_xdim, &g_ydim, &g_zdim, &norder, a);
	for (sum = i = 0; i < g_dimension; i++) sum += (g_img1[i] = g_img0[i]/g_imgg[i]);
	for (      i = 0; i < g_dimension; i++) g_img1[i] *= (sum0/sum);
/********************/
/* convergence test */
/********************/
	test = (drms > diffcrit);
	goto ITER;

COUNT:	fclose (datfp);
/*******************/
/* count neighbors */
/*******************/
	for (ir = 0; ir < g_nr; ir++) {
		for (g_nl = m = 0; m < g_region[ir].count; m++) {
			ix = g_region[ir].xyzn[m].ix;
			iy = g_region[ir].xyzn[m].iy;
			iz = g_region[ir].xyzn[m].iz;
			iv = ix + g_xdim*(iy + g_ydim*iz);
			nn = 0;
			k =  -1;	if (ix > 0)	     {if ((j = g_imgo[iv + k]) == ir + 1) nn++; else enter (j);}
			k *= -1;	if (ix < g_xdim - 1) {if ((j = g_imgo[iv + k]) == ir + 1) nn++; else enter (j);}
			k *= -g_xdim;	if (iy > 0)	     {if ((j = g_imgo[iv + k]) == ir + 1) nn++; else enter (j);}
			k *= -1;	if (iy < g_ydim - 1) {if ((j = g_imgo[iv + k]) == ir + 1) nn++; else enter (j);}
			k *= -g_ydim;	if (iz > 0)	     {if ((j = g_imgo[iv + k]) == ir + 1) nn++; else enter (j);}
			k *= -1;	if (iz < g_zdim - 1) {if ((j = g_imgo[iv + k]) == ir + 1) nn++; else enter (j);}
			if (nn > g_region[ir].nnmax) g_region[ir].nnmax = nn;
			g_imgb[iv] = g_region[ir].xyzn[m].nn = nn;
		}
		if (verbose) {
			printf ("%5d%10d%10d%5d%2d |", ir, g_region[ir].count, k, j, g_region[ir].nnmax);
			for (i = 0; i < g_nl; i++) printf (" %d", g_list[i]); printf ("\n");
		}
	}

/*********************************/
/* write gain field coefficients */
/*********************************/
	if (!(imgfp = fopen (gfcfile, "w"))) errw (program, gfcfile);
	printf ("Writing: %s\n", gfcfile);
	for (k = 0; k < nterm; k++) fprintf (imgfp, "%8.4f", a[k]); fprintf (imgfp, "\n");
	fclose (imgfp);

/*************************/
/* write corrected image */
/*************************/
WRITE:	for (sum0 = sum = i = 0; i < g_dimension; i++) {
		if (!gfreeze && fwhm > 0) g_img0[i] = img_orig[i];
		sum0 += g_img0[i];
		q = 1./g_imgg[i];
		t = (maxfac > 0.0) ? 2.*maxfac*tanh(.5*q/maxfac) : q;
		if (0) printf ("%10.6f %10.6f\n", q, t);
		sum  += (g_img1[i] = g_img0[i]*t);
	}
	for (i = 0; i < g_dimension; i++) g_img1[i] *= (sum0/sum);

	vmin = FLT_MAX; vmax = -vmin;
	for (i = 0; i < g_dimension; i++) {
		if (g_img1[i] < vmin) vmin = g_img1[i];
		if (g_img1[i] > vmax) vmax = g_img1[i];
	}
	sprintf (outfile, "%s_gfc.4dfp.img", outroot);
	fprintf (stdout, "Writing: %s\n", outfile);
	if (!(imgfp = fopen (outfile, "wb")) || ewrite (g_img1, g_dimension, control, imgfp)
	|| fclose (imgfp)) errw (program, outfile);
/***************/
/* ifh and hdr */
/***************/
	sprintf (ifhfile, "%s_gfc.4dfp.ifh", outroot);
	if (Writeifh (program, ifhfile, &ifh, control)) errw (program, ifhfile);
	sprintf (command, "ifh2hdr %s -r%.0f", ifhfile, vmax);
	system  (command);
/*******/
/* rec */
/*******/
	startrece (outfile, argc, argv, g_rcsid, control);
	sprintf (string, "%s gain field coefficients\n", imgroot);	printrec (string);
	powstring_ (&norder, &eight, string);				printrec (string);
	for (k = 0; k < nterm; k++) {sprintf (string, "%8.4f", a[k]);	printrec (string);} 	printrec ("\n");
	if (!gfreeze) {
		sprintf (recfile, "%s_gfc.4dfp.img.rec", outroot);
		if (!(imgfp = fopen (recfile, "a"))) errw (program, recfile); listreg (imgfp); fclose (imgfp);
	}
	catrec (imgfile);
	endrec ();
	if (gfreeze) goto FREE;

/*************************/
/* compute regions image */
/*************************/
	for (i = 0; i < g_dimension; i++) g_img1[i] = 0;
	for (ir = 0; ir < g_mr; ir++) {
		for (m = 0; m < g_region[ir].count; m++) {
			xyzn = g_region[ir].xyzn[m];
			iv = xyzn.ix + g_xdim*(xyzn.iy + g_ydim*xyzn.iz);
			g_img1[iv] = g_region[ir].u1;
		}
	}

/******************************/
/* write processing QA images */
/******************************/
	sprintf (volfile, "%s_gfc_vols.4dfp.img", outroot);
	printf ("Writing: %s\n", volfile);
	if (!(imgfp = fopen (volfile, "w"))) errw (program, volfile);
/*****************/
/* regions image */
/*****************/
	if (ewrite (g_img1, g_dimension, control, imgfp)) errw (program, volfile);
/****************/
/* border image */
/****************/
	for (i = 0; i < g_dimension; i++) g_img1[i] = g_imgb[i];
	if (ewrite (g_img1, g_dimension, control, imgfp)) errw (program, volfile);
/**************/
/* gain field */
/**************/
	for (i = 0; i < g_dimension; i++) g_img1[i] = g_imgg[i];
	if (ewrite (g_img1, g_dimension, control, imgfp)) errw (program, volfile);
	fclose (imgfp);
/***************/
/* ifh and hdr */
/***************/
	sprintf (ifhfile, "%s_gfc_vols.4dfp.ifh", outroot);
	ifh.matrix_size[3] = 3;
	if (Writeifh (program, ifhfile, &ifh, control)) errw (program, ifhfile);
	sprintf (string, "ifh2hdr %s -r%.0f", ifhfile, vmax);
	system  (string);
/*******/
/* rec */
/*******/
	startrece (volfile, argc, argv, g_rcsid, control);
	sprintf (string, "%s_gfc QA (3 volumes)\n", imgroot); printrec (string);
	catrec (outfile);
	printrec ("Volume 1: regions\nVolume 2: border\nVolume 3: gain field\n");
	endrec ();

/************/
/* clean up */
/************/
FREE:	free (g_img1); free (g_img0); free (g_imgg);
	free (g_imgo); free (g_imgb);
	if (img_orig) free (img_orig);
	for (ir = 0; ir < g_nr; ir++) free (g_region[ir].xyzn);
	exit (status);
}

static int rcompare (const void *p1, const void *p2) {
	REGION	*r1, *r2;

	r1 = (REGION *) p1;
	r2 = (REGION *) p2;
	return r2->count - r1->count;
}

static int rucompare (const void *p1, const void *p2) {
	REGION	*r1, *r2;

	r1 = (REGION *) p1;
	r2 = (REGION *) p2;
	        if (r1->u1 < g_gfc_thresh || r1->u1 > g_gfc_ceil) return +1;
	else	if (r2->u1 < g_gfc_thresh || r2->u1 > g_gfc_ceil) return -1;
	else	return r2->count - r1->count;
}

static void getmean () {
/********************************************/
/* compute region means and reassign g_imgo */
/********************************************/
	XYZN		xyzn;
	double		q, u1;
	int		ir, i, m, iv;

	for (ir = 0; ir < g_nr; ir++) {
		for (u1 = m = 0; m < g_region[ir].count; m++) {
			xyzn = g_region[ir].xyzn[m]; iv = xyzn.ix + g_xdim*(xyzn.iy + g_ydim*xyzn.iz);
			u1 += g_img1[iv];
		}
		g_region[ir].u1 = u1 / g_region[ir].count;
	}
	qsort ((void *) g_region, g_nr, sizeof (REGION), rcompare);
	while (!g_region[g_nr - 1].count) g_nr--;

	for (i = 0; i < g_dimension; i++) g_imgo[i] = 0;
	for (ir = 0; ir < g_nr; ir++) {
		for (m = 0; m < g_region[ir].count; m++) {
			xyzn = g_region[ir].xyzn[m]; iv = xyzn.ix + g_xdim*(xyzn.iy + g_ydim*xyzn.iz);
			g_imgo[iv] = ir + 1;
		}
	}
}

static void getxyzn0 () {
/*************************************/
/* compute region variance and xyzn0 */
/*************************************/
	double		q, u2, del2, del2min;
	XYZN		xyzn;
	int		ir, k, m, iv;
	int		onface;

	for (ir = 0; ir < g_nr; ir++) {
		if (!g_region[ir].count) {
			printf ("region %d is void\n", ir);
			continue;
		}
		g_region[ir].xyzn0 = g_region[ir].xyzn[0];	/* in case no other voxel is most typical */
		del2min = FLT_MAX;
		for (u2 = m = 0; m < g_region[ir].count; m++) {
			xyzn = g_region[ir].xyzn[m]; iv = xyzn.ix + g_xdim*(xyzn.iy + g_ydim*xyzn.iz);
			q = g_img1[iv] - g_region[ir].u1;
			del2 =  q*q;
			u2 += del2;
			onface = xyzn.ix == 0 || xyzn.ix == g_xdim - 1
			      || xyzn.iy == 0 || xyzn.iy == g_ydim - 1
			      || xyzn.iz == 0 || xyzn.iz == g_zdim - 1;
			if (onface) continue;
			k = 1;
			q = 0.5 * (g_img1[iv + k] - g_img1[iv - k]) * g_sconst / g_voxdim[0];	del2 += q*q; k *= g_xdim;
			q = 0.5 * (g_img1[iv + k] - g_img1[iv - k]) * g_sconst / g_voxdim[1];	del2 += q*q; k *= g_ydim;
			q = 0.5 * (g_img1[iv + k] - g_img1[iv - k]) * g_sconst / g_voxdim[2];	del2 += q*q;
			if (del2 < del2min) {
				del2min = del2;
				g_region[ir].xyzn0 = xyzn;
			}
		}
		g_region[ir].u2 = u2/g_region[ir].count;
	}
}

static void assign (int ir, XYZN xyzn) {
	int		i, j, k;
	int		debug = 0;

	if (debug && !(g_region[ir].count % 1000)) {
		printf ("assign: region=%4d voxel=%8d coords=%4d%4d%4d u=%8.2f\n",
		ir, g_region[ir].count, xyzn.ix, xyzn.iy, xyzn.iz, g_region[ir].u1);
	}
	if (!(g_region[ir].count < g_region[ir].nxyzn)) {
		g_region[ir].nxyzn += MEMB;
		if (g_region[ir].xyzn) {
			g_region[ir].xyzn = (XYZN *) realloc (g_region[ir].xyzn, g_region[ir].nxyzn * sizeof (XYZN));
		} else {
			g_region[ir].xyzn = (XYZN *)  malloc (g_region[ir].nxyzn * sizeof (XYZN));
		}
		if (!g_region[ir].xyzn) errm ("assign");
	}
	g_region[ir].xyzn[g_region[ir].count++] = xyzn;
	return;
}

static void listreg (FILE *fp) {
/******************************/
/* list top region statistics */
/******************************/
	int		ir;

	fprintf (fp, "region   voxels      mean      s.d.   ix   iy   iz\n");
	for (ir = 0; ir < g_mr; ir++) {
		fprintf (fp, "%-5d%10d%10.2f%10.4f%5d%5d%5d\n", ir + 1, g_region[ir].count, g_region[ir].u1,
		sqrt(g_region[ir].u2), g_region[ir].xyzn0.ix + 1, g_region[ir].xyzn0.iy + 1, g_region[ir].xyzn0.iz + 1);
	}
}
