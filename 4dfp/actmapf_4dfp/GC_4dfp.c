/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/GC_4dfp.c,v 1.7 2007/08/18 04:56:50 avi Exp $*/
/*$Log: GC_4dfp.c,v $
 * Revision 1.7  2007/08/18  04:56:50  avi
 * trap NaN -> 1.e-37 in write_4dfp()
 *
 * Revision 1.6  2007/07/30  03:50:01  avi
 * normalize lagged covariance (gamma) by true number of 2nd order samples/lag (ntot[l])
 *
 * Revision 1.5  2007/07/08  04:47:53  avi
 * Solaris 10 endian compliant
 *
 * Revision 1.4  2007/07/08  03:52:51  avi
 * option -F (write Fx,y Fx->y Fy->x Fx.y images)
 *
 * Revision 1.3  2006/08/06  03:34:21  avi
 * options -D -Z -a
 *
 * Revision 1.2  2006/08/05  02:41:11  avi
 * computes directed influences
 *
 * Revision 1.1  2006/08/03  06:38:36  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <unistd.h>			/* getpid () */
#include <assert.h>
#include <rec.h>
#include <Getifh.h>
#include <endianio.h>
#include <conc.h>

#define MAXO		7		/* maximum lag (model order) */ 
#define MAXL		256
#define MAXF		16384		/* maximum npts coded in format */
#define MAX(a,b)	(a>b? a:b)
#define MIN(a,b)	(a<b? a:b)

extern int	expandf (char *string, int len);					/* expandf.c */
extern void	compute_gc_  (int *dimx, int *dimy, float *Sx, float *Sy, float *Sq);	/* fGC.f */
extern void	compute_gc1_ (int *dimx, int *dimy, float *Sx, float *Sy, float *Sq,	/* fGC.f */
					float *cx, float *cy, float *cq, float *F);
extern void	solve_var_  (int *ncol, int *iorder, float *G, float *cc, float *ca,	/* fGC.f */
					float *cb, float *A, float *det, float *S);

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: GC_4dfp.c,v 1.7 2007/08/18 04:56:50 avi Exp $";
static char	program[MAXL];
static int	debug = 0;
static int	verbose = 0;
static int	conc_flag = 0;

float ***calloc_float3 (int n1, int n2, int n3) {
	unsigned int	i, j;
	float		***a;

	if (!(a = (float ***) malloc (n1 * sizeof (float **)))) errm (program);
	if (!(a[0] = (float **) malloc (n1 * n2 * sizeof (float *)))) errm (program);
	if (!(a[0][0] = (float *) calloc (n1 * n2 * n3, sizeof (float)))) errm (program);
	for (i = 0; i < n1; i++) {
		a[i] = a[0] + n2*i;
		for (j = 0; j < n2; j++) {
			a[i][j] = a[0][0] + n3*(n2*i + j);
		}
	}
	return a;
}

void free_float3 (float ***a) {
	free (a[0][0]);
	free (a[0]);
	free (a);
}

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

double rectifier (double x) {
	return (x > 0.0) ? x : 0.0;
}

int write_4dfp (char *imgroot, char *trailer, float *imgt, int *imgdim, char control, IFH *ifh, int argc, char *argv[]);

void usage (char *program) {
	fprintf (stderr, "Usage:\t%s <format> <4dfp|conc input> <order>\n", program);
	fprintf (stderr, "e.g.,\t%s \"4(4x190+)\" VB20579_rmsp_faln_dbnd_xr3d_atl.conc 2\n", program);
	fprintf (stderr, "\toption\n");
	fprintf (stderr, "\t-w<str>\tspecify timecourse profile file (one or more columns)\n");
	fprintf (stderr, "\t-i<int>\tuse only specified column (counting from 1) of timecourse profile\n");
	fprintf (stderr, "\t-a<str>\tappend specifed string to map output\n");
	fprintf (stderr, "\t-g\twrite lagged covariance (gamma) 4dfp stack\n");
	fprintf (stderr, "\t-D\twrite difference of directed influences (Fx->y - Fy->x) map\n");
	fprintf (stderr, "\t-Z\twrite Geweke \"N(0,2)\" measure (difference of square roots) map\n");
	fprintf (stderr, "\t-F\twrite Fx,y, Fx->y, Fy->x, Fx.y map stack\n");
	fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default input endian)\n");
	fprintf (stderr, "N.B.:\tconc files must have extension \"conc\"\n");
	fprintf (stderr, "N.B.:\teffective frame count is determined by <format>\n");
	fprintf (stderr, "N.B.:\t'x' frames in format are not counted\n");
	exit (1);
}

int main (int argc, char *argv[]) {
/**********************/
/* filename variables */
/**********************/
	CONC_BLOCK	conc_block;			/* conc i/o control block */
	FILE		*tmpfp, *profp, *imgfp;
	char		imgroot[MAXL], imgfile[MAXL];	/* input 4dfp stack filename */
	char		outfile[MAXL], profile[MAXL] = "";
	char		usrtrail[MAXL] = "";		/* appended to map outfile root */

/**********************/
/* image dimensioning */
/**********************/
	IFH		ifh;
	float		**imgt;				/* (iorder + 1) volumes */
	float		*imga;				/* volume mean */
	float		*imgd, *imgz, *imgb, **imgF;	/* output volumes */
	int		imgdim[4], vdim;		/* image dimensions */
	int		isbig;
	char		control = '\0';

/*************************/
/* timeseries processing */
/*************************/
	char		*str, format[MAXF];
	int		npts, nnez, ncol = 0, dimx = 0, icol = -1;
	int		ntot[MAXO + 1];
	float		**f;			/* f[icol][ipts] */
	float		***Gxx;			/* Gxx[l][jcol][icol] */
	float		***Gxy;			/* Gxy[l][jcol][ivox] */
	float		***Gyx;			/* Gyx[l][icol][ivox] */
	float		**Gyy;			/* Gyy[l][ivox] */
	float 		cscale = 1.0;

/*******************/
/* VAR (GC) arrays */
/*******************/
	float		***Gq, ***Aq;			/* full model gamma and AR coefficients */
	float		***Gx, ***Ax, ***Gy, ***Ay;	/* x and y    gamma and AR coefficients */
	float		**Sq, **Sx, **Sy;		/* sigma in Yule-Walker equations */
	float		*cc, *ca, *cb, det;		/* scratch arrays for AR solver */
	float		*xx, *yy, *qq, F[4];		/* scratch arrays for compute_gc1() */
	int		dimx1;				/* dimx + 1 */
	double		Gconst;				/* Geweke's (klp - 1)/3 */

/***********/
/* utility */
/***********/
	char		command[MAXL], string[MAXL], *ptr, *srgv[MAXL];
	int		c, i, j, k, l, m;
	int		one = 1;

/*********/
/* flags */
/*********/
	int		iorder = 0;
	int		status = 0;
	int		debug = 0;
	int		verbose = 0;
	int		write_gamma = 0;
	int		write_D = 0;
	int		write_Z = 0;
	int		write_B = 0;
	int		write_F = 0;

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
				case 'v': verbose++;			break;
				case 'g': write_gamma++;		break;
				case 'F': write_F++;			break;
				case 'D': write_D++;			break;
				case 'B': write_B++;			break;
				case 'Z': write_Z++;			break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
				case 'a': strcpy (usrtrail, ptr);	*ptr = '\0'; break;
				case 'i': icol = atof (ptr) - 1;	*ptr = '\0'; break;
				case 'w': strcpy (profile, ptr);	*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	strcpy  (format, argv[i]);	k++; break;
			case 1:	getroot (argv[i], imgroot);
				conc_flag = (strstr (argv[i], ".conc") == argv[i] + strlen (imgroot));
								k++; break;
			case 2:	iorder = atoi (argv[i]);	k++; break;
		}	
	}
	if (k < 3) usage (program);
	if (iorder > MAXO) {
		fprintf (stderr, "%s: model order cannot exceed %d\n", program, MAXO);
		usage (program);
	}

/****************/
/* parse format */
/****************/
	if (k = expandf (format, MAXF)) exit (k);
	if (debug) printf ("%s\n", format);
	npts = strlen (format);
	for (l = 0; l <= iorder; l++) ntot[l] = 0;
	for (nnez = k = 0; k < npts; k++) {
		if (format[k] != 'x') nnez++;
		for (l = 0; l <= iorder; l++) if (format[k] != 'x' && format[k + l] != 'x') ntot[l]++;
	}
	fprintf (stderr, "%s: time series defined for %d frames, %d exluded\n", program, npts, npts - nnez);
	for (l = 0; l <= iorder; l++) printf ("ntot[%d] = %d\n", l, ntot[l]);

/*****************************/
/* get 4dfp stack dimensions */
/*****************************/
	if (conc_flag) {
		conc_init (&conc_block, program);
		conc_open (&conc_block, imgroot);
		strcpy (imgfile, conc_block.lstfile);
		for (k = 0; k < 4; k++) imgdim[k] = conc_block.imgdim[k];
		if (Getifh (conc_block.imgfile0[0], &ifh)) errr (program, conc_block.imgfile0[0]);
		isbig = conc_block.isbig;
	} else {
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		if (Getifh (imgfile, &ifh)) errr (program, imgfile);
		for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];
		isbig = strcmp (ifh.imagedata_byte_order, "littleendian");
		if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
		printf ("Reading: %s\n", imgfile);
	}
	if (!control) control = (isbig) ? 'b' : 'l';
	vdim = imgdim[0] * imgdim[1] * imgdim[2];
	if (imgdim[3] < npts) {
		fprintf (stderr, "%s: format defined for more points than %s\n", program, imgroot);
		exit (-1);
	}
	if (imgdim[3] > npts) {
		fprintf (stderr, "%s warning: format defined for less points than %s\n", program, imgroot);
	}

/*****************/
/* parse profile */
/*****************/
	if (strlen (profile)) {
		fprintf (stderr, "Reading: %s\n", profile);
		if (!(profp = fopen (profile, "r"))) errr (program, profile);
		k = 0; while (fgets (string, MAXL, profp)) {
			if (!(m = split (string, srgv, MAXL))) continue;
			if (!k) {
				ncol = m;
			} else {
				if (m != ncol) {
					fprintf (stderr, "%s: %s format error\n", program, profile);
					exit (-1);
				}
			}
			k++;
		}
		if (k < npts) {
			fprintf (stderr, "%s: format defined for more points than %s\n",
					program, profile);
			exit (-1);
		}
		if (k > npts) {
			fprintf (stderr, "%s warning: format defined for less points than %s\n",
					program, profile);
		}
	}

	if (icol >= 0) {
		if (icol >= ncol) {
			fprintf (stderr, "%s: requested column of %s does not exist\n", program, profile);
			usage (program);
		}
		dimx = 1;
	} else {
		dimx = ncol;
	}
	fprintf (stderr, "npts=%d ncol=%d dimx=%d\n", npts, ncol, dimx);

/*******************/
/* allocate memory */
/*******************/
	f	= calloc_float2 (dimx, npts);
	if (!(imga = (float *) calloc (vdim, sizeof (float)))) errm (program);
	if (!(imgd = (float *) calloc (vdim, sizeof (float)))) errm (program);
	if (!(imgz = (float *) calloc (vdim, sizeof (float)))) errm (program);
	if (!(imgb = (float *) calloc (vdim, sizeof (float)))) errm (program);
	imgF	= calloc_float2 (         4, vdim);
	imgt	= calloc_float2 (iorder + 1, vdim);
	Gxx	= calloc_float3 (iorder + 1, dimx, dimx);
	Gxy	= calloc_float3 (iorder + 1, dimx, vdim);
	Gyx	= calloc_float3 (iorder + 1, dimx, vdim);
	Gyy	= calloc_float2 (iorder + 1,       vdim);

/****************/
/* read profile */
/****************/
	if (strlen (profile)) {
		rewind (profp);
		i = 0; while (i < npts) {
			fgets (string, MAXL, profp);
			if (!(m = split (string, srgv, MAXL))) continue;
			if (icol >= 0) {
				f[0][i] = atof (srgv[icol]);
			} else {
				for (j = 0; j < ncol; j++) f[j][i] = atof (srgv[j]);
			}
			i++;
		}
		fclose (profp);
		assert (i == npts);
	}

/*****************************/
/* compute mean image volume */
/*****************************/
	printf ("computing mean volume time point");
	for (k = 0; k < npts; k++) {printf(" %d", k + 1); fflush (stdout);
		if (conc_flag) {
			conc_read_vol (&conc_block, imgt[0]);
		} else {
			if (eread (imgt[0], vdim, isbig, imgfp)) errr (program, imgfile);
		}
		if (format[k] != 'x') for (i = 0; i < vdim; i++) imga[i] += imgt[0][i];
	}
	printf("\n");
	for (i = 0; i < vdim; i++) imga[i] /= nnez;

/*****************************************/
/* reposition read pointer to data start */
/*****************************************/
	if (conc_flag) {
		conc_rewind (&conc_block);
	} else {
		rewind (imgfp);
	}

/***************************************/
/* load first iorder volumes into imgt */
/***************************************/
	for (m = 0; m < iorder; m++) {	/* (pointer to lag 0 volume of imgt) % (iorder + 1) */
		if (conc_flag) {
			conc_read_vol (&conc_block, imgt[m]);
		} else {
			if (eread (imgt[m], vdim, isbig, imgfp)) errr (program, imgfile);
		}
		for (i = 0; i < vdim; i++) imgt[m][i] -= imga[i];
	}

/*****************************/
/* compute cross covariances */
/*****************************/
	printf ("computing covariances time point");
	for (k = 0; k < npts - l; k++) {printf(" %d", k + 1); fflush (stdout);
		if (conc_flag) {
			conc_read_vol (&conc_block, imgt[m]);
		} else {
			if (eread (imgt[m], vdim, isbig, imgfp)) errr (program, imgfile);
		}
		for (i = 0; i < vdim; i++) imgt[m][i] -= imga[i];
		m++; m %= iorder + 1;
		assert (m == k % (iorder + 1));
		for (l = 0; l <= iorder; l++) if (format[k] != 'x' && format[k + l] != 'x') {
			for (i = 0; i < vdim; i++) {
				Gyy[l][i] += imgt[(m + l) % (iorder + 1)][i] * imgt[m % (iorder + 1)][i];
				for (j = 0; j < dimx; j++) {
					Gyx[l][j][i] += imgt[(m + l) % (iorder + 1)][i] * f[j][k];
					Gxy[l][j][i] +=  f[j][k + l] * imgt[m % (iorder + 1)][i];
				}
			}
			for (j = 0; j < dimx; j++) {
			for (i = 0; i < dimx; i++) {
				Gxx[l][j][i] += f[i][k + l] * f[j][k];
			}}
		}
	}
	printf("\n");
	for (l = 0; l <= iorder; l++ ) {
		for (i = 0; i < vdim; i++) {
			Gyy[l][i] /= ntot[l];
			for (j = 0; j < dimx; j++) {
				Gyx[l][j][i] /= ntot[l];
				Gxy[l][j][i] /= ntot[l];
			}
		}
		for (j = 0; j < dimx; j++) {
		for (i = 0; i < dimx; i++) {
			Gxx[l][j][i] /= ntot[l];
		}}
	}

/*************/
/* GC memory */
/*************/
	k = dimx + 1;					/* dimensionality of full model */
	Gq = calloc_float3 (iorder + 1,    k,    k);	/* lagged cross covariance (gamma) */
	Aq = calloc_float3 (iorder    ,    k,    k);	/* AR coefficients */
	Sq = calloc_float2 (               k,    k);	/* sigma */
	Gx = calloc_float3 (iorder + 1, dimx, dimx);
	Ax = calloc_float3 (iorder    , dimx, dimx);
	Sx = calloc_float2 (            dimx, dimx);
	Gy = calloc_float3 (iorder + 1,  one,  one);
	Ay = calloc_float3 (iorder    ,  one,  one);
	Sy = calloc_float2 (             one,  one);
	dimx1 = dimx + 1;
	k = dimx1*dimx1*iorder;		/* solve_var_() scratch space of maximum dimension */
	if (!(cc = (float *) malloc (k * k * sizeof (float)))) errm (program);
	if (!(ca = (float *) malloc (k     * sizeof (float)))) errm (program);
	if (!(cb = (float *) malloc (k     * sizeof (float)))) errm (program);
	if (!(xx = (float *) malloc ( dimx *  dimx * sizeof (float)))) errm (program);	/* compute_gc1_() scratch space */
	if (!(yy = (float *) malloc (    1 *     1 * sizeof (float)))) errm (program);
	if (!(qq = (float *) malloc (dimx1 * dimx1 * sizeof (float)))) errm (program);

	Gconst = (dimx*1*iorder - 1)/3.;	/* Geweke's (klp - 1)/3 */

	if (debug) printf ("profile alone\n");
	solve_var_ (&dimx, &iorder, &Gxx[0][0][0], cc, ca, cb, &Ax[0][0][0], &det, &Sx[0][0]);	/* Gx == Gxx */

	for (l = 0; l <= iorder; l++) {
		for (j = 0; j < dimx; j++) {
		for (i = 0; i < dimx; i++) {
			Gq[l][j][i] = Gxx[l][j][i];
		}}
	}
/******************/
/* loop on voxels */
/******************/
	if (dimx) for (k = 0; k < vdim; k++) {
		if (debug) printf ("voxel offset = %d\n", k);
		for (l = 0; l <= iorder; l++) {
			Gq[l][dimx][dimx] = Gy[l][0][0] = Gyy[l][k];
			for (j = 0; j < dimx; j++) {
				Gq[l][dimx][j] = Gxy[l][j][k];
				Gq[l][j][dimx] = Gyx[l][j][k];
			}
		}
		if (debug) printf ("voxel alone\n");
		solve_var_ (&one,   &iorder, &Gy[0][0][0], cc, ca, cb, &Ay[0][0][0], &det, &Sy[0][0]);	/* voxel alone */
		if (debug) {
			printf ("first profile leads voxel by one frame %10.4f\n", Gxy[1][0][k]);
			printf ("voxel leads first profile by one frame %10.4f\n", Gyx[1][0][k]);
		}
		if (debug) printf ("full model\n");
		solve_var_ (&dimx1, &iorder, &Gq[0][0][0], cc, ca, cb, &Aq[0][0][0], &det, &Sq[0][0]);	/* full model */
		if (debug)	compute_gc_  (&dimx, &one, &Sx[0][0], &Sy[0][0], &Sq[0][0]);
				compute_gc1_ (&dimx, &one, &Sx[0][0], &Sy[0][0], &Sq[0][0], xx, yy, qq, F);
		imgd[k] = F[1] - F[2];
		imgb[k] = F[1]/(F[1] + F[2]);
/***************************************/
/* /sqrt(2.) converts N(0,2) to N(0,1) */
/***************************************/
/*		imgz[k] = sqrt (rectifier ((double) nnez*F[1] - Gconst))
			- sqrt (rectifier ((double) nnez*F[2] - Gconst));	*/
/*		imgz[k] = sqrt (fabs (nnez*F[1] - Gconst)) - sqrt (fabs (nnez*F[2] - Gconst));	*/
		imgz[k] = sqrt (nnez*F[1]) - sqrt (nnez*F[2]);
		for (i = 0; i < 4; i++) imgF[i][k] = F[i];
	}

/*********************/
/* write 4dfp output */
/*********************/
	if (write_gamma) {
		imgdim[3] = iorder + 1;
		write_4dfp (imgroot, "gamma", Gyy[0],  imgdim, control, &ifh, argc, argv);
	}
	if (write_F) {
		imgdim[3] = 4;
		sprintf (command, "GC%d_F", iorder);
		if (strlen (usrtrail)) sprintf (command, "%s_%s", command, usrtrail);
		write_4dfp (imgroot, command, imgF[0], imgdim, control, &ifh, argc, argv);
	}
	imgdim[3] = 1;
	if (write_D) {
		sprintf (command, "GC%d_D", iorder);
		if (strlen (usrtrail)) sprintf (command, "%s_%s", command, usrtrail);
		write_4dfp (imgroot, command, imgd,    imgdim, control, &ifh, argc, argv);
	}
	if (write_Z) {
		sprintf (command, "GC%d_Z", iorder);
		if (strlen (usrtrail)) sprintf (command, "%s_%s", command, usrtrail);
		write_4dfp (imgroot, command, imgz,    imgdim, control, &ifh, argc, argv);
	}
	if (0 && write_B) {
		sprintf (command, "GC%d_B", iorder);
		if (strlen (usrtrail)) sprintf (command, "%s_%s", command, usrtrail);
		write_4dfp (imgroot, command, imgb,    imgdim, control, &ifh, argc, argv);
	}


	if (conc_flag) {
		conc_free (&conc_block);
	} else {
		if (fclose (imgfp)) errr (program, imgfile);
	}

	free (imgd); free (imgz); free (imgb); free (imga); free_float2 (imgF);
	free_float2 (imgt);
	free_float2 (f);
	free_float3 (Gxx);
	free_float3 (Gxy);
	free_float3 (Gyx);
	free_float2 (Gyy);
	free (ca); free (cb); free (cc); free (xx); free (yy); free (qq);
	free_float3 (Gq); free_float3 (Aq);
	free_float3 (Gx); free_float3 (Ax);
	free_float3 (Gy); free_float3 (Ay);
	free_float2 (Sq); free_float2 (Sx); free_float2 (Sy);
	exit (status);
}

int write_4dfp (char *imgroot, char *trailer, float *imgt, int *imgdim, char control, IFH *ifh, int argc, char *argv[]) {
	FILE		*imgfp;
	char		*ptr, imgfile[MAXL], outfile[MAXL], command[MAXL];
	int		i, tdim, status;
	float		fmin, fmax;

	tdim = imgdim[0]*imgdim[1]*imgdim[2]*imgdim[3];
	fmin = FLT_MAX; fmax = -fmin;
	for (i = 0; i < tdim; i++) {
		if (isnan (imgt[i])) imgt[i] = (float) 1.e-37;
		fmin = MIN (fmin, imgt[i]);
		fmax = MAX (fmax, imgt[i]);
	}

/********************/
/* assemble outroot */
/********************/
	if (ptr = strrchr (imgroot, '/')) ptr++; else ptr = imgroot;
	sprintf (outfile, "%s_%s.4dfp.img", ptr, trailer);

/*********/
/* write */
/*********/
	printf ("Writing: %s\n", outfile);
	printf ("Max = %10.3f,\tMin = %10.3f\n", fmax, fmin);
	if (!(imgfp = fopen (outfile, "wb"))
	|| ewrite (imgt, tdim, control, imgfp)
	|| fclose (imgfp)) errw (program, outfile);

/*******/
/* ifh */
/*******/
	writeifhmce (program, outfile, imgdim,
		ifh->scaling_factor, ifh->orientation, ifh->mmppix, ifh->center, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr -r%dto%d %s", (int) (fmin-0.5), (int) (fmax+0.5), outfile);
	printf ("%s\n", command); status |= system (command);

/*******/
/* rec */
/*******/
	sprintf (imgfile, (conc_flag ? "%s.conc" : "%s.4dfp.img"), imgroot);
	startrecle (outfile, argc, argv, rcsid, control);
	catrec (imgfile);
	endrec ();

	return status;
}
