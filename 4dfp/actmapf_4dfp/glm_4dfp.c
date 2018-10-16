/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/glm_4dfp.c,v 1.26 2012/11/20 03:25:49 avi Exp avi $*/
/*$Log: glm_4dfp.c,v $
 * Revision 1.26  2012/11/20  03:25:49  avi
 * increase image filename buffers to 4*MAXL to accommodate long paths/filenames
 *
 * Revision 1.25  2011/12/29  21:17:56  avi
 * use intelligenet logic in original allocation of format buffer
 *
 * Revision 1.24  2011/04/12  22:19:56  avi
 * make use of conc_seek()
 *
 * Revision 1.23  2011/03/14  06:39:36  avi
 * option -b
 *
 * Revision 1.22  2011/02/07  04:47:40  avi
 * use imgb to speed inner loop computation of imgr and imgb
 *
 * Revision 1.21  2011/02/06  03:37:48  avi
 * squeeze skipped frames out of f on read
 * replace dmatinv_() with in-line C code
 *
 * Revision 1.20  2011/02/02  01:37:40  avi
 * imgr_flag
 *
 * Revision 1.19  2011/02/02  00:28:23  avi
 * execution speed and memory usage optimization
 *
 * Revision 1.18  2011/01/24  06:04:35  avi
 * invert design matrix only when necessary (finvt_flag)
 *
 * Revision 1.17  2011/01/19  04:23:35  avi
 * increase number of profile columns to 8192 and length of profile input line to 81920 chars
 *
 * Revision 1.16  2010/12/13  07:06:35  avi
 * trap "NaN" in input profiles
 *
 * Revision 1.15  2010/02/11  06:28:55  avi
 * profile input column limit -> 2048
 *
 * Revision 1.14  2009/12/01  04:32:23  avi
 * increase MAXC to 8192 and srgv to length 2*MAXL
 *
 * Revision 1.13  2008/11/26  02:07:48  avi
 * move cndnum test from dglm_4dfp.f to main c program
 *
 * Revision 1.12  2007/09/11  01:21:58  avi
 * format array length determined by profile length (eliminate MAXF)
 * clarify that lengths of format profile and input data must match
 * generalize profile reading to tolerate comment lines
 * eliminate f_init() and f_exit()
 *
 * Revision 1.11  2007/05/04  15:27:22  avi
 * enable partial correlation mapping
 *
 * Revision 1.10  2006/09/24  01:09:44  avi
 * Solaris 10
 *
 * Revision 1.9  2006/05/04  05:24:52  avi
 * option -Z inhibit removing constant from regressors
 *
 * Revision 1.8  2005/09/05  05:26:12  avi
 * increase profile input line buffer to 4096 chars
 * write frame number to stdout during processing
 *
 * Revision 1.7  2005/09/05  00:54:10  avi
 * double precision GLM inversion
 *
 * Revision 1.6  2005/05/24  06:39:46  avi
 * -C option (read previously computed coefficients from 4dfp image)
 *
 * Revision 1.5  2005/01/12  07:31:15  avi
 * MAXF 4096 -> 16384
 *
 * Revision 1.4  2004/11/27  05:45:35  avi
 * replace conc_io.c subroutines with conc.c subroutines
 * eliminate dependence on libmri (get4dfpdimN())
 *
 * Revision 1.3  2004/09/07  19:43:54  avi
 * optionally read conc files
 *
 * Revision 1.2  2004/06/10  00:15:22  avi
 * MAXF 1024 -> 4096
 *
 * Revision 1.1  2004/05/26  05:30:46  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <unistd.h>			/* getpid () */
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>
#include <conc.h>			/* defines MAXL */

#define MAXC		8192		/* maximum number of input profle columns */	
#define RTRAIL		"resid"
#define CTRAIL		"coeff"
#define BTRAIL		"tbeta"
#define PTRAIL		"pcorr"
#define TTRAIL		"tcorr"
#define MAX(a,b)	(a>b? a:b)
#define MIN(a,b)	(a<b? a:b)

int	expandf (char *string, int len);							/* expandf.c */
float	fimg_mode (float* fimg, int nval);							/* fimg_mode.c */
void	dmatinv_ (double *a, int *n, double *det);						/* dmatinv.f */
void    deigen_ (double *e, double *w, int *n);							/* deigen.f */

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: glm_4dfp.c,v 1.26 2012/11/20 03:25:49 avi Exp avi $";
static char	program[MAXL];
static int	debug = 0;

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

double dzeromean (double *f, int npts) {
	int		i;
	double		u;

	for (u = i = 0; i < npts; i++) u += f[i];
	u /= npts;
	for (i = 0; i < npts; i++) f[i] -= u;
	if (debug) for (i = 0; i < npts; i++) {
		printf ("%4d %10.6f\n", i + 1, f[i]);
	}
	return u;
}

double dunitvar (double *f, int npts) {
	int		i;
	double		v;

	for (v = i = 0; i < npts; i++) v += f[i]*f[i];
	v /= npts;
	for (i = 0; i < npts; i++) f[i] /= sqrt (v);
	if (debug) for (i = 0; i < npts; i++) {
		printf ("%4d %10.6f\n", i + 1, f[i]);
	}
	return v;
}

void dmatlst (double **a, int n) {
	int	i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) printf ("%11.7f", a[j][i]);
		printf ("\n");
	}
}

void usage (char *program) {
	fprintf (stderr, "Usage:\t%s <format> <profile> <4dfp|conc input>\n", program);
	fprintf (stderr, "e.g.,\t%s \"4x124+\" doubletask.txt b1_rmsp_dbnd_xr3d_norm\n", program);
	fprintf (stderr, "\toption\n");
	fprintf (stderr, "\t-Z\tsupress automatic removal of mean from input regressors\n");
	fprintf (stderr, "\t-C<str>\tread  partial beta coefficients from specified 4dfp image (default compute)\n");
	fprintf (stderr, "\t-o[str]\tsave  partial beta images with specified trailer (default = \"%s\")\n", CTRAIL);
	fprintf (stderr, "\t-R   compute  partial beta images as percent modulation\n");
	fprintf (stderr, "\t-b[str]\tsave  total   beta images with specified trailer (default = \"%s\")\n", BTRAIL);
	fprintf (stderr, "\t-p[str]\tsave  partial corr images with specified trailer (default = \"%s\")\n", PTRAIL);
	fprintf (stderr, "\t-t[str]\tsave  total   corr images with specified trailer (default = \"%s\")\n", TTRAIL);
	fprintf (stderr, "\t-r[str]\tsave  residual timeseries with specified trailer (default = \"%s\")\n", RTRAIL);
	fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default input endian)\n");
	fprintf (stderr, "N.B.:\tconc files must have extension \"conc\"\n");
	fprintf (stderr, "N.B.:\t<profile> lists temporal profiles (ASCII npts x ncol; '#' introduces comments)\n");
	fprintf (stderr, "N.B.:\t<profile> line limits are %d chars and %d fields\n", MAXC*10, MAXC);
	fprintf (stderr, "N.B.:\tabsent -C, options -o and -r require design matrix inversion; dimension limit 256\n");
	exit (1);
}

int main (int argc, char *argv[]) {
/**********************/
/* filename variables */
/**********************/
	FILE		*tmpfp, *imgfp, *coefp, *outfp, *profp;
	char		coeroot[4*MAXL], coefile[4*MAXL];	/* coefficients image */
	char		imgroot[4*MAXL], imgfile[4*MAXL];	/* input 4dfp stack filename */
	char		outroot[4*MAXL], outfile[4*MAXL], tmpfile[4*MAXL], profile[4*MAXL];
	char		ctrail[MAXL] = CTRAIL;
	char		rtrail[MAXL] = RTRAIL;
	char		ptrail[MAXL] = PTRAIL;
	char		btrail[MAXL] = BTRAIL;
	char		ttrail[MAXL] = TTRAIL;

/**********************/
/* image dimensioning */
/**********************/
	CONC_BLOCK	conc_block;		/* conc i/o control block */
	IFH		ifhimg, ifhcoe;
	float		**beta;			/* regression  coefficients */
	float		**imgr;			/* correlation coefficients */
	float		**imgb;			/* bold data with time as first index */
	float		*imgt;			/* one volume */
	float		*imga;			/* voxelwise mean */
	float		*imgv;			/* voxelwise var */
	int		*imgm;			/* mask of voxels defined at all time points */
	float		voxdim[3], mmppix[3], center[3];
	int		jndex, imgdim[4], vdim;	/* image dimensions */
	int		isbig, isbigr;
	char		control = '\0';

/*************************/
/* timeseries processing */
/*************************/
	char		*format;
	float		fmin, fmax, fmode;
	double		**f, **a, **ainv, **finvt, **e, **w, **omega, *sigma, det, cndnum, errmax;
	int		ppts, npts, nnez, ncol, ncolp1;

/***********/
/* utility */
/***********/
	char		command[4*MAXL], string[MAXC*10 + 4], *ptr, *srgv[MAXC];
	int		c, i, j, k, m;
	double		q;

/*********/
/* flags */
/*********/
	int		conc_flag = 0;
	int		zeromean_flag = 1;
	int		read_coeff = 0;
	int		save_pbeta = 0;
	int		save_tbeta = 0;
	int		save_resid = 0;
	int		save_pcorr = 0;
	int		save_tcorr = 0;
	int		scale_rel = 0;
	int		finvt_flag;		/* invert design matrix */
	int		imgr_flag;
	int		status = 0;

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
				case 'Z': zeromean_flag = 0;		break;
				case 'R': scale_rel++;			break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
				case 'C': read_coeff++; getroot (ptr, coeroot);	*ptr = '\0'; break;
				case 'o': save_pbeta++;
					if (strlen (ptr)) strcpy (ctrail, ptr);	*ptr = '\0'; break;
				case 'b': save_tbeta++;
					if (strlen (ptr)) strcpy (btrail, ptr);	*ptr = '\0'; break;
				case 'r': save_resid++;
					if (strlen (ptr)) strcpy (rtrail, ptr);	*ptr = '\0'; break;
				case 'p': save_pcorr++;
					if (strlen (ptr)) strcpy (ptrail, ptr);	*ptr = '\0'; break;
				case 't': save_tcorr++;
					if (strlen (ptr)) strcpy (ttrail, ptr);	*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	m = strlen (argv[i]) + 2; if (!(format = (char *) calloc (m, sizeof (char)))) errm (program);
				strcpy  (format,  argv[i]);		k++; break;
			case 1:	strcpy  (profile, argv[i]);		k++; break;
			case 2:	getroot (argv[i], imgroot);
				conc_flag = (strstr (argv[i], ".conc") == argv[i] + strlen (imgroot));
									k++; break;
		}	
	}
	if (k < 3) usage (program);

/****************/
/* read profile */
/****************/
	printf ("Reading: %s\n", profile);
	if (!(profp = fopen (profile, "r"))) errr (program, profile);
	ppts = 0; while (fgets (string, MAXC*10 + 1, profp)) {
		if (!strrchr (string, '\n')) {
			fprintf (stderr, "%s: %s exceeds line length limit (%d chars)\n", program, profile, MAXC*10);
			exit (-1);
		}
		if (!(m = split (string, srgv, MAXC))) continue;
		if (!ppts) {
			ncol = m;
			printf ("ncol=%d\n", ncol);
		} else {
			if (m != ncol) {
				fprintf (stderr, "%s: %s format error\n", program, profile);
				exit (-1);
			}
		}
		ppts++;
	}

/****************/
/* parse format */
/****************/
	if (!(format = (char *) realloc (format, (ppts + 1)*sizeof (char)))) errm (program);
	if (k = expandf (format, ppts + 1)) exit (-1);
	printf ("%s\n", format);
	npts = strlen (format);
	for (nnez = k = 0; k < npts; k++) if (format[k] != 'x') nnez++;
	printf ("%s: time series defined for %d frames, %d exluded\n", program, npts, npts - nnez);
	if (ppts != npts) {
		fprintf (stderr, "%s: %s/format length mismatch\n", program, profile);
		exit (-1);
	}

	if (save_pcorr && ncol > 64) {
		fprintf (stderr, "%s: partial correlation matrix dimension limit (64) exceeded\n", program);
		exit (-1);
	}
	finvt_flag = save_pbeta || save_resid;
	if (finvt_flag && ncol > 256) {
		fprintf (stderr, "%s: design matrix inversion dimension limit (256) exceeded\n", program);
		exit (-1);
	}

	f = calloc_double2 (ncol, nnez);
	rewind (profp);
	i = k = 0; while (fgets (string, MAXC*10 + 1, profp)) {
		if (!(m = split (string, srgv, MAXC))) continue;
		if (format[i] != 'x') {
			for (j = 0; j < ncol; j++) {
				if (isnan (q = atof (srgv[j]))) {
					fprintf (stderr, "%s: %s contains NaN\n", program, profile);
					exit (-1);
				}
				f[j][k] = q;
			}
			k++;
		}
		i++;
	}
	fclose (profp);
	assert (i == npts); assert (k == nnez);

/**********************************************/
/* make input profile zero mean unit variance */
/**********************************************/
	if (!(sigma = (double *) malloc (ncol * sizeof (double)))) errm (program);
	for (j = 0; j < ncol; j++) {
		if (zeromean_flag) dzeromean (f[j], nnez);
		sigma[j] = sqrt (dunitvar (f[j], nnez));
	}

	if (finvt_flag || save_pcorr) {
/**************************************/ 
/* compute profile correlation matrix */
/**************************************/ 
		a = calloc_double2 (ncol, ncol);
		for (i = 0; i < ncol; i++) {
			for (j = i; j < ncol; j++) {
				for (q = k = 0; k < nnez; k++) q += f[i][k]*f[j][k];
				a[i][j] = a[j][i] = q/nnez;
				if (debug) printf ("a[%d,%d]=%10.6f\n", i, j, a[i][j]);
			}
		}
	}

	if (finvt_flag) {
/************************/
/* invert design matrix */
/************************/
		finvt = calloc_double2 (ncol, nnez);
		ainv  = calloc_double2 (ncol, ncol);
		e     = calloc_double2 (ncol, ncol);
		w     = calloc_double2 (ncol, ncol);

		for (i = 0; i < ncol; i++) for (j = 0; j < ncol; j++) ainv[i][j] = e[i][j] = a[i][j];
		deigen_ (e[0], w[0], &ncol);
		cndnum = e[0][0]/e[ncol - 1][ncol - 1];
		printf ("condition_number=%12.6e\n", cndnum);
		if (cndnum < 0. || cndnum > 1.e4) {
			fprintf (stderr, "%s: ill-conditioned design matrix\n", program);
			exit (-2);
		}
		dmatinv_ (ainv[0], &ncol, &det);
		printf ("det=%12.6e\n", det);

		for (i = 0; i < nnez; i++) for (j = 0; j < ncol; j++) {
			for (finvt[j][i] = k = 0; k < ncol; k++) finvt[j][i] += f[k][i]*ainv[j][k];
		}
		errmax = 0.;
		for (i = 0; i < ncol; i++) for (j = 0; j < ncol; j++) {
			for (e[j][i] = k = 0; k < nnez; k++) e[j][i] += finvt[i][k]*f[j][k];
			e[j][i] /= nnez;
			q = fabs (e[j][i]); if (i == j) q -= 1.;
			if (q > errmax) errmax = q;
		}
		printf ("maximum_error=%12.6e\n", errmax);

		free_double2 (e); free_double2 (ainv); free_double2 (w);
	}

/*****************************/
/* get 4dfp stack dimensions */
/*****************************/
	if (conc_flag) {
		conc_init (&conc_block, program);
		conc_open (&conc_block, imgroot);
		strcpy (imgfile, conc_block.lstfile);
		for (k = 0; k < 4; k++) imgdim[k] = conc_block.imgdim[k];
		isbig = conc_block.isbig;
	} else {
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		if (Getifh (imgfile, &ifhimg)) errr (program, imgfile);
		for (k = 0; k < 4; k++) imgdim[k] = ifhimg.matrix_size[k];
		isbig = strcmp (ifhimg.imagedata_byte_order, "littleendian");
		if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
		printf ("Reading: %s\n", imgfile);
	}
	if (!control) control = (isbig) ? 'b' : 'l';
	vdim = imgdim[0] * imgdim[1] * imgdim[2];
	if (imgdim[3] < npts) {
		fprintf (stderr, "%s: more defined npts (%d) than frames (%d)\n", program, npts, imgdim[3]);
		exit (-1);
	}

	if (finvt_flag || read_coeff) {
		beta = calloc_float2 (ncol, vdim);
	}
	imgr_flag = save_tcorr || save_pcorr || save_tbeta;
	if (imgr_flag) {
		imgr = calloc_float2 (ncol, vdim);
	}
	if (!(imga = (float *) calloc (vdim, sizeof (float)))) errm (program);
	if (!(imgv = (float *) calloc (vdim, sizeof (float)))) errm (program);
	if (!(imgt = (float *) calloc (vdim, sizeof (float)))) errm (program);
	if (!(imgm = (int *)   calloc (vdim, sizeof (int))))   errm (program);
	for (jndex = 0; jndex < vdim; jndex++) imgm[jndex] = 1;

	if (read_coeff) {
/**************************************/
/* read regression coefficient images */
/**************************************/
		sprintf (coefile, "%s.4dfp.ifh", coeroot);
		printf ("Reading: %s\n", coefile);
		if (Getifh (coefile, &ifhcoe)) errr (program, coefile);
		for (k = 0; k < 3; k++) status |= (imgdim[k] != ifhcoe.matrix_size[k]);
		if (status) {
			fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot, coeroot);
			exit (-1);
		}
		if (ifhcoe.matrix_size[k] != ncol) {
			fprintf (stderr, "%s: %s %s column count mismatch\n", program, coeroot, profile);
			exit (-1);
		}
		isbigr = strcmp (ifhcoe.imagedata_byte_order, "littleendian");
		sprintf (coefile, "%s.4dfp.img", coeroot);
		printf ("Reading: %s\n", coefile);
		if (!(coefp = fopen (coefile, "rb")) || eread (beta[0], vdim*ncol, isbigr, coefp)
		|| fclose (coefp)) errr (program, coefile);
	} else {
/*****************************************/
/* compute regression coefficient images */
/*****************************************/
		imgb = calloc_float2 (vdim, nnez);
		printf ("reading frame");
		for (k = i = 0; i < npts; i++) {
			if (format[i] == 'x') continue;
			printf (" %d", i + 1); fflush (stdout);
			if (conc_flag) {
				conc_seek     (&conc_block, i);
				conc_read_vol (&conc_block, imgt);
			} else {
				if (fseek (imgfp, (long) i*vdim*sizeof(float), SEEK_SET)
				||  eread (imgt, vdim, isbig, imgfp)) errr (program, imgfile);
			}
			for (jndex = 0; jndex < vdim; jndex++) {
				if (isnan (imgt[jndex]) || imgt[jndex] == (float) 1.e-37) imgm[jndex] = 0;
				imga[jndex]	+= imgt[jndex];
				imgv[jndex]	+= imgt[jndex]*imgt[jndex];
				imgb[jndex][k]	=  imgt[jndex];
			}
			k++;
		} printf ("\n");  fflush (stdout);
		for (jndex = 0; jndex < vdim; jndex++) {
			imga[jndex] /= nnez;
			imgv[jndex] /= nnez;
			imgv[jndex] -= imga[jndex]*imga[jndex];
			if (imgv[jndex] < 1.e-6) imgm[jndex] = 0;
		}
		printf ("computing profile column");
		for (j = 0; j < ncol; j++) {
			printf (" %d", j + 1); fflush (stdout);
			for (jndex = 0; jndex < vdim; jndex++) {
				if (finvt_flag)	{
					for (k = 0; k < nnez; k++) beta[j][jndex] += imgb[jndex][k]*finvt[j][k];
					beta[j][jndex] /= nnez;
				}
				if (imgr_flag) {
					for (k = 0; k < nnez; k++) imgr[j][jndex] += imgb[jndex][k]*    f[j][k];
					imgr[j][jndex] /= nnez*sqrt(imgv[jndex]);
				}
			}
		} printf ("\n");  fflush (stdout);
		free_float2 (imgb); 
	}

/********************************/
/* create temp process log file */
/*******************************/
	sprintf (tmpfile, "temp%d", getpid ());
	if (!(tmpfp = fopen (tmpfile, "w"))) errw (program, tmpfile);
	fprintf (tmpfp, "Timepoint counts: counted=%d  skipped=%d\n", nnez, npts - nnez);
	fprintf (tmpfp, "%s\n", format);
	if (read_coeff) fprintf (tmpfp, "Regression coefficients read from %s\n", coeroot);
	fclose (tmpfp);

	if (save_resid) {
/*****************************/
/* remove profile components */
/*****************************/
		if (imgdim[3] != npts) {
			fprintf (stderr, "%s: profile length (%d)/frame count (%d) mismatch;",
				program, npts, imgdim[3]);
			fprintf (stderr, " residual computation not possible\n");
			exit (-1);
		}
		if (conc_flag) {
			conc_newe (&conc_block, rtrail, control);
			strcpy (outfile, conc_block.outfile);
			conc_rewind (&conc_block);
		} else {
			sprintf (outroot, "%s_%s", imgroot, rtrail);
			sprintf (outfile, "%s.4dfp.img", outroot);
			if (!(outfp =  fopen (outfile, "wb"))) errw (program, outfile);
			rewind (imgfp);
		}
		printf ("Writing: %s\n", outfile);
		printf ("computing residual frame"); fflush (stdout);
		for (k = i = 0; i < npts; i++) {
			printf (" %d", i + 1); fflush (stdout);
			if (conc_flag) {
				conc_read_vol (&conc_block, imgt);
			} else {
				if (eread (imgt, vdim, isbig, imgfp)) errr (program, imgfile);
			}
			if (format[i] != 'x') {
				for (jndex = 0; jndex < vdim; jndex++) {
					for (q = j = 0; j < ncol; j++) q += beta[j][jndex]*f[j][k];
					imgt[jndex] -= q;
				}
				k++;
			}
			if (conc_flag) {
				conc_write_vol (&conc_block, imgt);
			} else {
				if (ewrite (imgt, vdim, control, outfp)) errw (program, outfile);
			}
		} printf ("\n"); fflush (stdout);
		if (conc_flag) {
			status |= conc_ifh_hdr_rec (&conc_block, argc, argv, rcsid);
		} else {
			if (fclose (outfp)) errw (program, outfile);
			sprintf (command, "/bin/cp %s.4dfp.ifh %s.4dfp.ifh", imgroot, outroot);
			printf ("%s\n", command); status |= system (command);
			sprintf (command, "ifh2hdr -r4000 %s", outroot); system (command);
			printf ("%s\n", command); status |= system (command);
		}
		startrecle (outfile, argc, argv, rcsid, control);
		catrec (tmpfile);
		catrec (imgfile);
		endrec ();
	}

	if (save_pbeta) {
/****************************************/
/* save partial regression coefficients */
/****************************************/
		sprintf (outroot, "%s_%s", imgroot, ctrail);
		sprintf (outfile, "%s.4dfp.img", outroot);
		if (scale_rel) fmode = fimg_mode (imga, vdim);
		fmin = FLT_MAX; fmax = -fmin;
		for (j = 0; j < ncol; j++) {
			for (jndex = 0; jndex < vdim; jndex++) {
				if (scale_rel) {
					if (imga[jndex] > 0.25*fmode) {
						beta[j][jndex] *= 100./imga[jndex];
					} else {
						beta[j][jndex] = 1.e-37;
					}
				}
				fmin = MIN (fmin, beta[j][jndex]);
				fmax = MAX (fmax, beta[j][jndex]);
			}
		}
		printf ("Writing: %s\n", outfile);
		printf ("Max = %10.3f,\tMin = %10.3f\n", fmax, fmin);
		if (!(outfp = fopen (outfile, "wb")) || ewrite (beta[0], vdim*ncol, control, outfp)
		|| fclose (outfp)) errw (program, outfile);
		imgdim[3] = ncol;
		if (conc_flag) {
			writeifhmce (program, outfile, imgdim,
				conc_block.voxdim, conc_block.orient, conc_block.mmppix, conc_block.center, control);
		} else {
			writeifhmce (program, outfile, imgdim,
				ifhimg.scaling_factor, ifhimg.orientation, ifhimg.mmppix, ifhimg.center, control);
		}
		sprintf (command, "ifh2hdr -r%dto%d %s", (int) (fmin-0.5), (int) (fmax+0.5), outfile);
		printf ("%s\n", command); status |= system (command);
		startrecle (outfile, argc, argv, rcsid, control);
		if (scale_rel) {
			sprintf (command, "Regression coefficients expressed percent modulation\n"); printrec (command);
		}
		catrec (tmpfile);
		catrec (imgfile);
		endrec ();
	}

	if (save_tcorr) {
/***************************************/
/* save total correlation coefficients */
/***************************************/
		sprintf (outroot, "%s_%s", imgroot, ttrail);
		sprintf (outfile, "%s.4dfp.img", outroot);
		fmin = FLT_MAX; fmax = -fmin;
		for (jndex = 0; jndex < vdim; jndex++) for (j = 0; j < ncol; j++) {
			if (!imgm[jndex]) {
				 imgr[j][jndex] = 1.e-37;
			} else {
				fmin = MIN (fmin, imgr[j][jndex]);
				fmax = MAX (fmax, imgr[j][jndex]);
			}
		}
		printf ("Writing: %s\n", outfile);
		printf ("Max = %10.6f,\tMin = %10.6f\n", fmax, fmin);
		if (!(outfp = fopen (outfile, "wb")) || ewrite (imgr[0], vdim*ncol, control, outfp)
		|| fclose (outfp)) errw (program, outfile);
		imgdim[3] = ncol;
		if (conc_flag) {
			writeifhmce (program, outfile, imgdim,
				conc_block.voxdim, conc_block.orient, conc_block.mmppix, conc_block.center, control);
		} else {
			writeifhmce (program, outfile, imgdim,
				ifhimg.scaling_factor, ifhimg.orientation, ifhimg.mmppix, ifhimg.center, control);
		}
		sprintf (command, "ifh2hdr -r-1to1 %s", outfile);
		printf ("%s\n", command); status |= system (command);
		startrecle (outfile, argc, argv, rcsid, control);
		catrec (tmpfile);
		catrec (imgfile);
		endrec ();
	}

	if (save_tbeta) {
/**************************************************/
/* compute and save total regression coefficients */
/**************************************************/
		sprintf (outroot, "%s_%s", imgroot, btrail);
		sprintf (outfile, "%s.4dfp.img", outroot);
		if (!(outfp = fopen (outfile, "wb"))) errw (program, outfile);
		printf ("Writing: %s\n", outfile);
		fmin = FLT_MAX; fmax = -fmin;
		for (j = 0; j < ncol; j++) {
			for (jndex = 0; jndex < vdim; jndex++) {
				if (!imgm[jndex]) {
					 imgt[jndex] = 1.e-37;
				} else {
					imgt[jndex] = imgr[j][jndex]*sqrt(imgv[jndex])/sigma[j];
					fmin = MIN (fmin, imgt[jndex]);
					fmax = MAX (fmax, imgt[jndex]);
				}
			}
			if (ewrite (imgt, vdim, control, outfp)) errw (program, outfile);
		}
		if (fclose (outfp)) errw (program, outfile);
		printf ("Max = %10.6f,\tMin = %10.6f\n", fmax, fmin);
		imgdim[3] = ncol;
		if (conc_flag) {
			writeifhmce (program, outfile, imgdim,
				conc_block.voxdim, conc_block.orient, conc_block.mmppix, conc_block.center, control);
		} else {
			writeifhmce (program, outfile, imgdim,
				ifhimg.scaling_factor, ifhimg.orientation, ifhimg.mmppix, ifhimg.center, control);
		}
		sprintf (command, "ifh2hdr -r%.0fto%.0f %s", fmin, fmax, outfile);
		printf ("%s\n", command); status |= system (command);
		startrecle (outfile, argc, argv, rcsid, control);
		catrec (tmpfile);
		catrec (imgfile);
		endrec ();
	}

	if (save_pcorr) {
/*****************************************************/
/* compute and save partial correlation coefficients */
/*****************************************************/
		ncolp1 = ncol + 1;
		omega = calloc_double2 (ncolp1, ncolp1);
		sprintf (outroot, "%s_%s", imgroot, ptrail);
		sprintf (outfile, "%s.4dfp.img", outroot);
		fmin = FLT_MAX; fmax = -fmin;
		for (jndex = 0; jndex < vdim; jndex++) {
			if (!imgm[jndex]) {
				for (j = 0; j < ncol; j++) imgr[j][jndex] = 1.e-37;
			} else {
				omega[0][0] = 1.;
				for (j = 0; j < ncol; j++) {
					omega[j + 1][0] = omega[0][j + 1] = imgr[j][jndex];
					for (i = 0; i < ncol; i++) omega[j + 1][i + 1] = a[j][i];
				}
				dmatinv_ (omega[0], &ncolp1, &det);
				for (j = 0; j < ncol; j++) {
					imgr[j][jndex] = -omega[0][j + 1]/sqrt (omega[0][0]*omega[j + 1][j + 1]);
					fmin = MIN (fmin, imgr[j][jndex]);
					fmax = MAX (fmax, imgr[j][jndex]);
				}
			}
		}
		printf ("Writing: %s\n", outfile);
		printf ("Max = %10.6f,\tMin = %10.6f\n", fmax, fmin);
		if (!(outfp = fopen (outfile, "wb")) || ewrite (imgr[0], vdim*ncol, control, outfp)
		|| fclose (outfp)) errw (program, outfile);
		imgdim[3] = ncol;
		if (conc_flag) {
			writeifhmce (program, outfile, imgdim,
				conc_block.voxdim, conc_block.orient, conc_block.mmppix, conc_block.center, control);
		} else {
			writeifhmce (program, outfile, imgdim,
				ifhimg.scaling_factor, ifhimg.orientation, ifhimg.mmppix, ifhimg.center, control);
		}
		sprintf (command, "ifh2hdr -r-1to1 %s", outfile);
		printf ("%s\n", command); status |= system (command);
		startrecle (outfile, argc, argv, rcsid, control);
		catrec (tmpfile);
		catrec (imgfile);
		endrec ();
		free_double2 (omega);
	}

 	remove (tmpfile);
	if (conc_flag) {
		conc_free (&conc_block);
	} else {
		if (fclose (imgfp)) errr (program, imgfile);
	}
	free_double2 (f);
	if (finvt_flag || save_pcorr)	free_double2 (a);
	if (finvt_flag)			free_double2 (finvt);
	if (finvt_flag || read_coeff)	free_float2 (beta);
	if (imgr_flag)			free_float2 (imgr);
	free (imgt); free (imga); free (imgv); free (imgm);
	free (format); free (sigma);
	exit (status);
}
