/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/covariance.c,v 1.31 2012/05/29 04:26:17 avi Exp $*/
/*$Log: covariance.c,v $
 * Revision 1.31  2012/05/29  04:26:17  avi
 * option -t (remove trend)
 *
 * Revision 1.30  2011/02/21  00:35:13  avi
 * option -g (auxilliary profile regression)
 *
 * Revision 1.29  2011/01/19  05:56:16  avi
 * increase MAXC and MAXS to match glm_4dfp
 *
 * Revision 1.28  2010/07/30  23:40:19  avi
 * change 'readlable' command line option form '#' to 'L'
 *
 * Revision 1.27  2010/07/29  02:10:15  avi
 * option -#
 * end
 *
 * Revision 1.26  2009/04/05  23:43:12  avi
 * increase size of ROI label (MAXR) to 64 chars
 *
 * Revision 1.25  2009/01/12  03:49:07  avi
 * if mlag==0 output CC[RV] matrix
 * tune format to 3 or 6 decimals according to unitvar state
 *
 * Revision 1.24  2008/12/02  06:16:45  avi
 * eigen_() -> deigen_()
 * option -D (reduce dimensionality of profile)
 *
 * Revision 1.23  2008/05/22  04:38:04  avi
 * account for phase wrap when computing delay (tau)
 *
 * Revision 1.22  2008/05/22  01:43:45  avi
 * modify to support processing of very long input timeseries (e.g., as in ECoG)
 * check to keep mlag < npts
 *
 * Revision 1.21  2008/02/20  00:11:25  avi
 * list condition number with option -e
 *
 * Revision 1.20  2008/02/19  23:48:08  avi
 * linux conpliant (getroot -> getrootl, to avoid conflict with endianio.h)
 *
 * Revision 1.19  2008/02/16  01:07:51  avi
 * correct bug in *format alloc
 *
 * Revision 1.18  2008/02/15  22:32:45  avi
 * allocate format according to data
 * linux compliant
 *
 * Revision 1.17  2007/01/05  05:20:45  avi
 * format[] can be optionally expanded beyond default MAXF using -F
 *
 * Revision 1.16  2006/11/29  01:34:41  avi
 * correct malfunctioning unit_var()
 * format "%10.4f" -> " %9.3f"
 *
 * Revision 1.15  2006/11/26  00:05:28  avi
 * Solaris10
 *
 * Revision 1.14  2006/09/14  21:52:29  avi
 * trap case of more ROI lables than data columns
 *
 * Revision 1.13  2006/07/27  02:32:40  avi
 * change cross spectral analysis from FFT based to Bartlett-Tukey
 *
 * Revision 1.12  2006/07/09  04:48:51  avi
 * correct code although -r option not very useful because of inherently
 * noisy behavior of unsmoothed spectra
 *
 * Revision 1.11  2006/07/07  23:15:41  avi
 * output cross spectra
 * revise logic of extracting ROI names
 *
 * Revision 1.10  2006/02/19  05:36:13  avi
 * inhibit computation of cross-spectrum
 *
 * Revision 1.9  2005/12/16  02:55:33  avi
 * option -T (Tukey window smoothing of Bartlett power spectrum)
 *
 * Revision 1.8  2005/12/14  05:57:17  avi
 * output Bartlett window smoothed power spectrum
 *
 * Revision 1.7  2005/12/13  05:04:15  justinv
 * minor revisions
 *
 * Revision 1.6  2005/12/05  04:54:25  avi
 * maximum lag output increased by 1 (for same specified -m<int>)
 * include cross-spectral code but not yet output
 *
 * Revision 1.5  2005/05/19  04:05:22  avi
 * change meaning of mlag (make new -m32 equivalent to old -m33)
 *
 * Revision 1.4  2005/02/26  04:33:32  avi
 * MAXS (profile input line) -> 1024 chars
 *
 * Revision 1.3  2005/01/27  06:52:29  avi
 * -q (quiet mode)
 * -a (write ACV for all ROIs in one file)
 *
 * Revision 1.2  2005/01/25  07:08:18  avi
 * compute lagged covariance
 *
 * Revision 1.1  2005/01/22  02:09:29  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>				/* getpid () */
#include <assert.h>
#include <endianio.h>
#include <librms.h>

#define MAXR		64			/* maximum ROI label length */
#define MAXL		256			/* maximum filnename characters */
#define MAXC		8192			/* maximum number of profile fields */
#define MAXS		MAXC*10			/* maximum length of profile input line */
#define MLAG		32			/* default maximum lag */
#define DELTA		1.0			/* default delta */

typedef struct {
	float	r;
	float	i;
} CMPLX;

/*************/
/* externals */
/*************/
extern int	expandf (char *string, int len);						/* expandf.c */

static int split (char *string, char *srgv[], int maxp) {
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

static void getrootl (char *filespc, char *filroot) {
	char	*str;
	strcpy (filroot, filespc);
	while (str = strrchr (filroot, '.')) {
			if (!strcmp (str, ".rec"))	*str = '\0';
		else	if (!strcmp (str, ".dat"))	*str = '\0';
		else	if (!strcmp (str, ".txt"))	*str = '\0';
		else	break;
	}
}

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: covariance.c,v 1.31 2012/05/29 04:26:17 avi Exp $";
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

CMPLX **calloc_cmplx2 (int n1, int n2) {
	int	i;
	CMPLX	**a;

	if (!(a = (CMPLX **) malloc (n1 * sizeof (CMPLX *)))) errm (program);
	if (!(a[0] = (CMPLX *) calloc (n1 * n2, sizeof (CMPLX)))) errm (program);
	for (i = 1; i < n1; i++) a[i] = a[0] + i*n2;
	return a;
}

void free_cmplx2 (CMPLX **a) {
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

char **calloc_char2 (int n1, int n2) {
	int	i;
	char	**a;

	if (!(a = (char **) malloc (n1 * sizeof (char *)))) errm (program);
	if (!(a[0] = (char *) calloc (n1 * n2, sizeof (char)))) errm (program);
	for (i = 0; i < n1; i++) a[i] = a[0] + i*n2;
	return a;
}

void free_char2 (char **a) {
	free (a[0]);
	free (a);
}

double zero_mean (float *f, int npts, char *format) {
	int		i, n;
	double		u;

	for (u = n = i = 0; i < npts; i++) if (format[i] != 'x') {n++; u += f[i];}
	u /= n;
	for (i = 0; i < npts; i++) f[i] -= u;
	return u;
}

double trend_out (float *f, int npts, char *format) {
	int		i;
	double		x, sxx, sxy, beta;

	for (sxy = sxx = i = 0; i < npts; i++) if (format[i] != 'x') {
		x = -1. + 2.*i/(npts - 1);
		sxx += x*x;
		sxy += x*f[i];
	}
	beta = sxy/sxx; 
	for (i = 0; i < npts; i++) {
		x = -1. + 2.*i/(npts - 1);
		f[i] -= beta*x;
	}
	return 2.*beta/(npts - 1);
}

double unit_var (float *f, int npts, char *format) {
	int		i, n;
	double		v;

	for (v = n = i = 0; i < npts; i++) if (format[i] != 'x') {n++; v += f[i]*f[i];}
	v /= n;
	for (i = 0; i < npts; i++) {
		if (format[i] == 'x') {
			f[i] = 0.0;
		} else {
			f[i] /= sqrt (v);
		}
	}
	return v;
}

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n#%s\n", rcsid);
}

void usage (char *program) {
	fprintf (stderr, "Usage:\t%s <format> <profile>\n", program);
	fprintf (stderr, " e.g.,\t%s \"4x124+\" doubletask.txt\n", program);
	fprintf (stderr, "\toption\n");
	fprintf (stderr, "\t-q\tquiet mode\n");
	fprintf (stderr, "\t-t\toptionally remove trend (ramp) from input timeseries\n");
	fprintf (stderr, "\t-u\toptionally normalize all input timeseries to unit variance\n");
	fprintf (stderr, "\t-o\toutput lagged CCV dat files (CCR with -u)\n");
	fprintf (stderr, "\t-a\toutput lagged ACV dat file  (ACR with -u)\n");
	fprintf (stderr, "\t-r\toutput Bartlett smoothed cross spectra (spectral density with -u)\n");
	fprintf (stderr, "\t-p\toutput Bartlett smoothed auto  spectra (spectral density with -u)\n");
	fprintf (stderr, "\t-e\tcompute eigenvectors of lag 0 CCV\n");
	fprintf (stderr, "\t-L\tread ROI labels from <profile> (default ignore \'#\' commented lines)\n");
	fprintf (stderr, "\t-T<int>\tadditionally smooth spectra with Tukey window of specified width (in frames)\n");
	fprintf (stderr, "\t-d<flt>\tspecify frame TR in sec for Fourier analysis (default = %.4f)\n", DELTA);
	fprintf (stderr, "\t-m<int>\tspecify CCV function maxumum lag in frames (default = %d)\n", MLAG);
	fprintf (stderr, "\t-D<flt>\tSVD lag 0 CCV and output new profile with cndnum < specified value\n");
	fprintf (stderr, "\t-g<fil>\tregress timeseries in named file out of <profile>\n");
	fprintf (stderr, "N.B.:\tall input timeseries are made zero mean as a first step\n");
	fprintf (stderr, "N.B.:\tregion names can be specified on the first line of <profile> with '#' in the first column\n");
	fprintf (stderr, "N.B.:\toption -D inplies -e\n");
	exit (1);
}

int main (int argc, char *argv[]) {
/**********************/
/* filename variables */
/**********************/
	FILE		*tmpfp, *imgfp, *outfp, *profp;
	char		outroot[MAXL], outfile[MAXL], tmpfile[MAXL], profile[MAXL], proroot[MAXL];
	char		**ROI;

/************************/
/* auxiliary timeseries */
/************************/
	char		profil2[MAXL] = "";
	int		gcol, gpts;
 	double		ddet, **g, **gtg, **ginvt, **beta;

/*************************/
/* timeseries processing */
/*************************/
	char		*format;
	double		*sigma, *mean, *slope;
	float		*f, **C, *a, *b;
	int		ppts, npts, nnez, ncol, mlag = MLAG;

/*********************/
/* spectral analysis */
/*********************/
	int		one = 1, negone = -1, mlag2;
	int		Tukeym = 0;	/* Tukey window maximum lag */
	float		*Tukeyw;	/* Tukey window smoothing factor J&W p252 */
	float		delta = DELTA;	/* sampling interval in sec (frame TR) */
	float		**F;		/* F[lag][icol] spectral density array */
	float		*power;		/* power[icol] */
	double		p, q, tau, amp, phi, phi0, coh, delrad;
	CMPLX		**G;		/* FFT cross spectral array */

/*******/
/* SVD */
/*******/
	double		*P, *W, cndcrit = 0.0;

/***********/
/* utility */
/***********/
	char		command[MAXL], string[MAXS+4], *ptr, *srgv[MAXC];
	int		c, i, j, k, l, m;

/*********/
/* flags */
/*********/
	int		status = 0;
	int		unitvar	= 0;
	int		trendout = 0;
	int		diagonalize = 0;
	int		write_CRS = 0;
	int		write_CCV = 0;
	int		write_ACV = 0;
	int		write_autospec = 0;
	int		dospec;
	int		readlabel = 0;
	int		quiet = 0;
	int		debug = 0;

	fprintf (stderr, "%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);

	if (!(format = (char *) calloc (1024, sizeof (char)))) errm (program);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'q': quiet++;			break;
				case 'L': readlabel++;			break;
				case 'r': write_CRS++;			break;
				case 'o': write_CCV++;			break;
				case 'a': write_ACV++;			break;
				case 'p': write_autospec++;		break;
				case 'e': diagonalize++;		break;
				case 't': trendout++;			break;
				case 'u': unitvar++;			break;
				case 'd': delta = atof (ptr);		break;
				case 'g': strcpy (profil2, ptr);		*ptr = '\0'; break;
				case 'm': mlag = atoi (ptr);			*ptr = '\0'; break;
				case 'T': Tukeym = atoi (ptr);			*ptr = '\0'; break;
				case 'D': cndcrit = atof (ptr); diagonalize++;	*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	strncpy (format, argv[i], 1023);	k++; break;
			case 1:	strcpy (profile, argv[i]);		k++; break;
		}	
	}
	if (mlag < 0) mlag = 0;
	dospec = write_autospec || write_CRS;
	printf ("format=%s\n", format);
	if (k < 2) usage (program);

/****************/
/* read profile */
/****************/
	getrootl (profile, proroot);
	fprintf (stderr, "Reading: %s\n", profile);
	if (!(profp = fopen (profile, "r"))) errr (program, profile);
	ppts = 0; while (fgets (string, MAXS + 1, profp)) {
		if (!(m = split (string, srgv, MAXC))) continue;
		if (!ppts) {
			ncol = m;
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
	if (k = expandf (format, ppts + 1)) {
		fprintf (stderr, "%s: format defined for more points than %s\n", program, profile);
		exit (-1);
	}
	if (debug) printf ("%s\n", format);
	npts = strlen (format);
	for (nnez = k = 0; k < npts; k++) if (format[k] != 'x') nnez++;
	if (ppts > npts) {
		fprintf (stderr, "%s warning: format defined for less points than %s\n", program, profile);
	}
	fprintf (stderr, "npts=%d ncol=%d\n", npts, ncol);
	if (!(mean =	(double *) malloc (       ncol * sizeof (double)))) errm (program);
	if (!(slope =	(double *) malloc (       ncol * sizeof (double)))) errm (program);
	if (!(sigma =	(double *) malloc (       ncol * sizeof (double)))) errm (program);
	if (!(f =	(float *)  malloc (npts * ncol * sizeof (float))))  errm (program);
	C = calloc_float2 (mlag + 1, ncol * ncol);

if (0) {
	if (!(a =	(float *)  malloc (   2 * mlag * sizeof (float))))  errm (program);
	if (!(b =	(float *)  malloc (   2 * mlag * sizeof (float))))  errm (program);
	G = calloc_cmplx2 (mlag + 1, ncol * ncol);
}

/********************/
/* set region names */
/********************/
	ROI = calloc_char2 (ncol, MAXR + 1); for (j = 0; j < ncol; j++) sprintf (ROI[j], "ROI%d", j + 1);
	rewind (profp);
	fgets (string, MAXS + 1, profp);
	if (readlabel && string[0] == '#') {
		ptr = string; while (ptr = strpbrk (ptr, "#")) *ptr = '\x20'; /* convert '#' to space */
		if (m = split (string, srgv, MAXC)) {
			if (m != ncol) {
				fprintf (stderr, "%s warning: number of ROI labels does not match number of columns\n",
					program);
				if (m > ncol) m = ncol;
			}
			for (j = 0; j < m; j++) {
				if (strlen (srgv[j]) > MAXR) {
					fprintf (stderr, "% warning: ROI label %s truncated in output filenames\n",
					program, srgv[j]);
				}
				strncpy (ROI[j], srgv[j], MAXR);
			}
		}	
	}
	if (!quiet) for (j = 0; j < ncol; j++) fprintf (stderr, "ROI%d\t%s\n", j + 1, ROI[j]);

	rewind (profp);
	i = 0; while (i < npts) {
		fgets (string, MAXS + 1, profp);
		if (!(m = split (string, srgv, MAXC))) continue;
		for (j = 0; j < ncol; j++) (f + j*npts)[i] = atof (srgv[j]);
		i++;
	}
	fclose (profp);
	assert (i == npts);

/******************************************/
/* always make input timeseries zero mean */
/******************************************/
	for (j = 0; j < ncol; j++) mean[j] = zero_mean (f + j*npts, npts, format);
/**********************************/
/* optionally remove trend (ramp) */
/**********************************/
	if (trendout) for (j = 0; j < ncol; j++) {
		slope[j] = trend_out (f + j*npts, npts, format);
		if (debug) printf ("%10d% 9.4f\n", j + 1, slope[j]);
		if (0) for (i = 0; i < npts; i++) printf (" %10.6f", (f + j*npts)[i]);
	}
/*********************************************/
/* optionally impose unit variance condition */
/*********************************************/
	if (unitvar) for (j = 0; j < ncol; j++) {
		sigma[j] = unit_var (f + j*npts, npts, format);
		if (debug) printf ("%10d% 9.3f\n", j + 1, sigma[j]);
	}

	if (strlen (profil2)) {
/***********************/
/* read second profile */
/***********************/
		fprintf (stderr, "Reading: %s\n", profil2);
		if (!(profp = fopen (profil2, "r"))) errr (program, profil2);
		gpts = 0; while (fgets (string, MAXS + 1, profp)) {
			if (!(m = split (string, srgv, MAXC))) continue;
			if (!gpts) {
				gcol = m;
			} else {
				if (m != gcol) {
					fprintf (stderr, "%s: %s format error\n", program, profil2);
					exit (-1);
				}
			}
			gpts++;
		}
		if (gpts != ppts) {
			fprintf (stderr, "%s: %s %s length mismatch\n", program, profile, profil2);
			exit (-1);
		}
		g = calloc_double2 (gcol, nnez);
		rewind (profp);
		k = i = 0; while (i < npts) {
			fgets (string, MAXS + 1, profp);
			if (!(m = split (string, srgv, MAXC))) continue;
			if (format[i++] == 'x') continue;
			for (j = 0; j < gcol; j++) g[j][k] = atof (srgv[j]);
			k++;
		}
		fclose (profp);
		assert (i == npts);
		assert (k == nnez);

/*************************************************************/
/* always make auxilliary timeseries zero mean unit variance */
/*************************************************************/
		for (j = 0; j < gcol; j++) {
			for (q = i = 0; i < nnez; i++) q += g[j][i];
			q /= nnez;
			for (    i = 0; i < nnez; i++) g[j][i] -= q;
		}
		for (j = 0; j < gcol; j++) {
			for (q = i = 0; i < nnez; i++) q += g[j][i]*g[j][i];
			q /= nnez;
			for (    i = 0; i < nnez; i++) g[j][i] /= sqrt (q);
		}

/***********************************************/
/* compute gcol x ncol regression coefficients */
/***********************************************/
		gtg = calloc_double2 (gcol, gcol);
		for (i = 0; i < gcol; i++) for (j = 0; j < gcol; j++) {
			for (k = 0; k < nnez; k++) gtg[j][i] += g[i][k]*g[j][k];
			gtg[j][i] /= nnez;
		}
		dmatinv_ (gtg[0], &gcol, &ddet);
		ginvt = calloc_double2 (gcol, nnez);
		for (i = 0; i < gcol; i++) for (j = 0; j < nnez; j++) {
			for (k = 0; k < gcol; k++) ginvt[i][j] += gtg[i][k]*g[k][j];
			
		}
		beta = calloc_double2 (ncol, gcol);
		for (m = 0; m < ncol; m++) for (l = 0; l < gcol; l++) {
			for (k = i = 0; i < npts; i++) if (format[i] != 'x') {
				beta[m][l] += ginvt[l][k]*(f + m*npts)[i];
				k++;
			}
			assert (k == nnez);
			beta[m][l] /= nnez;
		}

if (0 && gcol == 1) {
	printf ("gtg[0][0]=%f\n", gtg[0][0]);
	printf ("ginvt[0][k]");
	for (k = 0; k < 5; k++) printf (" %10.6f", ginvt[0][k]); printf ("\n");
	printf ("f[k]       ");
	for (k = 0; k < 5; k++) printf (" %10.6f", f[k]); printf ("\n");
	for (m = 0; m < ncol; m++) printf ("beta[%d][0]=%10.6f\n", m, beta[m][0]);
}

/**************************************************/
/* regress auxilliary profile out of main profile */
/**************************************************/
		for (m = 0; m < ncol; m++) {
			for (k = i = 0; i < npts; i++) if (format[i] != 'x') {
				for (l = 0; l < gcol; l++) (f + m*npts)[i] -= g[l][k]*beta[m][l];
				k++;
			}
			assert (k == nnez);
		}
		free_double2 (g); free_double2 (gtg); free_double2 (ginvt); free_double2 (beta);

/***************************/
/* write regressed profile */
/***************************/
		sprintf (outfile, "%s_regr.dat", proroot);
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		fprintf (stderr, "Writing: %s\n", outfile);
		write_command_line (outfp, argc, argv);
		fprintf (outfp, "#%9s", ROI[0]);
		for (j = 1; j < ncol; j++) fprintf (outfp, "%10s", ROI[j]); fprintf (outfp, "\n");
		for (i = 0; i < npts; i++) {
			for (j = 0; j < ncol; j++) fprintf (outfp, " %9.4f", (f + j*npts)[i]);
			fprintf (outfp, "\n");
		}
		if (fclose (outfp)) errw (program, outfile);
	}	/* end auxilliary profile regression code */

/*************************************/
/* assemble lagged covariance matrix */
/*************************************/
	if (mlag >= npts) {
		mlag = npts - 1;
		fprintf (stderr, "%s: maximum lag adjusted to %d\n", program, mlag);
	}
	if (!quiet) printf ("computing covariance matrix lag");
	for (l = 0; l <= mlag; l++) {
		if (!quiet && !(l%1024)) {printf (" %d", l); fflush (stdout);}
		for (i = 0; i < ncol; i++) {
		for (j = 0; j < ncol; j++) {
			for (q = k = 0; k < npts - l; k++) {
				if (format[k] != 'x' && format[k + l] != 'x') q += (f + i*npts)[k] * (f + j*npts)[k + l];
			}
			C[l][i + j*ncol] = q/nnez;
			if (0) printf ("lag=%3d C[%d,%d]=%10.6f\n", l, i + 1, j + 1, C[l][i + j*ncol]);
		}}
	}
	if (!quiet) printf ("\n");

	if (dospec) {
		F = calloc_float2 (mlag + 1, ncol);
		m = (Tukeym > 0) ? Tukeym : mlag;
		if (m > mlag) {
			m = mlag;
			fprintf (stderr, "%s: Tukey width set to %d\n", program, mlag);
		}
		if (!(Tukeyw = (float *) calloc (m + 1, sizeof (float)))) errm (program);
/************************************************************/
/* prepare Tukey window and compute autospectrum J&W p. 310 */
/************************************************************/
		printf ("computing Tukey window"); fflush (stdout);
		for (k = 0; k <= m; k++) Tukeyw[k] = (Tukeym > 0) ? 0.5*(1. + cos(M_PI*k/m)) : 1.0;
		for (j = 0; j < ncol; j++) {
			for (i = 0; i <= m; i++) {
				q = C[0][j + j*ncol];
				for (k = 1; k <= m; k++) {
					q += 2*C[k][j + j*ncol]*Tukeyw[k]*cos(M_PI*k*i/m);
				}
				F[i][j] = 2*delta*q;
			}
		}
		printf ("\n");
	}

/**********************************************************************************/
/* compute cross spectral analysis (m = (Tukeym > 0) ? Tukeym : mlag;) J&W p. 419 */
/**********************************************************************************/
	if (write_CRS) {
		m = (Tukeym > 0) ? Tukeym : mlag;	/* Tukey window previously prepared */
		for (i =     0; i < ncol; i++) {
		for (j = i + 1; j < ncol; j++) {
			sprintf (outfile, "%s_%s_%s_CRS.dat", proroot, ROI[i], ROI[j]);
			if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
			fprintf (stderr, "Writing: %s\n", outfile);
			write_command_line (outfp, argc, argv);
			fprintf (outfp, "#%9s%10s%10s%10s%10s%10s%10s\n",
				"freq (Hz)", "cospec", "quadspec", "amplitude", "coherence", "phase", "lag (sec)");

			for (l = 0; l <= m; l++) {
				q = C[0][i + j*ncol];
				p = 0.;
				for (k = 1; k <= m; k++) {
					q += (C[k][i + j*ncol]+C[k][j + i*ncol])*Tukeyw[k]*cos(M_PI*k*l/m);
					p += (C[k][i + j*ncol]-C[k][j + i*ncol])*Tukeyw[k]*sin(M_PI*k*l/m);
				}
				amp = 2*delta*sqrt(q*q + p*p);
				coh = amp/sqrt(F[l][i]*F[l][j]);	/* autospectra previously computed */
				phi0 = phi;
				phi = atan2 (-p, q);			/* negative values indicate i precedes j */
				delrad = phi - phi0;
				if (fabs (delrad) > fabs (delrad - 2.*M_PI)) delrad -= 2.*M_PI;
				if (fabs (delrad) > fabs (delrad + 2.*M_PI)) delrad += 2.*M_PI;
				tau = (l > 0) ? m*delta*delrad/M_PI : 0.;
				fprintf (outfp, "%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n", l/(2*m*delta),
					2*delta*q, 2*delta*p, amp, coh, phi, tau);
			}

			if (fclose (outfp)) errw (program, outfile);
		}}
	}

/********************************/
/* compute cross spectra by FFT */
/********************************/
	if (0) {
		mlag2 = mlag*2;
		for (i = 0; i < ncol; i++) {
		for (j = 0; j < ncol; j++) {
			for (l = 0; l <  mlag; l++) {a[l] =          C[l][i + j*ncol]; b[l] =          0.0;}
			for (l = 1; l <= mlag; l++) {a[2*mlag - l] = C[l][j + i*ncol]; b[2*mlag - l] = 0.0;}
			fft_ (a, b, &one, &mlag2, &one, &negone);
			for (l = 0; l <= mlag; l++) {
				G[l][i + j*ncol].r = a[l]/mlag2;
				G[l][i + j*ncol].i = b[l]/mlag2;
			}
		}}
/**************************************/
/* output cross spectral FFT analysis */
/**************************************/
		for (i =     0; i < ncol; i++) {
		for (j = i + 0; j < ncol; j++) {
			sprintf (outfile, "%s_%s_%s_CRS.dat", proroot, ROI[i], ROI[j]);
			if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
			fprintf (stderr, "Writing: %s\n", outfile);
			write_command_line (outfp, argc, argv);
			fprintf (outfp, "#%9s%10s%10s%10s%10s%10s\n",
				"freq (Hz)", "cospec", "quadspec", "amplitude", "pha (sec)", "coherence");
			for (l = 0; l <= mlag; l++) {
				amp = sqrt (G[l][i + j*ncol].r*G[l][i + j*ncol].r + G[l][i + j*ncol].i*G[l][i + j*ncol].i);
				phi = atan2 (-G[l][i + j*ncol].i, G[l][i + j*ncol].r);
				coh = amp/sqrt(G[l][i + i*ncol].r*G[l][j + j*ncol].r);
				if (l == 0) {amp = phi = 0.; coh = 1.;}
				fprintf (outfp, "%10.6f%10.6f%10.6f%10.6f%10.4f%10.4f\n", l/(mlag2*delta),
					G[l][i + j*ncol].r, G[l][i + j*ncol].i, amp, delta*phi*mlag/M_PI, coh);
			}
			if (fclose (outfp)) errw (program, outfile);
		}}
	}

/********************************************/
/* write lagged covariance timecourse files */
/********************************************/
	if (write_CCV) {
		for (j = 0; j < ncol; j++) {
			sprintf (outfile, "%s_%s_CC%c.dat", proroot, ROI[j], (unitvar) ? 'R' : 'V');
			if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
			fprintf (stderr, "Writing: %s\n", outfile);
			fprintf (outfp, "#%s", program);
			for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]); fprintf (outfp, "\n");
			fprintf (outfp, "#%s\n", rcsid);
			fprintf (outfp, "#%9s", "lag");
			for (i = 0; i < ncol; i++) fprintf (outfp, "%10s", ROI[i]); fprintf (outfp, "\n");
			for (l = mlag; l > 0; l--) {
				fprintf (outfp, "%10d", -l);
				for (i = 0; i < ncol; i++) fprintf (outfp, ((unitvar) ? " %9.6f" : " %9.3f"), C[l][i + j*ncol]);
				fprintf (outfp, "\n");
			}
			for (l = 0; l <= mlag; l++) {
				fprintf (outfp, "%10d", l);
				for (i = 0; i < ncol; i++) fprintf (outfp, ((unitvar) ? " %9.6f" : " %9.3f"), C[l][j + i*ncol]);
				fprintf (outfp, "\n");
			}
			if (fclose (outfp)) errw (program, outfile);
		}
		if (mlag == 0) {
			sprintf (outfile, "%s_CC%c.dat", proroot, (unitvar) ? 'R' : 'V');
			fprintf (stderr, "Writing: %s\n", outfile);
			if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
			write_command_line (outfp, argc, argv);
			for (j = 0; j < ncol; j++) {
				for (i = 0; i < ncol; i++) fprintf (outfp, ((unitvar) ? " %9.6f" : " %9.3f"), C[0][j + i*ncol]);
				fprintf (outfp, "\n");
			}
			if (fclose (outfp)) errw (program, outfile);
		}
	}
	if (write_ACV) {
		sprintf (outfile, "%s_AC%c.dat", proroot, (unitvar) ? 'R' : 'V');
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		fprintf (stderr, "Writing: %s\n", outfile);
		fprintf (outfp, "#%s", program);
		for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]); fprintf (outfp, "\n");
		fprintf (outfp, "#%s\n", rcsid);
		fprintf (outfp, "#%9s", "lag");
		for (i = 0; i < ncol; i++) fprintf (outfp, "%10s", ROI[i]); fprintf (outfp, "\n");
		for (l = mlag; l > 0; l--) {
			fprintf (outfp, "%10d", -l);
			for (i = 0; i < ncol; i++) fprintf (outfp, ((unitvar) ? " %9.6f" : " %9.3f"), C[l][i + i*ncol]);
			fprintf (outfp, "\n");
		}
		for (l = 0; l <= mlag; l++) {
			fprintf (outfp, "%10d", l);
			for (i = 0; i < ncol; i++) fprintf (outfp, ((unitvar) ? " %9.6f" : " %9.3f"), C[l][i + i*ncol]);
			fprintf (outfp, "\n");
		}
		if (fclose (outfp)) errw (program, outfile);
	}

/***************************************************************/
/* write autospectrum file (m = (Tukeym > 0) ? Tukeym : mlag;) */
/***************************************************************/
	if (write_autospec) {
		if (!(power = (float *) calloc (ncol, sizeof (float)))) errm (program);
		sprintf (outfile, "%s_%s.dat", proroot, (unitvar) ? "PSD" : "PWR");
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		fprintf (stderr, "Writing: %s\n", outfile);
		fprintf (outfp, "#%s", program);
		for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]); fprintf (outfp, "\n");
		fprintf (outfp, "#%s\n", rcsid);
		fprintf (outfp, "#%9s", "Hz");
		for (i = 0; i < ncol; i++) fprintf (outfp, "%10s", ROI[i]); fprintf (outfp, "\n");
		for (j = 0; j <= m; j++) {
			fprintf (outfp, "%10.4f", j/(2*m*delta));
			for (i = 0; i < ncol; i++) {
				fprintf (outfp, ((unitvar) ? " %9.6f" : " %9.3f"), F[j][i]);
				power[i] += F[j][i];
			}
			fprintf (outfp, "\n");
		}
		for (i = 0; i < ncol; i++) power[i] /= 2*m*delta;
		fprintf (outfp, "#%-9s", "tot power");
		for (i = 0; i < ncol; i++) fprintf (outfp, ((unitvar) ? " %9.6f" : " %9.3f"), power[i]);
		fprintf (outfp, "\n");
		if (fclose (outfp)) errw (program, outfile);
		free (power);
	}

/******************************************/
/* diagonalize zero lag covariance matrix */
/******************************************/
	if (diagonalize) {
		if (!(P = (double *)  malloc (ncol * ncol * sizeof (double))))  errm (program);
		if (!(W = (double *)  malloc (ncol * ncol * sizeof (double))))  errm (program);
		for (i = 0; i < ncol*ncol; i++) P[i] = C[0][i];
		deigen_ (P, W, &ncol);
		printf ("eigenvectors:\n");
		for (j = 0; j < ncol; j++) {
			for (i = 0; i < ncol; i++) printf ("%10.6f", W[j*ncol + i]);
			printf (" eigval=%10.6f\n", P[j*ncol + j]);
		}
		printf ("condition number = %.4e\n", P[0]/P[ncol*ncol - 1]);

		if (cndcrit > 0.0) {
			m = 0; while (m < ncol) {
				if (P[0]/P[m*(ncol + 1)] > cndcrit) break;
				m++;
			}
			sprintf (outfile, "%s_SVD%d.dat", proroot, m);
			if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
			fprintf (stderr, "Writing: %s\n", outfile);
			if (0) write_command_line (outfp, argc, argv);
			for (i = 0; i < npts; i++) {
				for (k = 0; k < m; k++) {
					for (q = j = 0; j < ncol; j++) q += W[k*ncol + j]*(f + j*npts)[i];
					fprintf (outfp, "%9.4f ", q);
				}
				fprintf (outfp, "\n");
			}
			if (fclose (outfp)) errw (program, outfile);
		}
		free (P); free (W);
	}

	free (format);
	free_char2 (ROI);
	free (f); free (sigma); free (mean); free (slope);
	free_float2 (C);
	if (dospec) {free_float2 (F); free (Tukeyw);}
	if (0) {free (a); free (b); free_cmplx2 (G);}
	exit (status);
}



