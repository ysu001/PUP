/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/GC_dat.c,v 1.10 2007/07/30 04:00:45 avi Exp $*/
/*$Log: GC_dat.c,v $
 * Revision 1.10  2007/07/30  04:00:45  avi
 * Solaris 10
 *
 * Revision 1.9  2007/07/15  03:02:37  avi
 * correct normalization of covariance (nnez -> ntot)
 * write to files instead of stdout
 *
 * Revision 1.8  2006/11/29  02:09:43  avi
 * correct unit_var() bug
 *
 * Revision 1.7  2006/08/14  04:54:10  avi
 * remove /sqrt from "N(0,2)" result
 *
 * Revision 1.6  2006/08/08  03:41:06  avi
 * output Fx->y - Fy->x and Geweke's N(0,1) measure
 *
 * Revision 1.5  2006/07/30  04:56:11  avi
 * reorganize computation of post-VAR variance at zero lag (sigma)
 *
 * Revision 1.4  2006/07/30  03:26:03  avi
 * correct VAR formulation - works correctly
 *
 * Revision 1.3  2006/07/27  01:28:15  avi
 * compute GC indices
 * ? bug: GC indices may go slightly negative
 *
 * Revision 1.2  2006/07/26  04:14:39  avi
 * subroutines compute_covariance() and compute_AR()
 *
 * Revision 1.1  2006/07/26  02:24:06  avi
 * Initial revision
 **/

#include	<stdlib.h>
#include	<stdio.h>
#include	<math.h>
#include	<stdlib.h>
#include	<unistd.h>		/* R_OK */
#include	<string.h>
#include	<ctype.h>
#include	<assert.h>

#define MAXF		16384			/* maximum number of input timepoints */
#define MAXL		256
#define MAXS		4096			/* maximum length of profile input line */

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: GC_dat.c,v 1.10 2007/07/30 04:00:45 avi Exp $";
static char	program[MAXL];
static int	debug = 0;
static int	verbose = 0;

void compute_covariance (int npts, int ncol, int iorder, char *format, float *f, float ***G);	/* below */
void compute_AR (int npts, int ncol, int iorder, char* format, float *f, float ***A, float ***G, float **S);	/* below */

/*************/
/* externals */
/*************/
extern int	expandf (char *string, int len);				/* expandf.c */
extern void	format_test_ (char *format, int *npts);				/* fGC.f */
extern void	solve_var_  (int *ncol, int *iorder, float *G, float *cc, float *ca, float *cb, float *A, float *det, float *S);
extern void	apply_var_  (int *ncol, int *iorder, int *npts, float *A, char *format, float *f, float *z);
extern void	verify_var_ (int *ncol, int *iorder, float *G, float *A);
extern void	compute_gc_ (int *dimx, int *dimy, float *Sx, float *Sy, float *Sq);
extern void	compute_gc1_ (int *dimx, int *dimy, float *Sx, float *Sy, float *Sq,
					float *cx, float *cy, float *cq, float *F);

void errm (char* program) {
	fprintf (stderr, "%s: memory allocation error\n", program);
	exit (-1);
}

void errr (char* program, char* filespc) {
	fprintf (stderr, "%s: %s read error\n", program, filespc);
	exit (-1);
}

void errw (char* program, char* filespc) {
	fprintf (stderr, "%s: %s write error\n", program, filespc);
	exit (-1);
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

double zero_mean (float *f, int npts, char *format) {
	int		i, n;
	double		u;

	for (u = n = i = 0; i < npts; i++) if (format[i] != 'x') {n++; u += f[i];}
	u /= n;
	for (i = 0; i < npts; i++) f[i] -= u;
	return u;
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

double rectifier (double x) {
	return (x > 0.0) ? x : 0.0;
}

void getrootl (char *filespc, char *filroot) {
	char	*str;
	strcpy (filroot, filespc);
	while (str = strrchr (filroot, '.')) {
			if (!strcmp (str, ".rec"))	*str = '\0';
		else	if (!strcmp (str, ".dat"))	*str = '\0';
		else	if (!strcmp (str, ".txt"))	*str = '\0';
		else	break;
	}
}

void write_command_line (FILE *outfp, int argc, char **argv) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
}

void usage (char* program) {
	printf ("Usage:\t%s <format> <input_datafile> <order>\n", program);
	printf ("e.g.,\t%s \"4x106+\" ROI_timeseries.dat 2\n", program);
	printf ("\toption\n");
	printf ("\t-d\tdebug mode\n");
	printf ("\t-v\tverbose mode\n");
	printf ("\t-u\tnormalize all input timeseries to unit variance\n");
	printf ("\t-x<int>\tspecify dimensionality of x process (default = 1)\n");
	printf ("\t-m\tcreate text listing of AR model\n");
	printf ("\t-w\twrite residual after full AR modeling\n");
	printf ("\t-P\tformat residual output suitable for plotting (xyy)\n");
	exit (1);
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

int main (int argc, char **argv) {
	char		datfile[MAXL], profile[MAXL], proroot[MAXL], outfile[MAXL];
	FILE		*outfp, *profp;

/*************************/
/* timeseries processing */
/*************************/
	char		format[MAXF], *formax;
	double		*sigma, *mean;
	float		*f, *z;
	float		***Gq, ***Aq;			/* full model gamma and AR coefficients */
	float		***Gx, ***Ax, ***Gy, ***Ay;	/* x and y    gamma and AR coefficients */
	float		**Sq, **Sx, **Sy;		/* omega in Yule-Walker equations */
	float		*xx, *yy, *qq, F[4];		/* scratch arrays for compute_gc1() */
	int		iorder;				/* order of AR model */
	int		npts, nnez, ncol, dimx = 1, dimy;
	double		Gconst, q;

/***********/
/* utility */
/***********/
	char		*ptr, string[MAXL], command[MAXL], *srgv[MAXL];
	int		c, i, j, k, l, m;

/*********/
/* flags */
/*********/
	int		plot = 0;
	int		wflag = 0;
	int		mflag = 0;
	int		unitvar	= 0;
	int		status = 0;

/*	f_init ();*/
	printf ("%s\n", rcsid);
	setprog (program, argv);
/******************************/
/* get command line arguments */
/******************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'u': unitvar++;		break;
				case 'v': verbose++;		break;
				case 'd': debug++;		break;
				case 'w': wflag++;		break;
				case 'm': mflag++;		break;
				case 'P': plot++;		break;
				case 'x': dimx = atoi (ptr);	*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0:	strcpy (format, argv[i]);	k++; break;
			case 1:	strcpy (profile, argv[i]);	k++; break;
			case 2:	iorder = atoi (argv[i]);	k++; break;
		}
	}
	if (k < 3) usage (program);

/****************/
/* parse format */
/****************/
	if (k = expandf (format, MAXF)) exit (k);
	if (debug) printf ("%s\n", format);
	npts = strlen (format);
	for (nnez = k = 0; k < npts; k++) if (format[k] != 'x') nnez++;
	printf ("%s: time series defined for %d frames, %d exluded\n", program, npts, npts - nnez);
        if (0) format_test_ (format, &npts);

/*****************/
/* parse profile */
/*****************/
	printf ("Reading: %s\n", profile);
	if (!(profp = fopen (profile, "r"))) errr (program, profile);
	k = 0; while (fgets (string, MAXS, profp)) {
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
	dimy = ncol - dimx;
	if (dimx < 1 || dimy < 1) {
		fprintf (stderr, "(%s: illegal process dimensionality\n", program);
		usage (program);
	}
	if (k < npts) {
		fprintf (stderr, "%s: format defined for more points than %s\n", program, profile);
		exit (-1);
	}
	if (k > npts) {
		fprintf (stderr, "%s warning: format defined for less points than %s\n", program, profile);
	}

/*******************/
/* allocate memory */
/*******************/
	if (!(formax =	(char *)   calloc (npts + 1,     sizeof (char))))   errm (program);
	if (!(f =	(float *)  malloc (npts * ncol * sizeof (float))))  errm (program);
	if (!(z =	(float *)  malloc (npts * ncol * sizeof (float))))  errm (program);
	if (!(mean =	(double *) malloc (       ncol * sizeof (double)))) errm (program);
	if (!(sigma =	(double *) malloc (       ncol * sizeof (double)))) errm (program);
	Gq = calloc_float3 (iorder + 1, ncol, ncol);	/* lagged cross covariance (gamma) */
	Aq = calloc_float3 (iorder    , ncol, ncol);	/* AR coefficients */
	Sq = calloc_float2 (            ncol, ncol);	/* omega */
	Gx = calloc_float3 (iorder + 1, dimx, dimx);
	Ax = calloc_float3 (iorder    , dimx, dimx);
	Sx = calloc_float2 (            dimx, dimx);
	Gy = calloc_float3 (iorder + 1, dimy, dimy);
	Ay = calloc_float3 (iorder    , dimy, dimy);
	Sy = calloc_float2 (            dimy, dimy);
	if (!(xx = (float *) malloc (dimx * dimx * sizeof (float)))) errm (program);	/* compute_gc1_() scratch space */
	if (!(yy = (float *) malloc (dimy * dimy * sizeof (float)))) errm (program);
	if (!(qq = (float *) malloc (ncol * ncol * sizeof (float)))) errm (program);

/*************/
/* read data */
/*************/
	rewind (profp);
	i = 0; while (i < npts) {
		fgets (string, MAXS, profp);
		if (!(m = split (string, srgv, MAXL))) continue;
		for (j = 0; j < ncol; j++) (f + j*npts)[i] = atof (srgv[j]);
		i++;
	}
	fclose (profp);
	assert (i == npts);

/******************************************/
/* always make input timeseries zero mean */
/******************************************/
	for (j = 0; j < ncol; j++) mean[j] = zero_mean (f + j*npts, npts, format);

/*********************************************/
/* optionally impose unit variance condition */
/*********************************************/
	if (unitvar) for (j = 0; j < ncol; j++) sigma[j] = unit_var (f + j*npts, npts, format);

/*************************************************************/
/* compute format accounting for points excluded by AR order */
/*************************************************************/
	strcpy (formax, format);
	for (i = npts - 1; i >= 0; i--) {
		for (l = 1; l <= iorder; l++) {
			if (l > i) {
				formax[i] = 'x';
			} else {
				if (formax[i - l] == 'x') formax[i] = 'x';
			}
		}
	}
	if (verbose) printf ("%s\nformat post AR(%d)\n%s\n", format, iorder, formax);

/*****************************************************************/
/* compute and apply AR modeling to x, y, and combined processes */
/*****************************************************************/
	if (verbose) printf ("x dimensionality=%d\n", dimx);
	compute_AR (npts, dimx, iorder, format, f,             Ax, Gx, Sx);
	if (verbose) printf ("y dimensionality=%d\n", dimy);
	compute_AR (npts, dimy, iorder, format, f + dimx*npts, Ay, Gy, Sy);
	if (verbose) printf ("q dimensionality=%d\n", ncol);
	compute_AR (npts, ncol, iorder, format, f,             Aq, Gq, Sq);

	getrootl (profile, proroot);
	if (mflag) {
		sprintf (outfile, "%s_model.txt", proroot);
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		printf ("Writing: %s\n", outfile);
		write_command_line (outfp, argc, argv);
	} else {
		outfp = stdout;
	}
/*************************************/
/* compute Granger causality indices */
/*************************************/
/*	compute_gc_  (&dimx, &dimy, &Sx[0][0], &Sy[0][0], &Sq[0][0]);	*/
	compute_gc1_ (&dimx, &dimy, &Sx[0][0], &Sy[0][0], &Sq[0][0], xx, yy, qq, F);
	fprintf (outfp, "Fx,y            %10.6f\n", F[0]);
	fprintf (outfp, "Fx->y           %10.6f\n", F[1]);
	fprintf (outfp, "Fy->x           %10.6f\n", F[2]);
	fprintf (outfp, "Fx.y            %10.6f\n", F[3]);
	fprintf (outfp, "Fx->y - Fy->x   %10.6f\n", F[1] - F[2]);
/***************************************/
/* /sqrt(2.) converts N(0,2) to N(0,1) */
/***************************************/
	Gconst = (dimx*dimy*iorder - 1)/3.;	/* Geweke's (klp - 1)/3 */
/*	q = (sqrt (rectifier ((double) nnez*F[1] - Gconst))
	   - sqrt (rectifier ((double) nnez*F[2] - Gconst)))/sqrt(2.);	*/
	q = sqrt (nnez*F[1]) - sqrt (nnez*F[2]);
	fprintf (outfp, "Geweke's N(0,2) %10.6f\n", q);

/****************************************/
/* apply full model to input timeseries */
/****************************************/
	apply_var_ (&ncol, &iorder, &npts, &Aq[0][0][0], formax, f, z);

/**************************************/
/* write model to stdout or text file */
/**************************************/
	fprintf (outfp, "omega\n");
	for (i = 0; i < ncol; i++) {
		for (j = 0; j < ncol; j++) fprintf (outfp, "%10.4f", Sq[j][i]);
		fprintf (outfp, "\n");
	}
	for (l = 0; l < iorder; l++) {
		fprintf (outfp, "A(%d)\n", l + 1);
		for (i = 0; i < ncol; i++) {
			for (j = 0; j < ncol; j++) fprintf (outfp, "%10.6f", -Aq[l][j][i]);
			fprintf (outfp, "\n");
		}
	}
	if (mflag) fclose (outfp);

/******************************/
/* write whitened time series */
/******************************/
	if (wflag) {
		sprintf (outfile, "%s_prew.dat", proroot);
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		printf ("Writing: %s\n", outfile);
		write_command_line (outfp, argc, argv);
		for (i = 0; i < npts; i++) {
			if (plot) fprintf (outfp, "%10d", i);
			for (j = 0; j < ncol; j++) fprintf (outfp, "%10.4f", (z + j*npts)[i]);
			fprintf (outfp, "\n");
		}
		if (fclose (outfp)) errw (program, outfile);
	}

	free (formax);
	free (f); free (z); free (sigma); free (mean);
	free_float3 (Gq); free_float3 (Aq);
	free_float3 (Gx); free_float3 (Ax);
	free_float3 (Gy); free_float3 (Ay);
	free_float2 (Sq); free_float2 (Sx); free_float2 (Sy);
	free (xx); free (yy); free (qq);
/*	f_exit ();*/
	exit (0);
}

void compute_covariance (int npts, int ncol, int iorder, char *format, float *f, float ***G) {
	int		i, j, k, l, ntot;
	double		q;

	for (l = 0; l <= iorder; l++) {
		for (j = 0; j < ncol; j++) {
		for (i = 0; i < ncol; i++) {
			for (q = ntot = k = 0; k < npts - l; k++) if (format[k] != 'x' && format[k + l] != 'x') {
				q += (f + i*npts)[k + l] * (f + j*npts)[k];
				ntot++;
			}
			G[l][j][i] = q/ntot;
		}}
	}
	if (verbose) for (l = 0; l <= iorder; l++) {
		printf ("covariance (gamma) lag %d\n", l);
		for (i = 0; i < ncol; i++) {
			for (j = 0; j < ncol; j++) printf ("%10.4f", G[l][j][i]);
			printf ("\n");
		}
	}
}

void compute_AR (int npts, int ncol, int iorder, char* format, float *f, float ***A, float ***G, float **S) {
	float		*cc, *ca, *cb, det;		/* scratch arrays for AR solver */
	int		i, j, k;

	k = ncol*ncol*iorder;
	if (!(cc = (float *) malloc (k * k * sizeof (float)))) errm (program);
	if (!(ca = (float *) malloc (k     * sizeof (float)))) errm (program);
	if (!(cb = (float *) malloc (k     * sizeof (float)))) errm (program);
	if (verbose) printf ("before AR filtering\n");
	compute_covariance (npts, ncol, iorder, format, f, G);
	solve_var_             (&ncol, &iorder, &G[0][0][0], cc, ca, cb, &A[0][0][0], &det, &S[0][0]);
	if (debug) verify_var_ (&ncol, &iorder, &G[0][0][0], &A[0][0][0]);
	if (verbose) {
		printf ("det=%.4e\n", det);
		printf ("after AR(%d) filtering\n", iorder);
		printf ("omega\n");
		for (i = 0; i < ncol; i++) {
			for (j = 0; j < ncol; j++) printf ("%10.4f", S[j][i]);
			printf ("\n");
		}
	}
	free (ca); free (cb); free (cc);
}
