/*$Header: /data/petsun4/data1/src_solaris/qnt_4dfp/RCS/qntv_4dfp.c,v 1.6 2012/03/23 02:26:25 avi Exp $*/
/*$Log: qntv_4dfp.c,v $
 * Revision 1.6  2012/03/23  02:26:25  avi
 * improved Usage
 *
 * Revision 1.5  2011/10/04  22:43:53  avi
 * correct major voxel index / die index coding error
 *
 * Revision 1.4  2011/03/13  00:07:29  avi
 * correct computation of keyfile when input images are not in $cwd
 *
 * Revision 1.3  2011/03/10  04:01:45  avi
 * option -K
 *
 * Revision 1.2  2011/03/06  06:01:12  avi
 * correct O2 O3 output error
 *
 * Revision 1.1  2011/03/04  07:43:33  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>			/* getpid () */
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#ifndef _FEATURES_H		/* defined by math.h in Linux */
	#include <sunmath.h>	/* isnormal() */
	#include <ieeefp.h>
#endif
#include <endianio.h>
#include <conc.h>		/* includes definition of MAXL */

#define MIN(a,b)	((a) < (b) ? (a) : (b))
#define LENDIE		1
#define NCRIT		1
#define INCR		64	/* die voxel index memory increment size in items */
#define TOL		1.e-6	/* ratio of least to greatest eigenvalue */

/*************/
/* externals */
/*************/
extern void	 dsvdcmp0 (double **a, int m, int n, double *w, double **v);			/* dsvdcmp0.c */
extern int	ndsvdcmp0 (double **a, int m, int n, double *w, double **v, double tol);	/* dsvdcmp0.c */
extern int	expandf (char *string, int len);						/* expandf.c */

typedef struct {
	int	*p;		/* pointer to voxel indices */
	int	ma;		/* memory allocated for this many indices */
	int	n;		/* voxel count (ROI mask) defined */
	int	m;		/* voxel count (time series defined) */
	float	fx, fy, fz;	/* voxel coordinates */
	float	v;		/* value */
} DIEDAT;

/********************/
/* global variables */
/********************/
char	program[MAXL];
static char rcsid[]= "$Id: qntv_4dfp.c,v 1.6 2012/03/23 02:26:25 avi Exp $";

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

double dzeromean (double *f, int npts) {
	int		i;
	double		u;

	for (u = i = 0; i < npts; i++) u += f[i];
	u /= npts;
	for (i = 0; i < npts; i++) f[i] -= u;
	return u;
}

double dunitvar (double *f, int npts) {
	int		i;
	double		v;

	for (v = i = 0; i < npts; i++) v += f[i]*f[i];
	v /= npts;
	for (i = 0; i < npts; i++) f[i] /= sqrt (v);
	return v;
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
}

static void usage (char *program) {
	printf ("%s\n", rcsid);
	printf ("Usage:\t%s <(4dfp)|(conc) image> <(4dfp) ROI>\n", program);
	printf (" e.g.:\t%s TC30274_rmsp_faln_dbnd_xr3d_atl.conc iter10_roi_-02_-37_+27m_ROI\n", program);
	printf ("\toption\n");
	printf ("\t-H\tinclude header info in output\n");
	printf ("\t-V\tprint defined voxel counts per die\n");
	printf ("\t-K\tcreate die (voxel) coordinate listing\n");
	printf ("\t-Z\tcount zero voxels in <image> as defined\n");
	printf ("\t-O<int>\tselect output type (see below)\n");
	printf ("\t-f<str>\tspecify frames-to-count format (default count all frames)\n");
	printf ("\t-l<int>\tspecify length of die in voxels (default %d)\n", LENDIE);
	printf ("\t-n<int>\tspecify minimum die voxel count (default %d)\n", NCRIT);
	printf ("\t-t<flt>\tspecify svd output tolerance - ratio of least to greatest eigenvalue (default %.3g)\n", TOL);
	printf ("\t-o<str>\twrite output to specified text file (default stdout)\n");
	printf ("\toption -O output types\n");
	printf ("	1	timeseries directly extracted from dice\n");
	printf ("	2	timeseries extracted from dice with mean removed\n");
	printf ("	3	die timeseries passed through svd multiplied by eigenvalue\n");
	printf ("	4	die timeseries passed through svd (unit variance)\n");
	printf ("N.B.:\tconc files must have extension \"conc\"\n");
	printf ("N.B.:\tonly defined voxels (not 0.0 and not NaN and not 1.e-37 and finite) are counted\n");
	printf ("N.B.:\t%s ignores <(4dfp) ROI> ifh center and mmppix fields\n", program);
	printf ("N.B.:\tto obtain a GLM condition number = X specificy sqrt(1/X) as tol with option -t\n", program);
	exit (1);
}

int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	CONC_BLOCK	conc_block;			/* conc i/o control block */
	FILE            *imgfp, *ROIfp, *outfp, *keyfp;
	char            imgroot[MAXL], imgfile[MAXL], ROIroot[MAXL], ROIfile[MAXL], outfile[MAXL] = "", keyfile[MAXL], tmpfile[MAXL];

/********************/
/* iamge processing */
/********************/
	int             imgdim[4], ROIdim[4], diedim[3], vdim, imgori, ROIori, isbig, isbigm;
	int		lendie = LENDIE, ncrit = NCRIT, ndie, mdie, *mdielist;
	int		ix, iy, iz, jx, jy, jz;
	float           imgvox[3], ROIvox[3];
        float           *imgt, *imgr;
	DIEDAT		*diedat;

/*************************/
/* timeseries processing */
/*************************/
	char		*format;
	double		tol = TOL, **b, *lambda, **v;
	int		nnez, nred;

/***********/
/* utility */
/***********/
	char            *ptr, command[2*MAXL];
	double		q;
	int             c, i, j, k, l, m, n;

/*********/
/* flags */
/*********/
	int		outtype = 0;
	int		conc_flag = 0;
	int             status = 0;
	int		count_zero = 0;
	int		defined;			/* test on imgt and imgr voxels */
	int		test_svd = 1;
	int		H_flag = 0;
	int		K_flag = 0;
	int		V_flag = 0;

/**************/
/* debug code */
/**************/
	double		**bcopy, **e, errmax;

	setprog (program, argv);
	if (!(format = (char *) calloc (1024, sizeof (char)))) errm (program);
/************************/
/* process command line */
/************************/
	for (j = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'Z': count_zero++;			break; 
				case 'V': V_flag++;			break; 
				case 'H': H_flag++;			break; 
				case 'K': K_flag++;			break; 
				case 'O': outtype = atoi (ptr);		*ptr = '\0'; break; 
				case 'l': lendie = atoi (ptr);		*ptr = '\0'; break; 
				case 'f': strncpy (format, ptr, 1023);	*ptr = '\0'; break; 
				case 'n': ncrit = atoi (ptr);		*ptr = '\0'; break; 
				case 't': tol = atof (ptr);		*ptr = '\0'; break; 
				case 'o': strcpy (outfile, ptr);	*ptr = '\0'; break; 
				default:				break;
			}
		} else switch (j) {
			case 0:	getroot (argv[i], imgroot);
				conc_flag = (strstr (argv[i], ".conc") == argv[i] + strlen (imgroot));
                						j++; break;
                	case 1: getroot (argv[i], ROIroot);	j++; break;
			default:				break;
		}
	}
	if (j < 2) usage (program);
	if (lendie <= 0 || ncrit <= 0) usage (program);

	if (strlen (outfile)) {
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		printf ("Writing: %s\n", outfile);
	} else {
		outfp = stdout;
	}
	if (H_flag) write_command_line (outfp, argc, argv);

/**************************************/
/* get input 4dfp dims and open files */
/**************************************/
	if (conc_flag) {
		conc_init_quiet (&conc_block, program);
		conc_open_quiet (&conc_block, imgroot);
		strcpy (imgfile, conc_block.lstfile);
		for (k = 0; k < 4; k++) imgdim[k] = conc_block.imgdim[k];
		for (k = 0; k < 3; k++) imgvox[k] = conc_block.voxdim[k];
		imgori = conc_block.orient;
		isbig  = conc_block.isbig;
	} else {
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		if (get_4dfp_dimoe_quiet (imgfile, imgdim, imgvox, &imgori, &isbig) < 0) errr (program, imgfile);
		if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
	}

/*****************/
/* alloc buffers */
/*****************/
	vdim = imgdim[0]*imgdim[1]*imgdim[2];
	if (!(imgt = (float *) malloc (vdim*sizeof (float)))) errm (program);
	if (!(imgr = (float *) malloc (vdim*sizeof (float)))) errm (program);

/****************************/
/* prepare dice from ROIimg */
/****************************/
	sprintf (ROIfile, "%s.4dfp.img", ROIroot);
	if (get_4dfp_dimoe_quiet (ROIfile, ROIdim, ROIvox, &ROIori, &isbigm) < 0) errr (program, ROIfile);
	status = ROIori - imgori;
	for (k = 0; k < 3; k++) status |= (imgdim[k] != ROIdim[k]);
	for (k = 0; k < 3; k++) status |= (fabs (imgvox[k] - ROIvox[k]) > 0.0001);
	if (status) {
		fprintf (stderr, "%s: %s %s vdim mismatch\n", program, imgroot, ROIroot);
		exit (-1);
	}
	fprintf (stdout, "Reading: %s\n", ROIfile);
	if (!(ROIfp = fopen (ROIfile, "rb"))
	|| eread (imgr, vdim, isbigm, ROIfp)
	|| fclose (ROIfp)) errr (program, ROIfile);

	ndie = 1; for (k = 0; k < 3; k++) {
		diedim[k] = imgdim[k]/lendie;
		ndie *= diedim[k];
	}
	if (0) printf ("diedim %d %d %d ndie %d\n", diedim[0], diedim[1], diedim[2], ndie);
	if (!(diedat = (DIEDAT *) calloc (ndie, sizeof (DIEDAT)))) errm (program);
	j = 0;		/* index of die */
	for (iz = 0; iz < imgdim[2] - (imgdim[2] % lendie); iz += lendie) {
	for (iy = 0; iy < imgdim[1] - (imgdim[1] % lendie); iy += lendie) {
	for (ix = 0; ix < imgdim[0] - (imgdim[0] % lendie); ix += lendie) {
	if (0) printf ("ix, iy, iz, j %5d%5d%5d%10d\n", ix, iy, iz, j);
		diedat[j].fx = ix + lendie/2. - 0.5;
		diedat[j].fy = iy + lendie/2. - 0.5;
		diedat[j].fz = iz + lendie/2. - 0.5;
		j++;
	}}}
	assert (j == ndie);

	for (iz = 0; iz < imgdim[2] - (imgdim[2] % lendie); iz++) {jz = iz/lendie;
	for (iy = 0; iy < imgdim[1] - (imgdim[1] % lendie); iy++) {jy = iy/lendie;
	for (ix = 0; ix < imgdim[0] - (imgdim[0] % lendie); ix++) {jx = ix/lendie;
		j = jx + diedim[0]*(jy + diedim[1]*jz);		/* index of die */
		i = ix + imgdim[0]*(iy + imgdim[1]*iz);		/* index of voxel */
		defined = isnormal (imgr[i]) && imgr[i] != (float) 1.e-37;
		if (defined) {
			if (diedat[j].n == diedat[j].ma) {
				diedat[j].ma += INCR;
				if (!diedat[j].p) {
					diedat[j].p =  		    malloc (diedat[j].ma*sizeof (int));
				} else {
					diedat[j].p = realloc (diedat[j].p, diedat[j].ma*sizeof (int));
				}
				if (!diedat[j].p) errm (program);
			}
			diedat[j].p[diedat[j].n++] = i;
		}
		i++;
	}}}

/********************/
/* create mdie list */
/********************/
	for (mdie = j = 0; j < ndie; j++) {
		if (V_flag) printf ("j diedat[j].n %5d %5d\n", j, diedat[j].n);
		if (diedat[j].n >= ncrit) mdie++;
	}
	printf ("defined die count = %d\n", mdie);
	if (!(mdielist = (int *) calloc (mdie, sizeof (int)))) errm (program);
	for (mdie = j = 0; j < ndie; j++) if (diedat[j].n >= ncrit) mdielist[mdie++] = j;

	if (K_flag) {
/*********************************/
/* create voxel address key file */
/*********************************/
		sprintf (tmpfile, "temp%d", getpid ());
		if (!(keyfp = fopen (tmpfile, "w"))) errw (program, tmpfile);
		for (m = 0; m < mdie; m++) {
			j = mdielist[m];			/* index of die with at least ncrit voxels in ROI mask */
			fprintf (keyfp, "%10.2f%10.2f%10.2f\n", diedat[j].fx, diedat[j].fy, diedat[j].fz);
		}
		if (fclose (keyfp)) errw (program, tmpfile);
		strcpy (keyfile, imgroot);
		if (!(ptr = strrchr (keyfile, '/'))) ptr = keyfile; else ptr++;
		strcpy (keyfile, ptr);
		strcat (keyfile, "_");
		strcpy (command, ROIroot);
		if (!(ptr = strrchr (command, '/'))) ptr = command; else ptr++;
		strcat (keyfile, ptr);
		strcat (keyfile, "_indices.txt");
		if (!(keyfp = fopen (keyfile, "w"))) errw (program, keyfile);
		printf ("Writing: %s\n", keyfile);
		write_command_line (keyfp, argc, argv);
		if (fclose (keyfp)) errw (program, keyfile);
		sprintf (command, "index2atl %s %s >> %s",
			(conc_flag) ? conc_block.imgfile0[0] : imgroot, tmpfile, keyfile);
		printf ("%s\n", command);
		status |= system (command);
		remove (tmpfile);
		if (status) exit (-1);
	}

/****************/
/* parse format */
/****************/
	if (!(format = (char *) realloc (format, (imgdim[3] + 1)*sizeof (char)))) errm (program);
	if (strlen (format)) {
		if (k = expandf (format, imgdim[3] + 1)) exit (-1);
		if (strlen (format) != imgdim[3]) {
			fprintf (stderr, "%s: %s frame-to-count format length mismatch\n", program, imgroot);
			exit (-1); 
		}
		for (nnez = k = 0; k < imgdim[3]; k++) if (format[k] != 'x') nnez++;
	} else {
		nnez = imgdim[3];
		for (i = 0; i < imgdim[3]; i++) format[i] = '+';
	}
	printf ("%s\n", format);
	printf ("%s: time series defined for %d frames, %d exluded\n", program, imgdim[3], imgdim[3] - nnez);

/***************************/
/* allocate memory for svd */
/***************************/
	if (!(lambda = (double *) calloc (nnez, sizeof (double)))) errm (program);
	v = calloc_double2 (nnez, nnez);
	b = calloc_double2 (mdie, nnez);

/***********/
/* process */
/***********/
	printf ("Reading: %s\n", imgfile);
	for (n = k = 0; k < imgdim[3]; k++) {
		if (conc_flag) {
			conc_read_vol (&conc_block, imgt);
		} else {
			if (eread (imgt, vdim, isbig, imgfp)) errr (program, imgfile);
		}
		for (m = 0; m < mdie; m++) {
			j = mdielist[m];			/* index of die with at least ncrit voxels in ROI mask */
			diedat[j].v = diedat[j].m = 0;		/* initialize anew for each frame */
			for (l = 0; l < diedat[j].n; l++) {	/* loop on voxels */
				i = diedat[j].p[l];		/* index into imgt */
				defined = isnormal (imgt[i]) && imgt[i] != (float) 1.e-37;
				if (count_zero && imgt[i] == 0.) defined = 1;
				if (defined) {
					diedat[j].v += imgt[i];
					diedat[j].m++;
				}
			}
			if (diedat[j].m < ncrit) {
				fprintf (stderr, "%s: %s contains undefined voxels within %s mask\n",
								program, imgfile, ROIroot);
				exit (-1);
			}
			diedat[j].v /= diedat[j].m;
		}
		if (format[k] != 'x') {
			for (m = 0; m < mdie; m++) {
				j = mdielist[m];		/* index of die with at least ncrit voxels in ROI mask */
				b[m][n] = diedat[j].v;
			}
			n++;					/* increment non-skipped frames count */
		}
	}
	assert (n == nnez);
	free (imgt); free (imgr);
	if (conc_flag) {
		conc_free (&conc_block);
	} else {
		if (fclose (imgfp)) errr (program, imgfile);
	}
	if (!outtype) goto DONE;

/*****************************************************/
/* make timeseries zero mean (but not unit variance) */
/*****************************************************/
	if (outtype > 1) for (j = 0; j < mdie; j++) dzeromean (b[j], nnez);

/*******************************/
/* output extracted timeseries */
/*******************************/
	if (outtype == 1 || outtype == 2) {
		for (n = k = 0; k < imgdim[3]; k++) {
			for (m = 0; m < mdie; m++) fprintf (outfp, " %9.4f", (format[k] != 'x') ? b[m][n] : 0.0);
			if (format[k] != 'x') n++;
			if (mdie) fprintf (outfp, "\n");
		}
		assert (n == nnez);
		goto DONE;
	}

	if (test_svd) {
		bcopy = calloc_double2 (mdie, nnez);
   	   	    e = calloc_double2 (mdie, nnez);
		for (j = 0; j < mdie; j++) for (i = 0; i < nnez; i++) bcopy[j][i] = b[j][i];
	}
   	        dsvdcmp0 (b, mdie, nnez, lambda, v);
	nred = ndsvdcmp0 (b, mdie, nnez, lambda, v, tol);
	if (H_flag) {
		fprintf (outfp, "#reduced dimensionality = %d\n", nred);
		fprintf (outfp, "#lambda after ndsvdcmp0");
		for (i = 0; i < nred; i++) fprintf (outfp, " %.2f", lambda[i]); fprintf (outfp, "\n");
	}
	if (test_svd) {
		for (j = 0; j < mdie; j++) for (i = 0; i < nnez; i++) for (k = 0; k < nred; k++) {
			e[j][i] += v[i][k]*lambda[k]*b[j][k];
		}
		for (errmax = i = 0; i < nnez; i++) for (j = 0; j < mdie; j++) {
			q = fabs (e[j][i] - bcopy[j][i]);
			if (q > errmax) errmax = q;
		}
		if (H_flag) fprintf (outfp, "#maximum error = %.3g\n", errmax);
		free_double2 (bcopy); free_double2 (e);
	}
	for (k = 0; k < nnez; k++) for (j = 0; j < nred; j++) {
		v[k][j] *= (outtype == 3) ? lambda[j] : sqrt ((double) nnez);
	}
/******************************/
/* output computed timeseries */
/******************************/
	if (outtype == 3 || outtype == 4) {
		for (n = k = 0; k < imgdim[3]; k++) {
			for (m = 0; m < nred; m++) fprintf (outfp, " %9.4f", (format[k] != 'x') ? v[n][m] : 0.0);
			if (format[k] != 'x') n++;
			if (mdie) fprintf (outfp, "\n");
		}
		assert (n == nnez);
	}

/*********************/
/* clean up and exit */
/*********************/
DONE:	if (strlen (outfile)) if (fclose (outfp)) errw (program, outfile);
	free_double2 (v); free_double2 (b); free (lambda);
	free (mdielist);
	for (j = 0; j < ndie; j++) if (diedat[j].p) free (diedat[j].p);
	free (diedat);
	free (format);
	exit (status);
}
