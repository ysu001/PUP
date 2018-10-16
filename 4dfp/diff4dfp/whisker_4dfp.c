/*$Header: /home/usr/shimonyj/diff4dfp/RCS/whisker_4dfp.c,v 1.8 2008/09/16 02:54:57 avi Exp $*/
/*$Log: whisker_4dfp.c,v $
 * Revision 1.8  2008/09/16  02:54:57  avi
 * reorder thresh_flg and bvox_flg logic
 *
 * Revision 1.7  2008/09/16  02:26:31  avi
 * linux and little endian compliant
 * fimg_mode() -> fimg_modetn
 *
 * Revision 1.6  2005/12/08  03:21:46  avi
 * install calls to t4tolin() and affine_DTrot() to deal with affine transformed DWI data
 * prevent main voxel loop from skipping output (remove continue statements)
 * fix imgt address bug (i -> jndex)
 * properly free memory at main program termination
 *
 * Revision 1.5  2005/11/27  23:02:55  avi
 * comment typos and update utility subroutines
 *
 * Revision 1.4  2005/08/06  01:08:40  avi
 * -3 option (for creating diffusion ellipse displays)
 * eliminate all references to FORTRAN (flip?_)
 *
 * Revision 1.3  2002/04/30  23:39:42  avi
 * -E (output eignevalues) option
 *
 * Revision 1.2  2001/06/13  19:50:28  avi
 * Alpay Matlab => sign invert eigenvector z component
 *
 * Revision 1.1  2001/06/11  21:16:27  avi
 * Initial revision
 **/
/*****************************************/
/* diffusion calculations on 4dfp stacks */
/*****************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <JSSutil.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>

/***************/
/* fimg_mode.c */
/***************/
float	fimg_modetn (float* fimg, int nval, float histmin, int nsmooth);
float	fimg_mode (float*, int);

/***********/
/* cflip.c */
/***********/
extern void	flipx (float *imgf, int *pnx, int* pny, int *pnz);
extern void	flipy (float *imgf, int *pnx, int* pny, int *pnz);
extern void	flipz (float *imgf, int *pnx, int* pny, int *pnz);

/********************/
/* get_dti_params.c */
/********************/
int	dti_dimen (char *file);
int	get_dti_params_nrc (char *file, int n, float *b_vals, float **q_vals);

/*************/
/* t4tolin.c */
/*************/
void	t4tolin (char *t4file, int orient, float **lin);

/******************/
/* affine_DTrot.c */
/******************/
void	affine_DTrot (float **atm, float **evec, float **cevec);

/*********************************/
/* size of input and output sets */
/*********************************/
#define MAXL		256
#define I0_THRESH	0.1		/* threshold as franction of I0 image mode value */
#define HISTMIN		100.		/* minimum value for fimg_modetn() to find I0 image mode value */
#define NSMOOTH		4		/* fimg_modetn() smoothing parameter */
#define DELX		1
#define DELY		1
#define DELZ		1

static char rcsid[] = "$Id: whisker_4dfp.c,v 1.8 2008/09/16 02:54:57 avi Exp $";
int main (int argc, char *argv[]) {
/**************/
/* processing */
/**************/
	char	prmfile[MAXL] = "/home/usr/shimonyj/diff4dfp/tp7_params.dat";
	int 	n, bvox_flg, bvox_num;
	double	Asig, prolaticity, Dbar;
	float	I0_mode = 500., I0_thresh = I0_THRESH, histmin = HISTMIN;
	int	nsmooth = NSMOOTH;

/***************/
/* 4dfp images */
/***************/
	FILE	*fp;
	char	imgroot[MAXL], outroot[MAXL];
	char	outfile[MAXL], ifhfile[MAXL], imgfile[MAXL], t4file[MAXL] = "";

	float	*imgt;
	float 	voxdim[3];
	int	orient, imgdim[4], vdim, slcdim, isbig;

/*************************************/
/* memory for diffusion calculations */
/*************************************/
	float	*ival, *bval, *ang, *dvec, *eigenval;
	float	**qval, **eigenvec, **reigenvec, **dinv, **exp_arr, **lin;

/***********/
/* utility */
/***********/
	int	c, i, j, k, l, jndex;
	int	ix, iy, iz, delx = DELX, dely = DELY, delz = DELZ;
	char	*str,  program[MAXL], command[MAXL];

/*********/
/* flags */
/*********/
	int	thresh_flg;
	int	Eigvals = 0;
	int	axes_flg = 0;
	int	transform = 0;

	printf ("%s\n", rcsid);
	if (!(str = strrchr (argv[0], '/'))) str = argv[0]; else str++;
	strcpy (program, str);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') { 
			strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case '3': axes_flg++;			break;
				case 'E': Eigvals++;			break;
				case 'n': nsmooth = atoi (str);		*str = '\0'; break;
				case 'h': histmin = atof (str);		*str = '\0'; break;
				case 't': I0_thresh = atof (str);	*str = '\0'; break;
				case 'T': strcpy (t4file, str);		*str = '\0'; break;
				case 'd': j = *str++; switch (j) {
					case 'x': case 'X': delx = atoi (str);	*str = '\0'; break;
					case 'y': case 'Y': dely = atoi (str);	*str = '\0'; break;
					case 'z': case 'Z': delz = atoi (str);	*str = '\0'; break;
				} break;
			}
		} else switch (k) {
			case 0: strcpy  (prmfile, argv[i]); k++; break;
			case 1: getroot (argv[i], imgroot); k++; break;
		}
	}
	if (k != 2) {
		printf ("Usage:\t%s <prm_file> <file_4dfp>\n", program);
		printf (" e.g.,\t%s tp7_params.dat -dz3 /data/emotion/data3/track_sub3/track_sub3_DTI_avg\n", program);
		printf ("	option\n");
		printf ("	-h<flt>\tspecify minimum I0 mode (default=%.2f)\n", HISTMIN);
		printf ("	-n<int>\tspecify number of I0 histogram smoothings (default=%d)\n", NSMOOTH);
		printf ("	-t<flt>\tspecify mask threshold as fraction of I0 mode (default=%.4f)\n", I0_THRESH);
		printf ("	-T<str>\tspecify t4 file used to transform DWI data\n");
		printf ("	-E\tadditionally output eigenvalues\n");
		printf ("	-3\toutput 3 eigenvectors scaled by eigenvalue\n");
		printf ("	-d<x|y|z><int>\tspecify quiver spacing in pixels (default=1)\n");
		printf ("N.B.:	default output is first two eigenvectors scaled by Asigma\n");
		exit (1);
	}

	if (delx < 1) delx = 1;
	if (dely < 1) dely = 1;
	if (delz < 1) delz = 1;
	printf ("delx=%d dely=%d delz=%d\n", delx, dely, delz);

/******************************/
/* get number of measurements */
/******************************/
	n = dti_dimen (prmfile);

/**************************************/
/* allocate memory for diffusion calc */
/**************************************/
	ival =		vector(1, n);
	bval =		vector(1, n);
	eigenval =	vector(1, 3);
	ang =		vector(1, 9);
	dvec =		vector(1, 7);
	if (!ival || !bval || !eigenval || !ang || !dvec) errm (program);

	exp_arr = 	matrix(1, n, 1, 7);
	eigenvec =	matrix(1, 3, 1, 3);
	reigenvec =	matrix(1, 3, 1, 3);	/* rotated eigenvectors */
	lin =		matrix(1, 3, 1, 3);	/* linear transform derived from t4 file */
	dinv =		matrix(1, 3, 1, 3);
	qval =		matrix(1, n, 1, 3);
	if (!exp_arr || !eigenvec || !reigenvec || !lin || !dinv || !qval) errm (program);

/*****************************************/
/* get the b values and the unit vectors */
/*****************************************/
	get_dti_params_nrc (prmfile, n, bval, qval);

/********************************************/
/* set up the svd for diffusion calculation */
/********************************************/
	set_svd (n, bval, qval, exp_arr);

/********************************************/
/* get stack dimensions and optional t4file */
/********************************************/
	if (get_4dfp_dimoe (imgroot, imgdim, voxdim, &orient, &isbig)) errr (program, imgroot);
	if (imgdim[3] != n) {
		fprintf (stderr, "%s: %s %s dimension mismatch\n", program, prmfile, imgroot);
		exit (-1);
	}
	if (strlen (t4file)) t4tolin (t4file, orient, lin);

/***************************************/
/* allocate input memory and get stack */
/***************************************/
	vdim =		imgdim[0] * imgdim[1] * imgdim[2];
	slcdim =	imgdim[0] * imgdim[1];
	if (!(imgt = (float *) malloc (n*vdim*sizeof(float)))) errm (program);
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	printf ("Reading: %s\n", imgfile);
	if (!(fp = fopen (imgfile, "rb")) || eread (imgt, n*vdim, isbig, fp)
	|| fclose (fp)) errr (program, imgfile);

/************************/
/* flip 4dfp to analyze */
/************************/
	for (l = 0; l < n; l++) switch (orient) {
		case 4:	flipx (imgt + l*vdim, imgdim+0, imgdim+1, imgdim+2);
		case 3:	flipz (imgt + l*vdim, imgdim+0, imgdim+1, imgdim+2);
		case 2:	flipy (imgt + l*vdim, imgdim+0, imgdim+1, imgdim+2); break;
		default: fprintf (stderr, "%s: illegal %s orientation (%d)\n", program, imgfile, orient);
		exit (-1);
	}

/******************/
/* compute I0 mode*/
/******************/
	I0_mode = fimg_modetn (imgt, vdim, histmin, nsmooth);

/********************************************/
/* loop over the images and calc dti params */
/********************************************/
	sprintf (outfile, "%s_whisker.dat", imgroot);
	printf ("Writing: %s\n", outfile);
	if (!(fp = fopen (outfile, "w"))) errw (program, outfile);

	bvox_num = 0;
	printf ("processing slice:");
	for (iz = 0; iz < imgdim[2]; iz += delz) {printf (" %d", iz + 1); fflush (stdout);
	for (iy = 0; iy < imgdim[1]; iy += dely) {
	for (ix = 0; ix < imgdim[0]; ix += delx) {
		jndex = ix + imgdim[0] * (iy + imgdim[1] * iz);
		thresh_flg = (imgt[jndex] < I0_thresh * I0_mode);
		for (bvox_flg = j = 0; j < n; j++) {
			if ((ival[j+1] = (imgt + j*vdim)[jndex]) <= 0.0) bvox_flg++;
		}
		if (!bvox_flg && !thresh_flg) {
			calc_svd (n, ival, exp_arr, dvec);
			eigen_calc (dvec, dinv, eigenval, eigenvec, ang);
			Dbar = ang[4]; if (Dbar < 0.0) {bvox_flg++; bvox_num++;}
			Asig = ang[5]; if (Asig < 0.0) Asig = 0.0; if (Asig > 1.0) Asig = 1.0;

/***************************************/
/* accommodate affine transformed data */
/***************************************/
			if (strlen (t4file)) {
				affine_DTrot (lin, eigenvec, reigenvec);			
			} else {
				for (k = 1; k <= 3; k++) for (i = 1; i <= 3; i++) reigenvec[k][i] = eigenvec[k][i];
			}

/*******************************************************************/
/* sign invert z component based on Alpay's matlab display results */
/*******************************************************************/
			for (i = 1; i <= 3; i++) reigenvec[3][i] *= -1.0;
		}
		if (thresh_flg || bvox_flg) {
			for (i = 1; i <= 3; i++) eigenval[i] = 0.0;
			Asig = 0.0;
		}

		fprintf (fp, "%10d%10d%10d", ix + 1 , iy + 1, iz + 1);
		if (axes_flg) {
			for (i = 3; i > 0; i--) {
				for (k = 1; k <= 3; k++) {
					fprintf (fp, "%10.4f", eigenval[i]*reigenvec[k][i]);
				}
			}
		} else {
				for (k = 0; k < 3; k++) fprintf (fp, "%10.4f", Asig*reigenvec[k + 1][3]);
				for (k = 0; k < 3; k++) fprintf (fp, "%10.4f", Asig*reigenvec[k + 1][2]);
		}
		if (Eigvals)	for (k = 0; k < 3; k++) fprintf (fp, "%10.5f", eigenval[k + 1]);
		fprintf (fp, "\n");
	}}}
	printf ("\nbad voxel count %d out of %d\n", bvox_num, vdim);
	fclose (fp);

/*******************/
/* create rec file */
/*******************/
	startrec (outfile, argc, argv, rcsid);
	printrec ("dti parameter file\n");
	catrec (prmfile);
	if (strlen (t4file)) {
		printrec ("affine transform file\n");
		catrec (t4file);
	}
	sprintf (command, "delx=%d dely=%d delz=%d\n", delx, dely, delz); printrec (command);
	sprintf (command, "Bad voxel fraction %d/%d\n", bvox_num, vdim); printrec (command);
	sprintf (command, "Mask threshold = %.4f*I0 mode = %.2f\n", I0_thresh, I0_thresh * I0_mode); printrec (command);
	catrec (imgfile);
	endrec ();

/***************/
/* free memory */
/***************/
	free (imgt);
	free_vector (ival, 1, n);
	free_vector (bval, 1, n);
	free_vector (eigenval, 1, 3);
	free_vector (ang, 1, 9);
	free_vector (dvec, 1, 7);
	free_matrix (exp_arr, 1, n, 1, 7);
	free_matrix (eigenvec, 1, 3, 1, 3);
	free_matrix (reigenvec, 1, 3, 1, 3);
	free_matrix (lin, 1, 3, 1, 3);
	free_matrix (dinv, 1, 3, 1, 3);
	free_matrix (qval, 1, n, 1, 3);
	return 0;
}
