/*$Header: /home/usr/shimonyj/diff4dfp/RCS/diff_4dfp.c,v 1.32 2010/03/24 03:08:33 avi Exp $*/
/*$Log: diff_4dfp.c,v $
 * Revision 1.32  2010/03/24  03:08:33  avi
 * eliiminate complexities is computation of isbig
 *
 * Revision 1.31  2010/03/24  02:49:38  avi
 * Usage and minor bug fixes
 *
 * Revision 1.30  2010/01/25  05:40:00  shimonyj
 * options -a, -M, -Z, -B, -b, -R, -r, -s
 * output volumes reordered - Dtensor follows Asigma
 *
 * Revision 1.29  2008/09/16  03:09:47  avi
 * linux compliant
 * #include <JSSutil.h>
 *
 * Revision 1.28  2008/06/11  03:26:48  avi
 * install corrected Levenberg-Marquardt code
 * remove option -v (xyz encoded DWI)
 *
 * Revision 1.27  2008/05/26  03:48:05  avi
 * remove Bayes and xyz options
 *
 * Revision 1.26  2008/04/14  01:35:36  adrian
 * recode FA logic
 *
 * Revision 1.25  2008/03/13  02:05:24  avi
 * option -F (output FA)
 *
 * Revision 1.24  2007/10/15  21:18:56  avi
 * eliminate doubly-defined utilities in endianio.c
 *
 * Revision 1.23  2007/10/15  21:14:11  avi
 * endian-compliant i/o
 * LM code
 *
 * Revision 1.22  2006/03/18  07:30:27  avi
 * remove additive constant option -A
 * install better I0 threshold routine (fimg_mode() -> fimg_modetn())
 * get input 4dfp orient using get_4dfp_dimo() instead of Getifh()
 *
 * Revision 1.21  2005/08/10  02:43:23  avi
 * updated Acom computation
 * optionally output up to three eigenvectors
 *
 * Revision 1.20  2005/02/28  22:57:07  avi
 * Bayes working
 *
 * Revision 1.19  2004/12/09  23:39:15  avi
 * Bayes calls inserted but not debugged
 *
 * Revision 1.18  2004/12/07  22:11:02  avi
 * do not skip tensor computation if an input image has negative intensity
 * force 0 <= Asig <= 1.
 *
 * Revision 1.17  2002/04/10  03:44:03  avi
 * check for existence of dti parameter file
 *
 * Revision 1.16  2001/10/12  17:54:40  avi
 * also disable -E in "vector only" processing
 *
 * Revision 1.15  2001/09/26  04:24:17  avi
 * zero Asig values outside range (0, 3) rather than (0, 1)
 *
 * Revision 1.14  2001/07/24  02:36:56  avi
 * -E option (output eigenvalues)
 *
 * Revision 1.13  2001/07/04  20:35:03  avi
 * xyz_flag
 * new definition of prolaticity
 *
 * Revision 1.12  2001/06/23  02:54:37  avi
 * ensure normal exit status = 0
 *
 * Revision 1.11  2001/05/31  23:57:26  avi
 * correct transposed readout of eigenvectors
 *
 * Revision 1.10  2001/04/01  23:44:11  avi
 * options -V (principal eigenvector) and -P (output prolaticity)
 *
 * Revision 1.9  2001/01/13  22:14:09  avi
 * -D (D tensor output) option
 *
 * Revision 1.8  2001/01/07  05:02:08  avi
 * include input rec file in output recfile
 *
 * Revision 1.7  2001/01/07  03:26:48  avi
 * remove all FORTRAN references
 * read orientation from ifh file
 * use get_4d_dimensions () instead of Get4dfpDimN ()
 *
 * Revision 1.6  2000/12/14  22:33:04  shimonyj
 * prepare file for shipment to nih, eliminate header file
 *
 * Revision 1.5  2000/10/25  23:49:00  avi
 * generalize image input to accept multi-contrast stacks
 *
 * Revision 1.4  2000/10/20  07:21:02  avi
 * output prolaticity image
 * value limits put on output computed parameters
 * masking by thresholding I0 relative to the I0 image mode value
 *
 * Revision 1.3  2000/10/11  05:20:30  avi
 * 4dfp conventions and overhead
 *
 * Revision 1.2  2000/10/06  22:25:52  shimonyj
 * add interp_flg
 *
 * Revision 1.1  2000/10/06  22:14:04  shimonyj
 * Initial revision
 **/
/*****************************************/
/* diffusion calculations on 4dfp stacks */
/*****************************************/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>	/* R_OK */
#include <math.h>
#include <string.h>
#include <JSSutil.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>
#include <knee.h>

/***************/
/* fimg_mode.c */
/***************/
float	fimg_modetn (float* fimg, int nval, float histmin, int nsmooth);
float	fimg_mode (float*, int);

/********************/
/* get_dti_params.c */
/********************/
int	dti_dimen (char *file);
int	get_dti_params_nrc (char *file, int n, float *b_vals, float **q_vals);

/***************/
/* nonlinear.c */
/***************/
int nonlin_dti (int ndir, float *in, float *bval, float **bvec, int badflg, int CO_flag, float *eval, float *ang,
		float *parval, float *parstd);

#define MAXL		256
#define NOUT		2
#define I0_THRESH	0.1		/* threshold as fraction of I0 image mode value */
#define HISTMIN		100.		/* minimum value for fimg_modetn() to find I0 image mode value */
#define NSMOOTH		4		/* fimg_modetn() smoothing parameter */
#define MA		8		/* number of coefficients fitted by LM algprithm */

/**************************************/
/* bad encode subroutine in this file */
/**************************************/
void fit(float x[], float y[], int ndata, float sig[], int mwt,
	float *a, float *b, float *siga, float *sigb, float *chi2, float *q, float *lcc);
void eval2dmat(float eval1, float eval2, float eval3, float ang1, float ang2, float ang3, float **dmat);

/********************************************/
/* global memory for diffusion calculations */
/********************************************/
static char rcsid[] = "$Id: diff_4dfp.c,v 1.32 2010/03/24 03:08:33 avi Exp $";
int 	nmeas;
float	*ival, *bval, *ang, *dvec, *eigenval;
float	**qval, **eigenvec, **dinv, **exp_arr;

void usage (char *program) {
	printf ("Usage:\t%s <prm_file> <file_4dfp>\n", program);
	printf (" e.g.,\t%s tp7_params.dat /data/emotion/data3/track_sub3/track_sub3_DTI_avg\n", program);
	printf ("	\nCOMPUTATIONAL OPTIONS:\n");
	printf ("	-N\tcompute D using nonlinear Levenberg-Marquardt (default log linear LS)\n");
	printf ("	-Z\tuse nonlinear approach to repair bad voxels from log linear LS\n");
	printf ("	-c\testimate non-mobile diffusion term (applies only to Levenberg-Marquardt algorithm)\n");
	printf ("	-s<int>\tcompute only selected slice for debugging\n");
	printf ("	-v<int>\tcompute only selected voxel for debugging\n");
	printf ("	-B<flt>\tignore bad encoding at specified threshold (units = s.d.) (default=2.0)\n");
	printf ("	-b<flt>\textend -B functionality to background voxels\n");
	printf ("	\nMASKING OPTIONS:\n");
	printf ("	-M\tcompute threshold mask without holes\n");
	printf ("	-t<flt>\tspecify mask threshold as fraction of I0 mode (default=%.4f)\n", I0_THRESH);
	printf ("	-h<flt>\tspecify minimum I0 mode (default=%.2f)\n", HISTMIN);
	printf ("	-n<int>\tspecify number of I0 histogram smoothings (default=%d)\n", NSMOOTH);
	printf ("	\nOUTPUT OPTIONS:\n");
	printf ("	-a<str>\tappend specified trailer to output fileroots\n");
	printf ("	-p\tprint out pixel numbers for debugging\n");
	printf ("	-D\toutput D tensor\n");
	printf ("	-F\toutput FA (fractional anisotropy)\n");
	printf ("	-E\toutput eigenvalues\n");
	printf ("	-V[int]\toutput [specified number of (default=1)] eigenvectors (principal first)\n");
	printf ("	-P\toutput prolaticity\n");
	printf ("	-R\toutput single residue volume for model\n");
	printf ("	-r\toutput squared residue values for all encodings in a separate file\n");
	printf ("	-o\toutput extra full LM output files (applies only to LM algorithm)\n");
	printf ("	-d\tdebug mode, provide extra volume of output as needed \n");
	printf ("	-@<b|l>\toutput big or little endian (default input endian)\n");
	printf ("\nN.B.:	the first data volume must have high SNR from b=0 or low b value \n");
	printf ("N.B.:	optional output volumes are appended to MD and RA\n");
	printf ("N.B.:	output order: MD,RA,(Dxx,Dyy,Dzz,Dxy,Dxz,Dyz),(FA),(E123,RD),(CO),(Res),(Evecs),(Prol)\n");
	printf ("N.B.:	-b implies -B\n");
	printf ("N.B.:	-o implies -N\n");
	printf ("N.B.:	-B parameter useful range is 1.5 to 3\n");
	printf ("N.B.:	eigenvalue ordering is = Eval3 > Eval2 > Eval1\n");
	exit (1);
}

int main (int argc, char *argv[]) {
/**************/
/* processing */
/**************/
	char	prmfile[MAXL] = "/home/usr/shimonyj/diff4dfp/tp7_params.dat";
	int 	bvox_flg, bvox_num;
	double	Asig, FA, prolaticity, Dbar, Asig0 = 0.0;
	float	I0_max, I0_mode = 500., I0_thresh = I0_THRESH, histmin = HISTMIN;
	int	nsmooth = NSMOOTH;	/* I0 voxel value histogram */
	int	slicesel, voxsel; 	/* slice and voxel for special testing */

/***************/
/* 4dfp images */
/***************/
	FILE	*fp;
	IFH	ifh;
	char	imgroot[MAXL], outroot[MAXL], trailer[MAXL] = "";
	char	outfile[MAXL], ifhfile[MAXL], imgfile[MAXL];
	char	*outtype0[] = {"Dbar", "Asigma"};
	char	*outtypeF[] = {"Fractional Anisotropy"};
	char	*outtypeC[] = {"Additive Constant"};
	char	*outtypeR[] = {"Sqrt of sum of residuals squared"};
	char	*outtypeP[] = {"prolaticity"};
	char	*outtypeV[] = {"Lambda1 eigenvec.x", "Lambda1 eigenvec.y", "Lambda1 eigenvec.z",
			       "Lambda2 eigenvec.x", "Lambda2 eigenvec.y", "Lambda2 eigenvec.z",
			       "Lambda3 eigenvec.x", "Lambda3 eigenvec.y", "Lambda3 eigenvec.z"};
	char	*outtypeD[] = {"Dxx", "Dyy", "Dzz", "Dxy", "Dxz", "Dyz"};
	char	*outtypeE[] = {"Lambda1", "Lambda2", "Lambda3 = Axial Diff", "Radial Diff"};
	char	*labels[32];

	float	*in_img, *out_img;
	float 	voxdim[3];
	int	imgdim[4], vdim, slcdim, nout = NOUT, orient, isbig;
	char	control = '\0';

/**************************/
/* nonlinear LM algorithm */
/**************************/
	FILE	*fp_info;
	char	Info_file[MAXL];
	char 	*outtypeN[] = { "num Converge fits", "Lambda1", "Lambda2", "Lambda3", "Phi", "Theta", "Psi",
				"Amplitude", "Constant", "sqrt sum residuals sq",
				"stddev Lambda1", "stddev Lambda2", "stddev Lambda3", 
				"stddev Phi", "stddev Theta", "stddev Psi",
				"stddev Amplitude", "stddev Constant"};

	char	LMoutroot[MAXL], LMoutfile[MAXL];
	float	*LM_out_img;
	int 	indx, indx1, indx2, indx3, LMout = (MA+1)*2; /* 0-returned Q, 1:8-parval, 9:calc resid, 10:17-parstd */ 
	double	sa,ca,sb,cb,sg,cg,ex,ey,ez,res,res2,drot,ival2;
	float	parval[MA + 1], parstd[MA + 1], **dmat, **dmatd;

/***********************/
/* full residue output */
/***********************/
	float	*res_out_img;	
	char	RESoutroot[MAXL], RESoutfile[MAXL];

/*****************/
/* bad encodings */
/*****************/
	int	**badenc, *badvox, nbadenc, nbadtot, nslc, nnmeas, nmeas0, nbad;
	int	ibad1, ibad2, ibad3;
	float	*bval0, **qval0, *resid;
	float	**avgenc, **numenc, **abgenc, **nbgenc, *xvals, *yvals, *sigs;
	float	bstd_thresh = 2.0;
	float	afit, bfit, siga, sigb, chi2, qq, lcc, resid2;
	float   sres, s2res;
	short	*img0, *img1, *imsk;

/***********/
/* utility */
/***********/
	int	c, i, j, k, l, m, n, ix, iy, iz, jx, jy, jz, inew;
	int	nblob, amax, imax;
	char	*str,  program[MAXL], command[MAXL];
	float 	wavg, wgt;

/*********/
/* flags */
/*********/
	int	status = 0;
	int	Prolat = 0;
	int	Eigvecs = 0;
	int	Dtensor = 0;
	int	Eigvals = 0;
	int	FA_flag = 0;
	int	LM_flag = 0;
	int	LMout_flag = 0;
	int 	CO_flag = 0;
	int 	Z_flag = 0;
	int 	M_flag = 0;
	int 	B_flag = 0;
	int 	b_flag = 0;
	int 	s_flag = 0;
	int 	pvox_flag = 0;
	int 	debug = 0;
	int 	vox_flag = 0;
	int	Res_flag = 0;
	int	res_flag = 0;

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
				case 'd': debug++;  		nout +=1;	break;
				case 'N': LM_flag++;				break;
				case 'o': LMout_flag++;	LM_flag++;		break;
				case 'c': CO_flag++; LM_flag++;	nout += 1;	break;
				case 'P': Prolat++;		nout += 1;	break;
				case 'R': Res_flag++;		nout += 1;	break;
				case 'D': Dtensor++;		nout += 6;	break;
				case 'E': Eigvals++;		nout += 4;	break;
				case 'F': FA_flag++;		nout += 1;	break;
				case 'r': res_flag++;				break;
				case 'Z': Z_flag++; 				break;
				case 'M': M_flag++; 				break;
				case 'B': B_flag++; bstd_thresh = atof(str); 	break;
				case 'V': l = atoi (str); if (l < 0 || l > 3) usage (program); if (l == 0) l = 1;
					  Eigvecs = l;  nout += 3*l;		*str = '\0'; break;
				case 't': I0_thresh = atof (str);		*str = '\0'; break;
				case 'a': sprintf (trailer, "_%s", str);	*str = '\0'; break;
				case 'n': nsmooth = atoi (str);			*str = '\0'; break;
				case 's': s_flag++; slicesel = atoi (str);	*str = '\0'; break;
				case 'v': vox_flag++; voxsel = atoi (str); 	*str = '\0'; break;
				case 'p': pvox_flag++; 				*str = '\0'; break;
				case 'h': histmin = atof (str);			*str = '\0'; break;
				case '@': control = *str++;			*str = '\0'; break;
				case 'b': B_flag++; b_flag++; bstd_thresh = atof(str); 	break;
			}
		} else switch (k) {
			case 0: strcpy (prmfile, argv[i]);	k++; break;
			case 1: getroot (argv[i], imgroot);	k++; break;
		}
	}
	if (k != 2) usage (program);

/******************************/
/* get number of measurements */
/******************************/
	if (access (prmfile, R_OK) || (nmeas0 = dti_dimen (prmfile)) < 1) errr (program, prmfile);
	if ((LM_flag || B_flag || Z_flag) && nmeas0 < 12) {
		fprintf (stderr, "%s: Estimation with flags: N,B can't be performed when number of encodings (%d) < 12\n", program, nmeas0);
		usage (program);
	}
	nmeas = nmeas0;

/**************************************/
/* allocate memory for diffusion calc */
/**************************************/
	ival =		vector(1, nmeas0);
	bval =		vector(1, nmeas0);
	eigenval =	vector(1, 3);
	ang =		vector(1, 9);
	dvec =		vector(1, 7);
	if (!ival || !bval || !eigenval || !ang || !dvec) errm (program);
	exp_arr = 	matrix(1, nmeas0, 1, 7);
	eigenvec =	matrix(1, 3, 1, 3);
	dinv =		matrix(1, 3, 1, 3);
	qval =		matrix(1, nmeas0, 1, 3);
	if (!exp_arr || !eigenvec || !dinv || !qval) errm (program);

	dmat = 		matrix(1, 3, 1, 3);
	dmatd = 	matrix(1, 3, 1, 3);

/*****************************************/
/* get the b values and the unit vectors */
/*****************************************/
	get_dti_params_nrc (prmfile, nmeas0, bval, qval);

/********************************************/
/* set up the svd for diffusion calculation */
/********************************************/
	set_svd (nmeas0, bval, qval, exp_arr);
	
/************************/
/* get stack dimensions */
/************************/
	if (Getifh (imgroot, &ifh)) errr (program, imgroot);
	isbig = strcmp (ifh.imagedata_byte_order, "littleendian");
	printf ("isbig=%d\n", isbig);
	if (!control) control = (isbig) ? 'b' : 'l';

	for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];
	if (imgdim[3] != nmeas0) {
		fprintf (stderr, "%s: %s %s dimension mismatch\n", program, prmfile, imgroot);
		exit (-1);
	}

/***************************************/
/* allocate input memory and get stack */
/***************************************/
	vdim =		imgdim[0] * imgdim[1] * imgdim[2];
	slcdim =	imgdim[0] * imgdim[1];
	if (!(in_img = (float *) malloc (nmeas0*vdim*sizeof(float)))) errm (program);
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	printf ("Reading: %s\n", imgfile);
	if (!(fp = fopen (imgfile, "rb")) || eread (in_img, nmeas0*vdim, isbig, fp)
	|| fclose (fp)) errr (program, imgfile);

/*******************/
/* compute I0 mode */
/*******************/
	I0_mode = fimg_modetn (in_img, vdim, histmin, nsmooth);
	printf ("%s computed I0 mode %.2f\n", imgfile, I0_mode);
	for (I0_max = i = 0; i < vdim; i++) if (in_img[i] > I0_max) I0_max = in_img[i];

/*****************************************/
/* allocate bad encoding specific memory */
/*****************************************/
	if (B_flag) {
		avgenc = 	matrix(0, imgdim[2], 0, nmeas0);
		numenc = 	matrix(0, imgdim[2], 0, nmeas0);
		abgenc = 	matrix(0, imgdim[2], 0, nmeas0);
		nbgenc = 	matrix(0, imgdim[2], 0, nmeas0);
		xvals =		vector(1, nmeas0);
		yvals =		vector(1, nmeas0);
		sigs  =		vector(1, nmeas0);
		bval0 =		vector(1, nmeas0);
		qval0 =		matrix(1, nmeas0, 1, 3);
	}
	badenc = 	imatrix(0, imgdim[2], 0, nmeas0);
	badvox = ivector(0, vdim);
	resid = vector(0, vdim);
	img0 = (short int *) malloc(slcdim*sizeof(short int));
	img1 = (short int *) malloc(slcdim*sizeof(short int));
	imsk = (short int *) malloc(vdim*sizeof(short int));

/**************************/
/* allocate output memory */
/**************************/
	if (!(out_img = (float *) calloc (nout*vdim, sizeof(float)))) errm (program);
	if (LM_flag && LMout_flag) {
		if (!(LM_out_img = (float *) calloc (LMout*vdim, sizeof (float)))) errm (program); 
	}
	if (res_flag) {
		if (!(res_out_img = (float *) calloc (imgdim[3]*vdim, sizeof (float)))) errm (program); 
	}

/**************************/
/* set output name        */
/**************************/
	sprintf (outroot, "%s%s_dti", imgroot, trailer);

/*****************************************/
/* Mask based on different options       */
/* Add option for input file here!!!     */
/*****************************************/
	/* nblob option */
	if (M_flag) for (i = 0; i < imgdim[2]; i++) {  /* slice loop */

		/* perform the raw threshold */
		for (k = 0; k < slcdim; k++) {  /* pixel in slice loop */
			indx = i*slcdim + k;
			*(img0 + k) = (short int) (in_img[indx] >= I0_thresh * I0_mode);
			*(imsk + indx) = 0;
		}

		/* smooth and find blobs */
		dilate_erode_bin(img0, img1, 13, imgdim[0], slcdim);
		dilate_erode_bin(img1, img0, 4, imgdim[0], slcdim);
		nblob = find_blob(img0, img1, imgdim[0], slcdim);

		/* attempt one repair with higher threshold if needed */
		if (nblob < 0) {
			printf ("%s: attempting to repair mask\n", program);
			for (k = 0; k < slcdim; k++) {  
				indx = i*slcdim + k;
				*(img0 + k) = (short int) (in_img[indx] >= 2.0*I0_thresh*I0_mode);
			}
			dilate_erode_bin(img0, img1, 13, imgdim[0], slcdim);
			dilate_erode_bin(img1, img0, 4, imgdim[0], slcdim);
			nblob = find_blob(img0, img1, imgdim[0], slcdim);
		}

		if (nblob > 0) {
			for (amax=0, j=0; j < nblob; j++) {
				if (objects[j]->m000 > amax) { 
					imax=objects[j]->nblob;
					amax=objects[j]->m000; 
				}
				if (debug) printf ("nblob %d %d %d\n", j, objects[j]->m000,objects[j]->nblob); 
			}
			if (debug) printf ("slice=%d nblob=%d max %d area %d \n", i, nblob, imax, amax);

			/* create the refined threshold, no holes, protect from prior masking */
			for (k = 0; k < slcdim; k++) {  /* pixel in slice loop */
				indx = i*slcdim + k;
				if (*(img1 + k) == imax && in_img[indx] > 0.0) *(imsk + indx) = 1;
			}
		}
		else { /* failed message */
			if (debug) printf ("slice=%d nblob=%d masking failure \n", i, nblob);
		}
	}

	/* classic threshold option */
	else for (i = 0; i < vdim; i++) {
		*(imsk + i) = (short int) (in_img[i] >= I0_thresh * I0_mode);
	}

/*****************************************/
/* calculate the bad encoding file       */
/*****************************************/
	for (i = 0; i < imgdim[2]; i++) {  
		for (j = 0; j < nmeas0; j++)  badenc[i][j] = 0;
	}
	if (B_flag) {
		/* accumulate slice sums */
		for (i = 0; i < imgdim[2]; i++) {  /* slice loop */
			for (j = 0; j < nmeas0; j++) {
				avgenc[i][j] = 0.0;
				numenc[i][j] = 0.0;
				abgenc[i][j] = 0.0;
				nbgenc[i][j] = 0.0;
			}

			for (k = 0; k < slcdim; k++) {  /* pixel in slice loop */
				/* check the threshold in the b=0 image */
				indx = i*slcdim + k;
				if (imsk[indx]) {  /* foreground */
					for (j = 0; j < nmeas0; j++) {
						avgenc[i][j] += in_img[j*vdim + indx];
						numenc[i][j] += 1.0;
					}
				}
				else {       /* avg background */
					for (j = 0; j < nmeas0; j++) {
						abgenc[i][j] += in_img[j*vdim + indx];
						nbgenc[i][j] += 1.0;
					}
				}
			}
		}

		/* calculate slice averages */
		for (i = 0; i < imgdim[2]; i++) {
			for (j = 0; j < nmeas0; j++) {
				if (numenc[i][j] == 0 || avgenc[i][j] < 0.0) avgenc[i][j] = 1.0e-6;
				else avgenc[i][j] /= numenc[i][j];
				if (nbgenc[i][j] == 0 || abgenc[i][j] < 0.0) abgenc[i][j] = 1.0e-6;
				else abgenc[i][j] /= nbgenc[i][j];
			}
		}

		/* calculate linear regression and mark bad encodes for each slice */
		nbadtot = 0;
		for (i = 0; i < imgdim[2]; i++) {
			/* Foreground calcs */
			for (j = 1; j <= nmeas0; j++) {
				xvals[j] = bval[j];
				yvals[j] = log(avgenc[i][j-1]);
				sigs[j] = 1.0/avgenc[i][j-1];
			}

			/* y = a + bx */
			fit(xvals, yvals, nmeas0, sigs, 1.0, &afit, &bfit, &siga, &sigb, &chi2, &qq, &lcc);

			/* background calcs */
			sres = s2res = 0.0; nbad = 0;
			for (j = 1; j <= nmeas0; j++) {
				sres += abgenc[i][j-1];
				s2res += abgenc[i][j-1]*abgenc[i][j-1];
			}
			sres = sres / nmeas0;
			s2res = (s2res - nmeas0*sres*sres)/(nmeas0 - 1.0);
			s2res = sqrt(s2res);

			/* NOTE: loop on encodes starts with 2, so not to discard b=0 */
			if (debug) printf ("slice %3d FG benc: ", i);
			nbadenc = 0;
			for (j = 2; j <= nmeas0; j++) {
				resid2 = fabs((afit + bfit*xvals[j] - yvals[j])/sigs[j]);
				if (resid2 > bstd_thresh*sqrt(chi2/(nmeas0-2.0))) {
					badenc[i][j-1] = 1;
					nbadenc++;
					if (debug) printf("%2d ", j);
				}
			}

			/* NOTE: background threshold, also start with 2nd encode */
			if (b_flag) {
				if (debug) printf ("BG benc:");
				for (j = 2; j <= nmeas0; j++) {
					if (abgenc[i][j-1] > 2.0*sres) {
						badenc[i][j-1] = 2;
						if (debug) printf ("%2d ", j);
						nbadenc++;
					}
				}
			}
			nbadtot += nbadenc;

			if (debug) printf ("total %d out of %d\n", nbadenc, nmeas0);
			if (nmeas0 - nbadenc < 10) {
				fprintf (stderr, "%s: number of bad encodes on slice %d is too large\n", program, i);
				exit (-1);
			}
		}
		printf ("total bad encodes: %d with thresh of %f\n", nbadtot, bstd_thresh);

		/* back up original bval, qval */
		for (j = 1; j <= nmeas0; j++) {
			bval0[j] = bval[j];
			for (k = 1; k <= 3; k++) qval0[j][k] = qval[j][k];
		} 
	}

/********************************************/
/* loop over middle slice to estimate noise */
/********************************************/
	sres = s2res = 0.0; nbad = 0;
	nslc = imgdim[2]/2;
	for (i = 0; i < slcdim; i++) {
		indx = nslc*slcdim + i;

		/* remove bad encodings */
		if (B_flag) {
			nnmeas = 0;
			for (j = 1; j <= nmeas0; j++) {
				if (badenc[nslc][j-1]) continue;
				nnmeas++;
				bval[nnmeas] = bval0[j];
				for (k = 1; k <= 3; k++) 
					qval[nnmeas][k] = qval0[j][k];
				ival[nnmeas] = *(in_img + (j-1)*vdim + indx);
			}
			set_svd (nnmeas, bval, qval, exp_arr);
		}
		else {
			nnmeas = nmeas0;
			for (j = 0; j < nmeas0; j++) {
				ival[j+1] = *(in_img + j*vdim + indx);
			}
		}

		/* process only voxels in mask */
		if (!imsk[indx]) continue;

		/* calculate the eignevalues */
		calc_svd (nnmeas, ival, exp_arr, dvec);
		eigen_calc (dvec, dinv, eigenval, eigenvec, ang);
		Dbar = ang[4];
		Asig = ang[5];
		FA   = ang[6];

		res2 = 0.0;
		for (k = 1; k <= nnmeas; k++) {
			ex = qval[k][1]*dvec[1] + qval[k][2]*dvec[4] + qval[k][3]*dvec[5];
			ey = qval[k][1]*dvec[4] + qval[k][2]*dvec[2] + qval[k][3]*dvec[6];
			ez = qval[k][1]*dvec[5] + qval[k][2]*dvec[6] + qval[k][3]*dvec[3];
			drot = qval[k][1]*ex + qval[k][2]*ey + qval[k][3]*ez;
			res = (ival[k] - exp(dvec[7])*exp(-bval[k]*drot));
			res2 += res*res;
		} 
		res2 = sqrt (res2/nnmeas);

		nbad++;
		sres += res2;
		s2res += res2*res2;
	}

	/* get stats on residues of brain voxels */
	if (nbad > 1) {
		sres = sres / nbad;
		s2res = (s2res - nbad*sres*sres)/(nbad - 1.0);
		s2res = sqrt(s2res);
		printf ("mid slice %d voxels %d mean res %f std %f\n",nslc,nbad,sres,s2res);
	} else {
		fprintf (stderr, "%s: no brain voxels in middle slice %d\n", program, nslc);
		exit (-i);
	}

/**********************************************/
/* main loop over all voxels, calc dti params */
/**********************************************/
	printf ("processing slice:");
	bvox_num = 0; nbad = 0;
	for (i = 0; i < vdim; i++) {
		badvox[i] = 0; 	
		resid[i] = 0.0;
		nslc = i/slcdim;

		/* single voxel mode */
		if (vox_flag && (voxsel != i)) continue;

		/* single slice mode */
		if (s_flag && (slicesel != nslc)) continue;

		/* print new slice count */
		if (!(i % slcdim)) {
			printf (" %d", i/slcdim + 1); fflush (stdout);
		}

		/* remove bad encodings */
		if (B_flag) {
			nnmeas = 0;
			for (j = 1; j <= nmeas0; j++) {
				if (badenc[nslc][j-1]) continue;
				nnmeas++;
				bval[nnmeas] = bval0[j];
				for (k = 1; k <= 3; k++) 
					qval[nnmeas][k] = qval0[j][k];
				ival[nnmeas] = *(in_img + (j-1)*vdim + i);
			}
			set_svd (nnmeas, bval, qval, exp_arr);
		}
		else {
			nnmeas = nmeas0;
			for (j = 0; j < nmeas0; j++) {
				ival[j+1] = *(in_img + j*vdim + i);
			}
		}

		/* process only voxels in mask */
		if (!imsk[i]) continue;

		/* print out voxel number for debug */
		if (pvox_flag) {
			printf (" %d", i); fflush (stdout);
		}

		/* calculate the eignevalues */
		calc_svd (nnmeas, ival, exp_arr, dvec);
		eigen_calc (dvec, dinv, eigenval, eigenvec, ang);
		Dbar = ang[4];
		Asig = ang[5];
		FA   = ang[6];

		bvox_flg = 0;
		for (j=1; j<=3; j++) if (eigenval[j] < 0.0) bvox_flg++;	
		if (Asig > 1.0 || Asig < 0.0) bvox_flg++;
		if (Dbar < 0.0) bvox_flg++;
		if (bvox_flg > 1) bvox_num++;
		badvox[i] = bvox_flg;

		res2 = 0.0;
		indx = 0;
		for (j = 1; j <= nmeas0; j++) {
			if (badenc[nslc][j-1]) continue;
			indx++;
			ex = qval[indx][1]*dvec[1] + qval[indx][2]*dvec[4] + qval[indx][3]*dvec[5];
			ey = qval[indx][1]*dvec[4] + qval[indx][2]*dvec[2] + qval[indx][3]*dvec[6];
			ez = qval[indx][1]*dvec[5] + qval[indx][2]*dvec[6] + qval[indx][3]*dvec[3];
			drot = qval[indx][1]*ex + qval[indx][2]*ey + qval[indx][3]*ez;
			res = (ival[indx] - exp(dvec[7])*exp(-bval[indx]*drot));
			res2 += res*res;
			if (res_flag) res_out_img[(j-1)*vdim + i] = fabs(res);
		} 
		res2 = sqrt (res2/nnmeas);
		resid[i] = res2;


/**********************************************/
/* voxel debug section, print out information */
/**********************************************/
		if (vox_flag) {
			printf("voxel debug = %d\n", voxsel);
			sprintf (Info_file, "%s.txt", outroot);
			if (!(fp_info = fopen (Info_file, "w"))) errw (program, Info_file);

			for (k = 0; k < argc; k++) fprintf (fp_info,"%s ", argv[k]);
			fprintf (fp_info, "\ngradient table:\n");
			for (k = 0; k < nnmeas; k++) fprintf (fp_info, "%10.4f\t%10.4f\t%10.4f\t%10.4f\n",
							bval[k+1], qval[k+1][1], qval[k+1][2], qval[k+1][3]);

			fprintf (fp_info, "\n\nVoxel: %d\n", i);
			k = i;
			iz = k / (imgdim[0]*imgdim[1]); k = k % (imgdim[0]*imgdim[1]);
			iy = k /  imgdim[0]; 		k = k % imgdim[0];
			ix = k;
 			fprintf (fp_info, "\tix= %d\n\tiy= %d\n\tiz= %d\n",ix,iy,iz);

			fprintf (fp_info, "\nIvals:\n");
			for (k = 0; k < nnmeas; k++) fprintf (fp_info, "%3d %10.4f\n",k+1,ival[k+1]);

			fprintf (fp_info, "\n\n\tInitial DTI calc (LLLS)\n");
			fprintf (fp_info,"\tEigenval 1 = %10.4f\n\tEigenval 2 = %10.4f\n\tEigenval 3 = %10.4f\n",
						eigenval[1],eigenval[2],eigenval[3]);
			fprintf (fp_info, "\tMD = %10.4f\n\n",Dbar);
			fprintf (fp_info, "\tRA = %10.4f\n\n",Asig);
			fprintf (fp_info, "\tFA = %10.4f\n\n",FA);

			fflush (fp_info);
			if ((fclose (fp_info))) errw (program, Info_file);	
			printf("Writing info file %s\n", Info_file);
		}

/***************************************************************/
/* LM module: use for LM mode or in repair mode for bad voxels */
/***************************************************************/
		if (LM_flag || (Z_flag && badvox[i] > 0 && nnmeas >= 10)) { 

			if (nonlin_dti (nnmeas, ival, bval, qval, bvox_flg, CO_flag, eigenval, ang, parval, parstd)) {
				printf ("LM failure\n"); fflush (stdout);
				badvox[i] = 9;  /* code for LM failure */
			} else {   /* good exit from LM processor */
				Dbar	= (parval[1] + parval[2] + parval[3])/3.0;
				Asig0	= Asig;
				Asig	= sqrt(((parval[1]-Dbar)*(parval[1]-Dbar)
					+       (parval[2]-Dbar)*(parval[2]-Dbar)
					+	(parval[3]-Dbar)*(parval[3]-Dbar))/6.0)/Dbar;

				for (k = 1; k <= 3; k++) {
					eigenval[k]	= parval[k];
					ang[k]		= parval[3 + k];
				}
				FA = sqrt(((parval[1]-Dbar)*(parval[1]-Dbar)+
					   (parval[2]-Dbar)*(parval[2]-Dbar)+
					   (parval[3]-Dbar)*(parval[3]-Dbar))/(2./3.)) 
					/sqrt (parval[1]*parval[1] + parval[2]*parval[2] + parval[3]*parval[3]);

				/* store LM results in a separate array */
				if (LM_flag && LMout_flag) {
					for (k = 0; k <= MA; k++) (LM_out_img + vdim*k)[i] = parval[k];
					for (k = 0; k <= MA; k++) (LM_out_img + vdim*(k+MA + 1))[i] = parstd[k];
				}
			}

/****************************************/
/* calculate residuals for the LM model */
/****************************************/
			eval2dmat(eigenval[1],eigenval[2],eigenval[3],ang[1],ang[2],ang[3],dmat);
			if (Dtensor) {
				dvec[1] = dmat[1][1];
				dvec[2] = dmat[2][2];
				dvec[3] = dmat[3][3];
				dvec[4] = dmat[1][2];
				dvec[5] = dmat[1][3];
				dvec[6] = dmat[2][3];
			}

			res2 = 0.0;
			indx = 0;
			for (j = 1; j <= nmeas0; j++) {
				if (badenc[nslc][j-1]) continue;
				indx++;
				ex = qval[indx][1]*dmat[1][1] + qval[indx][2]*dmat[1][2] + qval[indx][3]*dmat[1][3];
				ey = qval[indx][1]*dmat[2][1] + qval[indx][2]*dmat[2][2] + qval[indx][3]*dmat[2][3];
				ez = qval[indx][1]*dmat[3][1] + qval[indx][2]*dmat[3][2] + qval[indx][3]*dmat[3][3];
				drot = qval[indx][1]*ex + qval[indx][2]*ey + qval[indx][3]*ez;
				res = (ival[indx] - parval[8] - parval[7]*exp(-bval[indx]*drot));
				res2 += res*res;
				if (res_flag) res_out_img[(j-1)*vdim + i] = fabs(res);
			} 
			res2 = sqrt (res2/nnmeas);
			resid[i] = res2;

			/* store residue for nonlinear output */
			if (LM_flag && LMout_flag) (LM_out_img + vdim*(MA+1))[i] = res2;
		} /* end LM section */

/********************************************/
/* plug calculated values into output array */
/********************************************/
		j = 0;
		(out_img + vdim*j++)[i] = Dbar;
		(out_img + vdim*j++)[i] = Asig;
		if (Dtensor)	for (k = 0; k < 6; k++)	(out_img + vdim*j++)[i] = dvec[k + 1];
		if (FA_flag) 				(out_img + vdim*j++)[i] = FA;
		if (Eigvals) {
			for (k = 0; k < 3; k++)	(out_img + vdim*j++)[i] = eigenval[k + 1];
		 				(out_img + vdim*j++)[i] = 0.5*(eigenval[1]+eigenval[2]);
		}
		if (CO_flag)				(out_img + vdim*j++)[i] = parval[8]/(parval[8]+parval[7]);
		if (Res_flag) 				(out_img + vdim*j++)[i] = resid[i];
		for (l = 0; l < Eigvecs; l++) {
				for (k = 0; k < 3; k++)	(out_img + vdim*j++)[i] = eigenvec[k + 1][3 - l];
		}
		if (Prolat) {	
			prolaticity = (-(eigenval[3]+eigenval[1]) + 2.0*eigenval[2])/(eigenval[3]-eigenval[1]);
			(out_img + vdim*j++)[i] = prolaticity;
		}
		if (debug) (out_img + vdim*j++)[i] = 0.0;  /* to be set as needed */

	} /**** end main pixel loop ******/

	printf ("\n");
	printf ("number of bad voxels %d out of %d\n", bvox_num, vdim);

/**************************/
/* output computed images */
/**************************/
	sprintf (outfile, "%s.4dfp.img", outroot);
	printf ("Writing: %s\n", outfile);
	if (!(fp = fopen (outfile, "wb")) ||
	ewrite (out_img, nout*vdim, control, fp) || fclose (fp)) errw (program, outfile);

/*********************************************/
/* write ifh file and create analyze 7.5 hdr */
/*********************************************/
	imgdim[3] = nout;
	writeifhmce (program, outroot, imgdim, ifh.scaling_factor, ifh.orientation, ifh.mmppix, ifh.center, control);
	sprintf (command, "ifh2hdr -r-1to1 %s", outroot); status = system (command); 

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, rcsid, control);
	printrec ("dti parameter file\n");
	catrec (prmfile);
	sprintf (command, "Bad   voxel count %d/%d\n", bvox_num,    vdim); printrec (command);
	sprintf (command, "Mask threshold = %.4f*I0 mode = %.2f\n", I0_thresh, I0_thresh*I0_mode); printrec (command);
	catrec (imgfile);
	j = 0;		for (k = 0; k < 2; k++) labels[j++] = outtype0[k];
	if (Dtensor)	for (k = 0; k < 6; k++) labels[j++] = outtypeD[k];
	if (FA_flag)	for (k = 0; k < 1; k++) labels[j++] = outtypeF[k];
	if (Eigvals)	for (k = 0; k < 4; k++) labels[j++] = outtypeE[k];
	if (CO_flag)	for (k = 0; k < 1; k++) labels[j++] = outtypeC[k];
	if (Res_flag)	for (k = 0; k < 1; k++) labels[j++] = outtypeR[k];
	for (l = 0; l < Eigvecs; l++) {
			for (k = 0; k < 3; k++) labels[j++] = outtypeV[3*(2 - l) + k];
	}
	if (Prolat) 	for (k = 0; k < 1; k++) labels[j++] = outtypeP[k];
	for (k = 0; k < nout; k++) {
		sprintf (command, "Volume %2d %s\n", k + 1, labels[k]);
		printrec (command);
	}
	endrec ();

/***************************************/
/* output computed LM-specific images */
/***************************************/
	if (LM_flag && LMout_flag) {
		sprintf (LMoutroot, "%s%s_LMdti", imgroot, trailer);
		sprintf (LMoutfile, "%s.4dfp.img", LMoutroot);
		printf ("Writing: %s\n", LMoutfile);
		if (!(fp = fopen (LMoutfile, "wb")) 
		|| ewrite (LM_out_img, LMout*vdim, control, fp) 
		|| fclose (fp)) errw (program, LMoutfile);

	        imgdim[3] = LMout;
		writeifhmce (program, LMoutroot, imgdim, ifh.scaling_factor, ifh.orientation, ifh.mmppix, ifh.center, control);
		sprintf (command, "ifh2hdr -r1000 %s", LMoutroot); status = system (command);

		startrece (LMoutfile, argc, argv, rcsid, control);
		catrec (outfile);
		for (k = 0; k < 2*(MA + 1); k++) {
			sprintf (command, "Volume %2d %s\n", k + 1, outtypeN[k]);
			printrec (command);
		}
		endrec ();
	}

/**********************************/
/* output computed residue images */
/**********************************/
	if (res_flag) {
		sprintf (RESoutroot, "%s%s_res", imgroot, trailer);
		sprintf (RESoutfile, "%s.4dfp.img", RESoutroot);
		printf ("Writing: %s\n", RESoutfile);
		if (!(fp = fopen (RESoutfile, "wb")) 
		|| ewrite (res_out_img, nmeas0*vdim, control, fp) 
		|| fclose (fp)) errw (program, RESoutfile);

	        imgdim[3] = nmeas0;
		writeifhmce (program, RESoutroot, imgdim, ifh.scaling_factor, ifh.orientation, ifh.mmppix, ifh.center, control);
		sprintf (command, "ifh2hdr -r1000 %s", RESoutroot); status = system (command);

		startrece (RESoutfile, argc, argv, rcsid, control);
		catrec (outfile);
		sprintf (command, "Absolute value residuals for each encoding");
		printrec (command);
		endrec ();
	}

	free (in_img); free (out_img);
	free (img0); free (img1); free (imsk);
	free_vector (ival, 1, nmeas0);
	free_vector (bval, 1, nmeas0);
	free_vector (eigenval, 1, 3);
	free_vector (ang, 1, 9);
	free_vector (dvec,1, 7);
	free_matrix (exp_arr, 1, nmeas0, 1, 7);
	free_matrix (eigenvec, 1, 3, 1, 3);
	free_matrix (dinv, 1, 3, 1, 3);
	free_matrix (qval, 1, nmeas0, 1, 3);
	free_ivector(badvox, 0, vdim);
	free_vector(resid, 0, vdim);
	free_matrix(dmat, 1, 3, 1, 3);
	free_matrix(dmatd, 1, 3, 1, 3);	

	if (LM_flag && LMout_flag) {
		free (LM_out_img);
	}
	if (res_flag) {
		free (res_out_img);
	}

	free_imatrix(badenc, 0, imgdim[2], 0, nmeas0);
	if (B_flag) {
		free_matrix (avgenc, 0, imgdim[2], 0, nmeas0);
		free_matrix (numenc, 0, imgdim[2], 0, nmeas0);
		free_matrix (abgenc, 0, imgdim[2], 0, nmeas0);
		free_matrix (nbgenc, 0, imgdim[2], 0, nmeas0);
		free_vector (xvals, 1, nmeas0);
		free_vector (yvals, 1, nmeas0);
		free_vector (sigs, 1, nmeas0);
		free_vector (bval0, 1, nmeas0);
		free_matrix (qval0, 1, nmeas0, 1, 3);
	}
	exit (status);
}

/* recreate the diffution tensor matrix for residual calculations */
void eval2dmat(float eval1, float eval2, float eval3, float ang1, float ang2, float ang3, float **dmat)
{
	int k,l,m,n;
	float sa,sb,sg,ca,cb,cg;
	float eigenvec[4][4], dmatd[4][4];

	sa = sin((double) ang1); ca = cos((double) ang1);
	sb = sin((double) ang2); cb = cos((double) ang2);
	sg = sin((double) ang3); cg = cos((double) ang3);
	
	/* form rotation matrix Arfken p.180 */
	eigenvec[1][1] = cg*cb*ca - sg*sa;
	eigenvec[2][1] = cg*cb*sa + sg*ca;
	eigenvec[3][1] = -cg*sb;
	eigenvec[1][2] = -sg*cb*ca - cg*sa;
	eigenvec[2][2] = -sg*cb*sa + cg*ca;
	eigenvec[3][2] = sg*sb;
	eigenvec[1][3] = sb*ca;
	eigenvec[2][3] = sb*sa;
	eigenvec[3][3] = cb;	

	/* form the diagonal matrix */
	for (m = 1; m <= 3; m++) for (n = 1; n <= 3; n++) dmatd[m][n] = 0.0;
	dmatd[1][1] = eval1;
	dmatd[2][2] = eval2;
	dmatd[3][3] = eval3;
	
	/* rotate the tensor into the lab frame */
	for (m = 1; m <= 3; m++) for (n = 1; n <= 3; n++) {
		dmat[m][n] = 0.0;
		for (k = 1; k <= 3; k++) for (l = 1; l <= 3; l++) {
			dmat[m][n] += eigenvec[m][k]*dmatd[k][l]*eigenvec[n][l];
		}
	}

}

/* fit to a straight line, NR page 665 added linear cor coeff */
/* linear cor coeff is not the same as the Pearson R coefficient */
void fit(float x[], float y[], int ndata, float sig[], int mwt, float *a,
	float *b, float *siga, float *sigb, float *chi2, float *q, float *lcc)
{
	/* modified to return dummy q
	float gammq(float a, float x);
	*/
	int i;
	float wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
	
	*b=0;
	if (mwt) {
		ss=0.0;
		for (i=1;i<=ndata;i++) {
			wt=1.0/(sig[i]*sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	}
	else {
		for (i=1;i<=ndata;i++) {
			sx += x[i];
			sy += y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt) {
		for (i=1;i<=ndata;i++) {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
		}
	}
	else {
		for (i=1;i<=ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
	}
	*b /= st2;
	*a = (sy-sx*(*b))/ss;
	*siga = sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb = sqrt(1.0/st2);
	*lcc = -sx/(ss*st2*(*siga)*(*sigb));
	*chi2 = 0.0;
	if (mwt==0) {
		for (i=1;i<=ndata;i++)
			*chi2 += (y[i]-(*a)-(*b)*x[i])*(y[i]-(*a)-(*b)*x[i]);
		*q=1.0;
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *= sigdat;
		*sigb *= sigdat;
	}
	else {
		for (i=1;i<=ndata;i++) 
			*chi2 += ((y[i]-(*a)-(*b)*x[i])/sig[i])*((y[i]-(*a)-(*b)*x[i])/sig[i]);
			/*
			*q = gammq(0.5*(ndata-2),0.5*(*chi2)); 
			*/
			*q = 1.0;
	}
}
