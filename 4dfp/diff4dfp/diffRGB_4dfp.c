/*$Header: /home/usr/shimonyj/diff4dfp/RCS/diffRGB_4dfp.c,v 1.13 2010/06/16 23:50:31 avi Exp $*/
/*$Log: diffRGB_4dfp.c,v $
 * Revision 1.13  2010/06/16  23:50:31  avi
 * change gbr option flag to 'G' to maintain compatibility with diff_4dfp
 *
 * Revision 1.12  2010/05/10  03:42:25  avi
 * option -B (change color code to BGR) - compensates for B<->R reversal in analyze ps output
 *
 * Revision 1.11  2010/03/30  22:36:47  avi
 * correct error checking in case of tensor data
 *
 * Revision 1.10  2008/09/16  03:16:01  avi
 * linux compliant
 * #include <JSSutil.h>
 *
 * Revision 1.9  2006/09/15  06:26:44  adrian
 * option -D (bypass D tensor calculation)
 *
 * Revision 1.8  2006/09/04  06:25:08  avi
 * updated TRX; eliminate local use of ifh; correctly free_matrix and free_vector;
 *
 * Revision 1.7  2006/02/21  05:06:12  adrian
 * install calls to t4tolin() and affine_DTrot() to deal with affine transformed DWI data
 *
 * Revision 1.6  2004/12/07  22:09:16  avi
 * do not skip tensor computation if an input image has negative intensity
 *
 * Revision 1.5  2004/06/10  03:40:04  avi
 * substitute fflip.f -> cflip.c routines
 * eliminate FORTRAN calls
 *
 * Revision 1.4  2004/06/09  00:22:18  avi
 * more sophisticated I0 histogram mode finding
 * use float fimg_modetn (float* fimg, int nval, float histmin, int nsmooth);
 *
 * Revision 1.3  2003/10/02  02:06:16  avi
 * -q (sqrt(Asigma)) option
 *
 * Revision 1.2  2003/10/01  05:51:07  avi
 * correct usage
 *
 * Revision 1.1  2003/09/17  22:19:54  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <JSSutil.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>

/*********************************/
/* size of input and output sets */
/*********************************/
#define MAXL		256
#define I0_THRESH	0.1		/* threshold as fraction of I0 image mode value */
#define HISTMIN		100.		/* minimum value for fimg_modetn() to find I0 image mode value */
#define NSMOOTH		4		/* fimg_modetn() smoothing parameter */

/***************/
/* fimg_mode.c */
/***************/
float	fimg_modetn (float* fimg, int nval, float histmin, int nsmooth);
float	fimg_mode  (float*, int);

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

extern void	flipx (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */
extern void	flipy (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */
extern void	flipz (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */
extern int	Inithdr (struct dsr *phdr, int *imgdim, float *voxdim, char *proto_header);	/* Inithdr.c */

static int	logical (int l) {
	return (l) ? 1 : 0;
}

static char rcsid[] = "$Id: diffRGB_4dfp.c,v 1.13 2010/06/16 23:50:31 avi Exp $";
int main (int argc, char *argv[]) {
/**************/
/* processing */
/**************/
	char		prmfile[MAXL];
	int		n, bvox_flg, bvox_num;
	double		Asig, prolaticity, Dbar, q;
	float		I0_mode = 500., I0_thresh = I0_THRESH, histmin = HISTMIN;
	float		scalei = 1.;
	int		nsmooth = NSMOOTH;

/*******************/
/* i/o 4dfp images */
/*******************/
	struct dsr	hdr;		/* ANALYZE hdr */
	FILE		*fp, *fpimg, *fpout;
	char		imgfile[MAXL], outfile[MAXL], t4file[MAXL] = "";
	char		imgroot[MAXL], outroot[MAXL];

/****************/
/* image arrays */
/****************/
	unsigned char	*imgo;
	float		*imgf;
	float		voxdim[3];
	float		ValueIn = 0., ValueOut = 0., Vmax = 500.;
	int		imgdim[4], vdim;
	int		orient, isbig;
	char		control = '\0';
	
/*************************************/
/* memory for diffusion calculations */
/*************************************/
	float		*ival, *bval, *ang, *dvec, *eigenval;
	float		**qval, **eigenvec, **reigenvec, **dinv, **exp_arr, **lin;

/***********/
/* utility */
/***********/
	char		*str, command[MAXL], program[MAXL];
	int		c, i, j, k, l;
	int		ix, iy, iz, iv, rgb;

/*********/
/* flags */
/*********/
	int		thresh_flg;
	int		sqrt_flag = 0;	/* when set intesity scaled by sqrt(Asigma) instead of Asigma */
	int 		tensor_data = 0;
	int		swab_flag;
	int		bgr_flag = 0;

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
				case 'D': tensor_data++; 		break;
				case 'q': sqrt_flag++;			break;
				case 'G': bgr_flag++;			break;
				case 'n': nsmooth = atoi (str);		*str = '\0'; break;
				case 'h': histmin = atof (str);		*str = '\0'; break;
				case 't': I0_thresh = atof (str);	*str = '\0'; break;
				case 'T': strcpy (t4file, str);		*str = '\0'; break;
				case 'c': scalei = atof (str);		*str = '\0'; break; 
				case '@': control = *str++;		*str = '\0'; break;
			}
		} else switch (k) {
			case 0: getroot (argv[i], prmfile); k++; break;
			case 1:	getroot (argv[i], imgroot); k++; break;
		}	
	}
	if (k != 2) {
		printf ("Usage:\t%s <prm_file> <file_4dfp>\n", program);
		printf (" e.g.,\t%s -t0.5 -qc1.7 tp7_params.dat /data/DTI_avg\n", program);
		printf ("\toption\n");
		printf ("\t-q\tscale intesity by sqrt(Asig) instead of Asig\n");
		printf ("\t-G\tchange color coding to bgr (default rgb)\n");
		printf ("\t-c<flt>\tspecify the intensity scale value (default=%.4f)\n", scalei);
		printf ("\t-t<flt>\tspecify mask threshold as fraction of I0 mode (default=%.4f)\n", I0_THRESH);
		printf ("\t-T<str>\tspecify t4 file used to transform DWI data\n");
		printf ("\t-h<flt>\tspecify minimum I0 mode (default=%.2f)\n", HISTMIN);
		printf ("\t-n<int>\tspecify number of I0 histogram smoothings (default=%d)\n", NSMOOTH);
		printf ("\t-D\tinput <file_4dfp> is 8 volume diff_4dfp -D output (Dbar, Asigma, D tensor)\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:\t<prm_file> is ignored with -D option\n");
		exit (1);
	}

/******************************/
/* get number of measurements */
/******************************/
	if (tensor_data) {
		n = 8;
	} else {
		if ((n = dti_dimen (prmfile)) <= 0) errr (program, prmfile);
	}

/**************************************/
/* allocate memory for diffusion calc */
/**************************************/
	ival		= vector(1, n);
	bval		= vector(1, n);
	eigenval	= vector(1, 3);
	ang		= vector(1, 9);
	dvec		= vector(1, 7);
	if (!ival || !bval || !eigenval || !ang || !dvec) errm (program);

	exp_arr		= matrix(1, n, 1, 7);
	eigenvec	= matrix(1, 3, 1, 3);
	reigenvec	= matrix(1, 3, 1, 3); 	/* rotated eigenvectors */
	lin		= matrix(1, 3, 1, 3);	/* linear transform derived from t4 file */
	dinv		= matrix(1, 3, 1, 3);
	qval		= matrix(1, n, 1, 3);
	if (!exp_arr || !eigenvec || !reigenvec || !lin || !dinv || !qval) errm (program);

	if (!tensor_data) {
/*****************************************/
/* get the b values and the unit vectors */
/*****************************************/
		get_dti_params_nrc (prmfile, n, bval, qval);

/********************************************/
/* set up the svd for diffusion calculation */
/********************************************/
		set_svd (n, bval, qval, exp_arr);
	}

/*******************************************/
/* get stack dimensions and optional t4file*/
/*******************************************/
	if (get_4dfp_dimoe (imgroot, imgdim, voxdim, &orient, &isbig)) errr (program, imgroot);
	if (!control) control = (isbig) ? 'b' : 'l';
	if (imgdim[3] != n && !tensor_data) {
		fprintf (stderr, "%s: %s %s dimension mismatch\n", program, prmfile, imgroot);
		exit (-1);
	}
	if (imgdim[3] < 8 && tensor_data) {
		fprintf (stderr, "%s: %s has less than 8 volumes\n", program, imgroot);
		exit (-1);
	}
	if (strlen (t4file)) t4tolin (t4file, orient, lin);

/*********************************************************************/
/* get input 4dfp dimensions and allocate input memory and get stack */
/*********************************************************************/
	vdim = imgdim[0]*imgdim[1]*imgdim[2];
	imgf = (float *) calloc (n*vdim, sizeof (float));
	imgo = (unsigned char *) calloc (3*vdim, sizeof (unsigned char));
	if (!imgf || !imgo) errm (program);
	
	printf ("stack dimensions %10d%10d%10d%10d\n", imgdim[0], imgdim[1], imgdim[2], imgdim[3]);
	printf ("voxel size (mm)  %10.6f%10.6f%10.6f\n", voxdim[0], voxdim[1], voxdim[2]);

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	printf ("Reading: %s\n", imgfile);
	if (!(fpimg = fopen (imgfile, "rb")) || eread (imgf, n*vdim, isbig, fpimg)
	|| fclose (fpimg)) errr (program, imgfile);

/************************/
/* flip 4dfp to analyze */
/************************/
	for (iv = 0; iv < n; iv++) switch (orient) {
		case 4:	flipx (imgf + iv*vdim, imgdim + 0, imgdim + 1, imgdim + 2); /* sagittal */
		case 3:	flipz (imgf + iv*vdim, imgdim + 0, imgdim + 1, imgdim + 2); /* coronal */
		case 2:	flipy (imgf + iv*vdim, imgdim + 0, imgdim + 1, imgdim + 2); /* transverse */
		break;
		default: 
		fprintf (stderr, "%s: illegal %s orientation (%d)\n", program, imgfile, orient);
		exit (-1);
	}

	if (!tensor_data) {
/*******************/
/* compute I0 mode */
/*******************/
		I0_mode = fimg_modetn (imgf, vdim, histmin, nsmooth);
		printf ("%s computed I0 mode %.2f\n", imgfile, I0_mode);
	}

	sprintf (outroot, "%s_%s", imgroot, (bgr_flag) ? "BGR" : "RGB");
	sprintf (outfile, "%s.img", outroot);
	printf ("Writing: %s\n", outfile);
	if (!(fpout = fopen (outfile, "w"))) errw (program, outfile);
/****************************************/
/* loop over voxels and calc dti params */
/****************************************/
	bvox_num = 0;
	printf ("processing slice:");
	for (iz = 0; iz < imgdim[2]; iz++) {printf (" %d", iz+1); fflush (stdout);
	for (iy = 0; iy < imgdim[1]; iy++) {
	for (ix = 0; ix < imgdim[0]; ix++) {
		i = ix + imgdim[0]*(iy + imgdim[1]*iz); 
		if (!tensor_data) {	
			for (bvox_flg = j = 0; j < n; j++) {
				if ((ival[j+1] = *(imgf + j*vdim + i)) <= 0.) bvox_flg++;
			}
			if (bvox_flg) bvox_num++;
			thresh_flg = (imgf[i] < I0_thresh * I0_mode);
			if (thresh_flg) continue;
			calc_svd (n, ival, exp_arr, dvec);
			eigen_calc (dvec, dinv, eigenval, eigenvec, ang);	
			Dbar = ang[4]; if (Dbar < 0.) {++bvox_num; continue;}
			Asig = ang[5]; if (Asig < 0.) Asig = 0.; if (Asig > 1.) Asig = 1.;
		} else {
			Dbar = (imgf + 0*vdim)[i];
			Asig = (imgf + 1*vdim)[i];
			for (k = 2; k < 8; k++) dvec[k - 1] = (imgf + vdim*k)[i];
			eigen_calc (dvec, dinv, eigenval, eigenvec, ang);
		}

/***************************************/
/* accommodate affine transformed data */
/***************************************/
		if (strlen (t4file)) {
			affine_DTrot (lin, eigenvec, reigenvec);			
		} else {
			for (k = 1; k <= 3; k++) for (i = 1; i <= 3; i++) reigenvec[k][i] = eigenvec[k][i];
		}

/***************/
/* process rgb */
/***************/
		for (rgb = 0; rgb < 3; rgb++) {
			ValueIn = reigenvec[rgb+1][3];
			q = (sqrt_flag) ? sqrt (Asig) : Asig;
			ValueOut = fabs (255*q*scalei*ValueIn); if (ValueOut > 255.) ValueOut = 255.;
			k = (bgr_flag) ? 2 - rgb : rgb;
			imgo[ix + imgdim[0]*(iy + imgdim[1]*(3*iz + k))] = (unsigned char) ValueOut;
		}
	}}}
	printf ("\n");
	if (fwrite (imgo, sizeof (unsigned char), 3*vdim, fpout) != 3*vdim) errw (program, outfile);
	fclose (fpout);
	printf ("bad voxel count %d out of %d\n", bvox_num, vdim);

/**********************/
/* create ANALYZE hdr */
/**********************/
	Inithdr (&hdr, imgdim, voxdim, "");
	hdr.dime.dim[4] = 1;
	hdr.dime.bitpix = 24;		/* 3 * 8 */
	hdr.dime.datatype = 128;	/* RGB */
	hdr.hist.orient = orient - 2;
	hdr.dime.glmin = 0;
	hdr.dime.glmax = 255;

	printf ("control=%c\n", control);
	swab_flag = ((CPU_is_bigendian() != 0) && (control == 'l' || control == 'L'))
		 || ((CPU_is_bigendian() == 0) && (control == 'b' || control == 'B'));
	printf ("swab_flag=%d\n", swab_flag);
	if (swab_flag) swab_hdr (&hdr);

	sprintf (outfile, "%s.hdr", outroot);
	printf ("Writing: %s\n", outfile);
	if (!(fpout = fopen (outfile, "wb")) || fwrite (&hdr, sizeof (struct dsr), 1, fpout) != 1) errw (program, outfile);
	fclose (fpout);

/*******************/
/* create rec file */
/*******************/
 	sprintf (outfile, "%s.img", outroot);
	startrece (outfile, argc, argv, rcsid, control);
	printrec ("dti parameter file\n");
	catrec (prmfile);
	if (strlen (t4file)) {
		printrec ("affine transform file\n");
		catrec (t4file);
	}
	sprintf (command, "Bad voxel franction %d/%d\n", bvox_num, vdim); printrec (command);
	sprintf (command, "Mask threshold = %.4f*I0 mode = %.2f\n", I0_thresh, I0_thresh*I0_mode); printrec (command);
	catrec (imgfile);
	endrec ();

	free_vector (ival,	1, n);
	free_vector (bval,	1, n);
	free_vector (eigenval,	1, 3);
	free_vector (ang,	1, 9);
	free_vector (dvec,	1, 7);

	free_matrix (exp_arr,	1, n, 1, 7);
	free_matrix (eigenvec,	1, 3, 1, 3);
	free_matrix (reigenvec,	1, 3, 1, 3);
	free_matrix (lin,	1, 3, 1, 3);
	free_matrix (dinv,	1, 3, 1, 3);
	free_matrix (qval,	1, n, 1, 3);

	free (imgo);
	free (imgf);
	exit (0);
}
