/****************************************************************************************/
/* Copyright 2003, 2004, 2005, 2006, 2007, 2008, 2009					*/
/* Washington University, Mallinckrodt Institute of Radiology.   			*/
/* All Rights Reserved.									*/
/* This software may not be reproduced, copied, or distributed without written		*/
/* permission of Washington University. For further information contact A. Z. Snyder.	*/
/****************************************************************************************/
/*$Header: /data/petsun4/data1/src_solaris/dwi_xalign_4dfp/RCS/dwi_cross_xalign3d_4dfp.c,v 1.13 2009/02/20 01:15:44 avi Exp $*/
/*$Log: dwi_cross_xalign3d_4dfp.c,v $
 * Revision 1.13  2009/02/20  01:15:44  avi
 * option -I
 *
 * Revision 1.12  2008/07/08  22:41:12  avi
 * capture volume-specific total transforms in output data file instead of stdout
 * do not include data file in rec
 *
 * Revision 1.11  2008/07/01  03:31:53  avi
 * pass on 2D vs. 3D status using modes[] instead of &mo3d in call to resample_()
 *
 * Revision 1.10  2008/03/13  00:46:22  avi
 * option -p (2D alignment mode)
 *
 * Revision 1.9  2007/08/25  03:30:04  avi
 * correct failure to zero imgt (accumulator) in append mode
 *
 * Revision 1.8  2007/07/24  05:56:57  avi
 * optioins -@ (output endian control) and -a (append)
 *
 * Revision 1.7  2007/07/22  03:00:47  avi
 * Solaris 10
 *
 * Revision 1.6  2005/12/30  00:25:49  avi
 * command line control of parameter search radius (options -i and -j; fimgeta_() -> fimgetae_())
 *
 * Revision 1.5  2005/11/29  08:08:33  avi
 * list of input files (with optional weights) capability
 * MAXF -> 64
 *
 * Revision 1.4  2004/12/14  01:32:50  avi
 * msk0 computation accepts 1.0 as in mask
 *
 * Revision 1.3  2003/10/02  23:52:37  avi
 * zoom (-z) functionality
 *
 * Revision 1.2  2003/10/01  02:53:45  avi
 * -nflag (zero negative values)
 *
 * Revision 1.1  2003/05/04  00:08:33  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>
#include <librms.h>

#define MAXL		256
#define MAXF		64		/* maximum number of input dwi files */
#define FSMOOTH		0.0		/* pre-filter */
#define DEL		5.0		/* sampling interval in mm */
#define EPSMM		3.0		/* fimgetae_() translation parameter search radius in mm */
#define EPSRAD		40.0		/* fimgetae_() eps head radius in mm */

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

void errd (char *program, char *filesp1, char *filesp2) {
	fprintf (stderr, "%s: %s %s dimension, orientation or endianness mismatch\n", program, filesp1, filesp2);
	exit (-1);
}

void errn (char *program) {
	fprintf (stderr, "%s: maximum number of input files (%d) exceeded\n", program, MAXF);
	exit (-1);
}

/************/
/* eternals */
/************/
extern void	t4_init_  (float *a);							/* fomega.f */
extern void	t4_omega_ (long *mode, float *t4, float *omega);			/* fomega.f */
extern void	splineza_ (float *imgt, int *nx, int *ny, int *nz, float *d2zi);	/* spline3dvgh.f */
extern void	splinezf_ (float *imgt, int *nx, int *ny, int *nz, float *d2zi);	/* spline3dvgh.f */
extern void	splinex_  (float *imgt, int *nx, int *ny, int *nz, float *d2xi);	/* spline3dvgh.f */
extern void	spliney_  (float *imgt, int *nx, int *ny, int *nz, float *d2yi);	/* spline3dvgh.f */
extern void	fimgetae_ (float *img1, float *d2x1, float *d2y1, float *d2z1,
			   int *nx1, int *ny1, int *nz1, float *mmppix1, float *center1, short *msk1, 
			   long *mode, float *del, float *eps,
			   float *img2, float *d2x2, float *d2y2, float *d2z2,
			   int *nx2, int *ny2, int *nz2, float *mmppix2, float *center2, short *msk2,
			   float *t4);							/* fimgetae.f */
extern void	etagen_   (float *img1, float *d2x1, float *d2y1, float *d2z1,
			   int *nx1, int *ny1, int *nz1, float *mmppix1, float *center1, short *msk1,
			   long *mode, float *omega, float *del,
			   float *img2, float *d2x2, float *d2y2, float *d2z2,
			   int *nx2, int *ny2, int *nz2, float *mmppix2, float *center2, short *msk2,
			   float *eta, float *q);						/* etadeta.f */
extern void	resample_ (float *img1, float *d2x1, float *d2y1, float *d2z1,
			   int *nx1, int *ny1, int *nz1, float *mmppix1, float *center1, short *msk1,
			   long *mode, float *t4, float *img2,
			   int *nx2, int *ny2, int *nz2, float *mmppix2, float *center2);	/* etadeta.f */


void	filter3d (char *program, float *imgt, int *imgdim, float *voxdim, float f0, int wrap);	/* below */
void	spline3d (char *program, float *imgf, int *imgdim);					/* below */
void	volalign (float *tar, short *tarm, char *tarstr, float *src, short *srcm, char *srcstr,	/* below */
		  long *modes, float *t4, float *peta, float *pq);
void	load_filter_spline (char *srcroot, int ivol, float *imag);				/* below */

typedef struct {
	float	t4[16];
	float	eta;
	float	q;
	float	det;
	float	weight;
} VOLDAT;

/********************/
/* global variables */
/********************/
char		program[MAXL];
long		mo3d = 259;			/* enable 3D cross_modal */
int		wrap = 0;			/* flag used by filter3d */
float		fsmooth = FSMOOTH;		/* 3D pre-filter half freq (1/mm) */
float		del = DEL;			/* sampling interval */
float		eps[12];			/* fimgetae_() parameter search radii */
float		mmppixs[3], centers[3];
int		srcsiz, srcdim[4], orients, isbigs;
char		control = '\0';

static char rcsid[] = "$Id: dwi_cross_xalign3d_4dfp.c,v 1.13 2009/02/20 01:15:44 avi Exp $";
int main (int argc, char *argv[]) {
/************/
/* 4dfp I/O */
/************/
	FILE		*fp_out, *fp_mat, *fp_dat, *fp_rec, *fp_src, *fp_lst;
	char		filroo0[MAXL], filroot[MAXL], mskroot[MAXL] = "", srcroot[MAXF + 1][MAXL];
	char		mskfile[MAXL], srcfile[MAXL], lstfile[MAXL] = "";
	char		outfile[MAXL], geofile[MAXL], matfile[MAXL], datfile[MAXL], recfile[MAXL];

/****************/
/* image arrays */
/****************/
	float		*imgs, *img0, *img1, *imgw, *imgt;
	short		*msk0, *all1;
	float		mmppixm[3], mmppixt[3], mmppixo[3];
	float		zoomfac[3] = {1.0, 1.0, 1.0};
	int		outsiz, outdim[4], mskdim[4], tmpdim[4], orientm, orientt;
	int		isbigm, isbigt;

/***********/
/* utility */
/***********/
	char		*str, string[MAXL], command[MAXL], *srgv[MAXL];
	int		c, i, j, k, m, irun, nrun, ivol, four = 4;
	int		nzero = 0;

/**************/
/* processing */
/**************/
	long		modes[8];
	float		xenc_t4[16], total_t4[16];
	float		det;
	float		epsmm = EPSMM, epsrad = EPSRAD;
	VOLDAT		rundat[MAXF];
	int		I0vol = 0;

/*********/
/* flags */
/*********/
	int		append = 0;
	int		stretch = 0;
	int		status = 0;
	int		debug = 0;
	int		use_geo = 0;
	int		planar = 0;
	int		nflag = 0;	/* zero negative output */

	printf ("%s\n", rcsid);
	if (!(str = strrchr (argv[0], '/'))) str = argv[0]; else str++;
	strcpy (program, str);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (string, argv[i]); str = string;
			while (c = *str++) switch (c) {
				case 'g': use_geo++;			break;
				case 's': stretch++;			break;
				case 'n': nflag++;			break;
				case 'p': planar++;			break;
				case 'a': append++;			break;
				case 'w': mo3d |= 1024; wrap++; 	break;
				case 'l': strcpy (lstfile, str);
/*************************/
/* parse input list file */
/*************************/
				if (!(fp_lst = fopen (lstfile, "r"))) errr (program, lstfile);
				while (fgets (command, MAXL, fp_lst)) {
					if (!(m = split (command, srgv, MAXL))) continue;
					rundat[k].weight = 1.;
					getroot (srgv[0], srcroot[k]);
					for (j = 1; j < m; j++) {
						if (strstr (srgv[j], "weight=")) rundat[k].weight = atof (srgv[j] + 7);
					}
					if (++k >= MAXF) errn (program);
				}
				fclose (fp_lst);			*str = '\0'; break;
				case 'd': del = atof (str);		*str = '\0'; break;
				case 'i': epsmm = atof (str);		*str = '\0'; break;
				case 'j': epsrad = atof (str);		*str = '\0'; break;
				case 'f': fsmooth = atof (str);		*str = '\0'; break;
				case 'm': getroot (str, mskroot);	*str = '\0'; break;
				case 'z': switch (*str++) {
						case 'x': zoomfac[0] = atof (str); break;
						case 'y': zoomfac[1] = atof (str); break;
						case 'z': zoomfac[2] = atof (str); break;
					  }				*str = '\0'; break;
				case 'I': I0vol = atoi (str) - 1;	*str = '\0'; break;
				case '@': control = *str++;		*str = '\0'; break;
			}
		} else {
			if (k >= MAXF) errn (program);
			rundat[k].weight = 1.;
			getroot (argv[i], srcroot[k++]);
		}	
	}
	nrun = k - 1;
	if (nrun < 1) {
		printf ("Usage:\t%s <(4dfp) dwi1> <(4dfp) dwi2> <(4dfp) dwin> ... <(4dfp) dwi_out>\n", program);
		printf ("e.g.:\t%s -sgmjo_sub2-dwi1_mskt jo_sub2-dwi1 jo_sub2-dwi2 jo_sub2-dwi_all\n", program);
		printf ("e.g.:\t%s -sgmjo_sub2-dwi1_mskt -ljo_sub2_dwi.lst jo_sub2-dwi_all\n", program);
		printf ("\toption\n");
		printf ("\t-p\tplanar (2D; disable cross-slice) alignment\n");
		printf ("\t-w\tenable wrap addressing\n");
		printf ("\t-s\tenable cross DWI voxel size adjust (principal axis stretch)\n");
		printf ("\t-n\tzero negative values in output image\n");
		printf ("\t-z<x|y|z><flt>\tzoom output x y or z dimension by specified factor\n");
		printf ("\t-g\tuse group geometric mean (*_geom) volumes for cross-run registration\n");
		printf ("\t-a\tappend successive runs in output (default average)\n");
		printf ("\t-m<(4dfp) mask>\tspecify first volume mask\n");
		printf ("\t-I<int>\tspecify volume number of I0 counting from 1 (default 1)\n");
		printf ("\t-f<flt>\tspecify pre-blur filter half freq (1/mm) (default none)\n");
		printf ("\t-d<flt>\tspecify sampling interval in mm (default=%.4f)\n", DEL);
		printf ("\t-i<flt>\tspecify displacment search radius in mm (default=%.4f)\n", EPSMM);
		printf ("\t-j<flt>\tspecify parameter search object radius in mm (default=%.4f)\n", EPSRAD);
		printf ("\t-l<lst>\tread input file names from specified file (use before naming output)\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:\toption -I (non-default I0 volume) must be matched according to use of option -g\n");
		exit (1);
	}

	if (get_4dfp_dimoe (srcroot[0], srcdim, mmppixs, &orients, &isbigs)) errr (srcroot[0], program);
	if (!control) control = (isbigs) ? 'b' : 'l';
	for (outsiz = srcsiz = 1, k = 0; k < 3; k++) {
		centers[k] = mmppixs[k] * (srcdim[k] / 2);
		mmppixo[k] = mmppixs[k] / zoomfac[k];
		outdim[k] = (int) (0.5 + srcdim[k] * zoomfac[k]);
		srcsiz *= srcdim[k];
		outsiz *= outdim[k];
	}
	outdim[3] = srcdim[3];
	if (!(img1 = (float *) malloc (srcsiz * sizeof (float) * 4)))	errm (program);
	if (!(img0 = (float *) malloc (srcsiz * sizeof (float) * 4)))	errm (program);	/* I0 */
	if (!(imgs = (float *) malloc (srcsiz * sizeof (float))))	errm (program);	/* DWI volume */
	if (!(msk0 = (short *) malloc (srcsiz * sizeof (short))))	errm (program);
	if (!(all1 = (short *) malloc (srcsiz * sizeof (short))))	errm (program);
	for (i = 0; i < srcsiz; i++) all1[i] = 1;

/************/
/* set mask */
/************/
        if (!strlen (mskroot)) {
		for (i = 0; i < srcsiz; i++) msk0[i] = 1;
	} else {
		if (get_4dfp_dimoe (mskroot, mskdim, mmppixm, &orientm, &isbigm)) errr (mskroot, program);
		status = orientm != orients;
		for (k = 0; k < 3; k++) status |= (mskdim[k] != srcdim[k] || fabs (mmppixs[k] - mmppixm[k]) > 1.e-4);
		if (status) errd (program, srcroot[0], mskroot);
		load_4dfp_frame (mskroot, srcdim, 0, isbigm, imgs);
		for (i = 0; i < srcsiz; i++) msk0[i] = (imgs[i] >= 1.0) ? 1 : 0;
	}
	strcpy (filroo0, srcroot[0]); if (use_geo) strcat (filroo0, "_geom");
	load_4dfp_frame (filroo0, srcdim, 0, isbigs, imgs);
	for (i = 0; i < srcsiz; i++) if (imgs[i] < 1.0) msk0[i] = 0;

/**************************************/
/* allocate and initialize transforms */
/**************************************/
	for (irun = 0; irun < nrun; irun++) t4_init_ (rundat[irun].t4);
	rundat[0].eta = rundat[0].q = rundat[0].det = 1.0;

/***********************/
/* set fimgetae_() eps */
/***********************/
	for (i = 0; i < 3;  i++) eps[i] = epsmm;
	for (i = 3; i < 12; i++) eps[i] = epsmm/epsrad;

/********************************/
/* compute cross-run alignments */
/********************************/
	i = 0; modes[i++] = mo3d;
	if (stretch) modes[i++] = mo3d | (16 + 32 + 64);	/* enable xyz stretch */
	modes[i++] = 0;
	if (planar) for (i = 0; i < 8; i++) modes[i] &= ~(2 + 64);	/* turn off bits 2 and 64 */
	printf ("modes"); for (i = 0; i < 8; i++) printf (" %x", modes[i]); printf ("\n");
	load_filter_spline (filroo0, I0vol, img0);
	for (irun = 1; irun < nrun; irun++) {
		if (get_4dfp_dimoe (srcroot[irun], tmpdim, mmppixt, &orientt, &isbigt)) errr (srcroot[irun], program);
		status = tmpdim[3] != srcdim[3] || orientt != orients || isbigt != isbigs;
		for (k = 0; k < 3; k++) {
			status |= tmpdim[k] != srcdim[k];
			status |= fabs (mmppixt[k] - mmppixs[k]) > 1.e-4;
		}
		if (status) errd (program, srcroot[0], srcroot[irun]);
		strcpy (filroot, srcroot[irun]); if (use_geo) strcat (filroot, "_geom");
		load_filter_spline (filroot, I0vol, img1);
		volalign (img0, msk0, filroo0, img1, all1, filroot, modes,
			rundat[irun].t4, &rundat[irun].eta, &rundat[irun].q);
		determ12_ (rundat[irun].t4, &rundat[irun].det);
	}

/**************************/
/* allocate output memory */
/**************************/
	free (imgs);
	if (!(imgs = (float *) malloc (outsiz * sizeof (float))))		errm (program);
	if (!(imgt = (float *) calloc (outsiz * outdim[3], sizeof (float))))	errm (program);
	if (!(imgw = (float *) calloc (outsiz * outdim[3], sizeof (float))))	errm (program);

/***************/
/* open output */
/***************/
	sprintf (outfile, "%s_xenc.4dfp.img", srcroot[nrun]);
	printf ("Writing: %s\n", outfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);

/******************/
/* start dat file */
/******************/
	sprintf (datfile, "%s_xenc.4dfp.dat", srcroot[nrun]);
        if (!(fp_dat = fopen (datfile, "w"))) errw (program, datfile);
	printf ("Writing: %s\n", datfile);
/****************/
/* loop on runs */
/****************/
	for (irun = 0; irun < nrun; irun++) {
/*********************************************/
/* open matfile created by dwi_xalign3d_4dfp */
/*********************************************/
		sprintf (matfile, "%s_xenc.mat", srcroot[irun]);
		if (!(fp_mat = fopen (matfile, "r"))) errr (program, matfile);

/*******************/
/* loop on volumes */
/*******************/
		for (ivol = 0; ivol < srcdim[3]; ivol++) {
			fscanf (fp_mat, "%*s %*s %d", &j);
			if (j != ivol + 1) errr (program, matfile);
			for (j = 0; j < 4; j++) {
				fscanf (fp_mat, "%f %f %f %f",
					xenc_t4+0+j, xenc_t4+4+j, xenc_t4+8+j, xenc_t4+12+j);
			}
/*********************************/
/* compute transform composition */
/*********************************/
			matmul_ (xenc_t4, rundat[irun].t4, total_t4, &four);
			fprintf (fp_dat, "volume %2d total_t4\n", ivol + 1);
			for (j = 0; j < 4; j++) {
				fprintf (fp_dat, "%10.6f%10.6f%10.6f%10.4f\n",
					total_t4[0+j], total_t4[4+j], total_t4[8+j], total_t4[12+j]);
			}
			determ12_ (total_t4, &det);
			fprintf (fp_dat, "det=%f\tweight=%f\n", det, rundat[irun].weight);
			load_4dfp_frame (srcroot[irun], srcdim, ivol, isbigs, img1);
			spline3d (program, img1, srcdim);
			j = srcsiz;
			resample_ (img1,img1+1*j,img1+2*j,img1+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,all1,modes,
						     total_t4,imgs,outdim+0,outdim+1,outdim+2,mmppixo,centers);
			for (i = 0; i < outsiz; i++) if (imgs[i] != 0.0) {
				(imgt + ivol*outsiz)[i] += rundat[irun].weight * imgs[i] * det;
				(imgw + ivol*outsiz)[i] += rundat[irun].weight;
			}
		}	/* ivol */
		fclose (fp_mat);
		if (append) {
			if (ewrite (imgt, outsiz*outdim[3], control, fp_out)) errw (program, outfile);
			for (i = 0; i < outsiz*outdim[3]; i++) imgt[i] = 0.0;
		}
	}	/* irun */

	if (append) {
		outdim[3] *= nrun;
	} else {
/**************************************************/
/* average all volumes over runs and write result */
/**************************************************/
		for (i = 0; i < outsiz * outdim[3]; i++) {
			if (imgw[i]) imgt[i] /= imgw[i];
			if (nflag && imgt[i] < 0.0) {imgt[i] = 0.0; nzero++;}
		}
		if (ewrite (imgt, outsiz*outdim[3], control, fp_out)) errw (program, outfile);
	}
	if (fclose (fp_out)) errw (program, outfile);

/***********************/
/* outfile ifh and hdr */
/***********************/
	writeifhe (program, outfile, outdim, mmppixo, orients, control);
	sprintf (string, "ifh2hdr %s -r2000", outfile); system (string);

/*********************/
/* complete dat file */
/*********************/
	fprintf (fp_dat, "#      run       eta         q       det    weight\n");
	for (irun = 0; irun < nrun; irun++) {
		fprintf (fp_dat, "%10d%10.6f%10.6f%10.6f%10.6f\n", irun + 1,
		rundat[irun].eta, rundat[irun].q, rundat[irun].det, rundat[irun].weight);
	}
	fclose (fp_dat);

/***********/
/* recfile */
/***********/
	startrecle (outfile, argc, argv, rcsid, control);
	sprintf (recfile, "%s.rec", outfile);
	if (!(fp_rec = fopen (recfile, "a"))) errw (program, recfile);
	fprintf (fp_rec, "Cross-run alignment computed on basis of volume %d\n", I0vol + 1);
	if ((int) fsmooth) fprintf (fp_rec, "Run-to-run alignment pre-filter fhalf = %.4f (1/mm)\n", fsmooth);
	fprintf (fp_rec, "Cross-run alignment sampling interval = %.4f mm\n", del);
	fprintf (fp_rec, "Translation parameter search radius = %.4f mm\n", epsmm);
	fprintf (fp_rec, "Parameter search object radius = %.4f mm\n", epsrad);
	for (irun = 0; irun < nrun; irun++) {
		fprintf (fp_rec, "t4 run %d to run 1\n", irun + 1);
		for (j = 0; j < 4; j++) {
			fprintf (fp_rec, "%10.6f%10.6f%10.6f%10.4f\n",
			rundat[irun].t4[0+j], rundat[irun].t4[4+j], rundat[irun].t4[8+j], rundat[irun].t4[12+j]);
		}
	}
	if (nflag) fprintf (fp_rec, "%d output voxels zeroed to eliminate negative values\n", nzero);
	fprintf (fp_rec, "Zoom factors = {%.4f %.4f %.4f}\n", zoomfac[0], zoomfac[1], zoomfac[2]);
	fclose (fp_rec);
/*	catrec (datfile);	*/
	if (strlen (mskroot)) {
		sprintf (mskfile, "%s.4dfp.img", mskroot);
		catrec (mskfile);
	}
	for (irun = 0; irun < nrun; irun++) {
		sprintf (srcfile, "%s.4dfp.img", srcroot[irun]); catrec (srcfile);
	}
	endrec ();

/************/
/* clean up */
/************/
	free (imgs); free (img1); free (img0);
	free (msk0); free (all1);
	free (imgt); free (imgw);
	exit (status);
}

void filter3d (char *program, float *imgt, int *imgdim, float *mmppix, float f0, int wrap) {
	float		val, *imgp;
	int		k, margin, nxp, nyp, nzp;

	if (wrap) {
		nxp = imgdim[0];
		nyp = imgdim[1];
	} else {
		val = (0.5 + (2.0 * 0.1874 / (mmppix[0] * f0))); margin = val; nxp = npad_ (imgdim + 0, &margin);
		val = (0.5 + (2.0 * 0.1874 / (mmppix[1] * f0))); margin = val; nyp = npad_ (imgdim + 1, &margin);
	}
		val = (0.5 + (4.0 * 0.1874 / (mmppix[2] * f0))); margin = val; nzp = npad_ (imgdim + 2, &margin);
	printf ("image dimensions %d %d %d padded to %d %d %d \n", imgdim[0], imgdim[1], imgdim[2], nxp, nyp, nzp);
        if (!(imgp = (float *) calloc (nxp * nyp * nzp, sizeof (float)))) errm (program);
	imgpad_ (imgt, imgdim + 0, imgdim + 1, imgdim + 2, imgp, &nxp, &nyp, &nzp);
	gauss3d_  (imgp, &nxp, &nyp, &nzp, mmppix, &f0);
	imgdap_ (imgt, imgdim + 0, imgdim + 1, imgdim + 2, imgp, &nxp, &nyp, &nzp);
	free (imgp);
}

void spline3d (char *program, float *imgf, int *imgdim) {
	float		*imgp;
	int		k, nzt, four = 4, vdim;

	vdim = 1; for (k = 0; k < 3; k++) vdim *= imgdim[k];
	splinex_  (imgf, imgdim + 0, imgdim + 1, imgdim + 2, imgf + vdim*1);
	spliney_  (imgf, imgdim + 0, imgdim + 1, imgdim + 2, imgf + vdim*2);
	if (imgdim[2] < 48) {
		splineza_ (imgf,  imgdim + 0, imgdim + 1, imgdim + 2, imgf + vdim*3);
	} else {
		nzt = npad_ (imgdim + 2, &four);
		k = imgdim[0] * imgdim[1] * nzt;
		if (!(imgp = (float *) malloc (k * sizeof (float) * 2))) errm (program);
		imgpad_   (imgf,          imgdim + 0, imgdim + 1, imgdim + 2, imgp,     imgdim + 0, imgdim + 1, &nzt);
		splinezf_ (imgp,          imgdim + 0, imgdim + 1, &nzt,       imgp + k);
		imgdap_   (imgf + vdim*3, imgdim + 0, imgdim + 1, imgdim + 2, imgp + k, imgdim + 0, imgdim + 1, &nzt);
		free (imgp);
	}
}

void volalign (float *tar, short *tarm, char *tarstr, float *src, short *srcm, char *srcstr,
		long *modes, float *t4, float *peta, float *pq) {
	int		i, j, k, msrc, mtar;
	float		omega[12];

	printf ("aligning %s to %s\n", srcstr, tarstr);
	for (msrc = mtar = i = 0; i < srcsiz; i++) {
		if (tarm[i]) mtar++;
		if (srcm[i]) msrc++;
	} 
	printf ("target voxels=%d\tsource voxels=%d\n", mtar, msrc);

	j = srcsiz;
	i = 0; while (modes[i]) {
		fimgetae_ (tar,tar+1*j,tar+2*j,tar+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,tarm,modes+i,&del,eps,
			   src,src+1*j,src+2*j,src+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,srcm,t4);
		i++;
	} if (i) i--;
	t4_omega_(modes+i, t4, omega);
	etagen_  (tar,tar+1*j,tar+2*j,tar+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,tarm,modes+i,omega,&del,
		  src,src+1*j,src+2*j,src+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,srcm,peta,pq);
}

void load_filter_spline (char *srcroot, int ivol, float *imag) {
	load_4dfp_frame (srcroot, srcdim, ivol, isbigs, imag);
	if ((int) fsmooth) filter3d (program, imag, srcdim, mmppixs, fsmooth, wrap);
	spline3d (program, imag, srcdim);
}

