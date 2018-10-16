/*$Header: /data/petsun4/data1/src_solaris/cross_realign3d_4dfp/RCS/cross_realign3d_4dfp.c,v 1.23 2013/02/14 03:12:25 avi Exp $*/
/*$Log: cross_realign3d_4dfp.c,v $
 * Revision 1.23  2013/02/14  03:12:25  avi
 * remove limit on frames/run
 *
 * Revision 1.22  2010/12/23  06:37:41  avi
 * voxel size mismatch tolerance 1.e-5 -> 1.e-4
 *
 * Revision 1.21  2007/08/08  01:04:51  avi
 * startrec() -> startrece()
 *
 * Revision 1.20  2007/08/08  00:25:40  avi
 * masks now int instead of short (imag2mask -> img2lmask)
 * integer*2 -> integer*4 in all FORTRAN subroutines
 *
 * Revision 1.19  2007/04/16  04:31:22  avi
 * correct failure to initialize input stack read in 'make mode' (when .mat and _r3d_avg.4dfp files exist)
 *
 * Revision 1.18  2006/09/28  21:46:29  avi
 * Solaris 10
 *
 * Revision 1.17  2006/03/27  04:51:28  avi
 * correct omission of calling filter() after reading in first run _r3d
 * when first mat and _r3d files already exist
 *
 * Revision 1.16  2006/03/27  04:24:49  avi
 * endian invariant i/o
 *
 * Revision 1.15  2005/09/29  04:08:50  avi
 * option -c
 *
 * Revision 1.14  2005/08/28  00:02:06  avi
 * options -m and -b
 *
 * Revision 1.13  2004/12/13  21:11:17  rsachs
 * Replaced 'Get4dfpDimN' with 'get_4dfp_dimo'. Got rid of 'ifh_flag' & the -V option.
 *
 * Revision 1.12  2004/05/09  00:58:42  avi
 * output undefined voxels as 1.0e-37
 * override with -Z option
 *
 * Revision 1.11  2003/03/04  05:02:14  avi
 * eliminate separate cross_realign3d_4dfp.h include file
 * MAXR increase to 160
 * use modern utility subroutines (errr() errw() errm() getroot())
 *
 * Revision 1.10  2000/02/08  23:34:49  avi
 * make cross-run dimension match check work by including math.h
 *
 * Revision 1.9  1999/03/05  05:52:29  avi
 * list file fmt= ref= fields
 * process runs with non-functional frames distributed arbitrarily
 *
 * Revision 1.8  1999/01/22  04:53:29  avi
 * global_matstat
 *
 * Revision 1.7  1999/01/21  08:25:31  avi
 * self-realigned average image name r3d_avg
 * improved target and t4 info in rec files
 *
 * Revision 1.6  1999/01/21  06:59:53  avi
 * works saving self-aligned average
 *
 * Revision 1.5  1999/01/20  03:29:39  avi
 * Revision 1.4  1999/01/20  03:17:05  avi
 * further simplification (4 loops -> 2)
 * mat file written as soon as it is computed
 *
 * Revision 1.3  1999/01/19  05:48:38  avi
 * major cleaning prior to adding new functionality
 * sh shell eliminated
 *
 * Revision 1.2  1997/09/11  06:22:38  avi
 * norm_flag and optional intensity normalize disenable
 **/
/*_________________________________________________________________________
Module:		cross_realign3d_4dfp.c
Date:		April 11, 1996
Authors:	Avi Snyder, Tom Yang 
Description:	Read a group of 4dfp stacks, self-realign each stack,
		cross-realign each to reference (first named) run,
		and write each realigned stack as a xr3d.4dfp file.
_________________________________________________________________________*/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>				/* R_OK */
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>
#include <librms.h>

#define CRIT		0.2	/* img2lmask inclusion threshold */
#define FBLUR		0.06	/* half frequency in reciprocal mm for gauss3diz */
#define FHALF		0.02	/* half frequency in reciprocal mm for img2lmask */
#define MAXL		256	/* maximum string length */
#define MAXR		160	/* maximum number of runs */
#define MAXF		2048	/* expected maximum number of frames per fMRI run */

/*************/
/* mode word */
/*************/
#define ENABLE		1
#define ALIGN3D		2
#define STRETCH		4
#define GRADCOMP	8
#define CROSSMODAL	256
#define QUIET		512
#define WRAP		1024
#define ENDSLICE	2048

/*************/
/* externals */
/*************/
extern int	expandf (char *format, int nframe);			/* expandf.c */
extern void	alignmiz_ (int *nx, int *ny, int *nz, float *img1, float *img2, float *d2img2, int *mask,
			float *param, float *s4, float *mmppix, int *mode, float *err);			/* frmsfmri.f */ 
extern void	alignmrixyz_ (int *nx, int *ny, int *nz, float *imag, float *d2xi, float *d2yi, float *d2zi,
			float *imgr, int *mdef, float *param, float *s4, float *mmppix, int *mode);	/* alignmrixyz.f */
extern void	splinex_ (float *imag, int *nx, int *ny, int *nz, float *d2xi);				/* splinexyz.f */
extern void	spliney_ (float *imag, int *nx, int *ny, int *nz, float *d2yi);				/* splinexyz.f */
extern void	splinez_ (float *imag, int *nx, int *ny, int *nz, float *d2zi, float *scratch);		/* splinexyz.f */

void errd (char *program, char *file0, char *file1) {
	fprintf (stderr, "%s: %s %s dimension mismatch\n", program, file0, file1);
	exit (-1);
}

void write_min_max (float *imag, int n) {
	int	k;
	float	amax = -FLT_MAX;
	float	amin =  FLT_MAX;

	for (k = 0; k < n; k++) {
		if (imag[k] > amax) amax = imag[k];
		if (imag[k] < amin) amin = imag[k];
	}
	printf ("min=%10.4f      max=%10.4f\n", amin, amax);
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

/********************/
/* global variables */
/********************/
static float	*d2xi, *d2yi, *d2zi, *scratch;	/* computed second partials for spline interpolation */
static float	*img2;				/* unaligned frame */
static float	*imgf;				/* reference frame */
static float	*imgr;				/* aligned frame */
static float	*imgs;				/* scratch frame */
static float	*imga, *imge;			/* self-aligned functional frame stack averages */
static int	*mask;				/* region in which to evaulate and compute realignment */
static int	*mdef;				/* defined voxels returned by alignmrixyz */
static int	*imgn;				/* cumulative number of defined voxels returned by alignmrixyz */
static float	crit = CRIT;			/* img2lmask inclusion threshold */
static float	fblur = FBLUR;			/* half frequency in reciprocal mm for gauss3diz */
static float	fhalf = FHALF;			/* half frequency in reciprocal mm for img2lmask */

static char rcsid[] = "$Id: cross_realign3d_4dfp.c,v 1.23 2013/02/14 03:12:25 avi Exp $";
int main (int argc, char *argv[]) {
	int	spline (int threeD, int *imgdim, float *img2, float *d2xi, float *d2yi, float *d2zi);	/* below */
	void	filter (int threeD, float *imgs, int *imgdim, float *mmppix, float *pfblur);		/* below */
	void	WintraT (FILE *fp, int frame, float *totalT, float* tp);				/* below */
	void	WinterT (FILE *fp, char* imgroot, float *totalT, float* tp);				/* below */

	FILE		*imgfp, *mskfp, *avgfp, *matfp, *outfp;
	char		control = '\0';			/* controls endian state of image output */
	char		imgroot[MAXR][MAXL], mskroot[MAXR][MAXL];
	char		format[MAXR][MAXL];		/* anat/functional frame format */
	int		nfram[MAXR];			/* frames per run */
	int		ref_frame[MAXR];		/* frame numbers count from zero */
	int		matstat[MAXR];			/* set if mat file exists */
	int		orient[MAXR];			/* input orientation */
	int		isbig[MAXR], msbig[MAXR];	/* endian state of (4dfp) stacks and corresponding masks */
	int		asbig;				/* endian state of _r3d file */
	char		*bigfmt;			/* long string to hold expanded run format */
	char		lstfile[MAXL] = "", imgfile[MAXL], outfile[MAXL];
	char		matfile[MAXL], avgfile[MAXL], ifhfile[MAXL], mskfile[MAXL];

	float		*rp, *intraR, *interR;			/* alignment parameters */
	float		*intraS, *interS;			/* intensity scaling parameters */
	float		intraT[16], interT[16], totalT[16];	/* t4 matrices */
	float		*sp, *tp;				/* intensity scaling parameters */
	float		up[4] = {1., 0., 0., 0.};		/* unit intensity scaling parameters */

	float		voxdim[3], mmppix[3];
	int		imgdim[4], mskdim[4], avgdim[4];
	int		mode = ENABLE | ALIGN3D | ENDSLICE;	/* ordinary realignment mode */
	int		modef;					/* cross-modality realignment mode */
	int		ref_frame0 = 0;				/* frame numbers count from zero */
	int 		base_frame;				/* used as intraR pointer */
	int		nrun = 0, nframtot = 0, nf_anat0 = 0, nf_anat;
	int		xdim, ydim, zdim, volume_dim;

/*********/
/* flags */
/*********/
	int		e37_flag = 1;			/* output undefined voxels as 1.0e-37 instead of 0.0 */
	int		force_flag = 0;			/* force recomputing xr3d output files */
	int		norm_flag = 1;			/* enable intensity normalization */
	int		global_matstat = 1;
	int		global_mask = 0;
	int		debug = 0;

/***********/
/* utility */
/***********/
	float		err;
	int		c, i, j, k, m, ierr;
	int		four = 4;
	char		program[MAXL], *ptr, command[MAXL], *srgv[MAXL];

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
	for (i = 0; i < MAXR; i++) {
		format[i][0] = mskroot[i][0] = '\0';
		ref_frame[i] = 0;
	}

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'Z': e37_flag = 0;			break;
				case 'd': debug++;			break;
				case 'f': force_flag++;			break;
				case 'c': mode |= CROSSMODAL;		break;
				case 'g': mode |= GRADCOMP;		break;
				case 'p': mode &= ~ALIGN3D;		break;
				case 'q': mode |= QUIET;		break;
				case 's': mode |= STRETCH;		break;
				case 'w': mode |= WRAP; 		break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
				case 'l': strcpy (lstfile, ptr);	*ptr = '\0'; break;
				case 'm': getroot (ptr, mskroot[0]);
					  global_mask++;		*ptr = '\0'; break;
				case 'n': nf_anat0 = atoi (ptr);	*ptr = '\0'; break;
				case 'b': fblur = atof (ptr);		*ptr = '\0'; break;
				case 'r': ref_frame0 = atoi (ptr);	*ptr = '\0'; break;
				case 'v': norm_flag = atoi (ptr);	*ptr = '\0'; break;
			}
		} else {
			getroot (argv[i], imgroot[nrun++]);
		}	
	}

/*******************/
/* parse list file */
/*******************/
	if (strlen (lstfile)) {
		if (!(imgfp = fopen (lstfile, "r"))) errr (program, lstfile);
		while (fgets (command, MAXL, imgfp)) {
			if (nrun >= MAXR) {
				fprintf (stderr, "%s: maximum run count (%d) exceeded\n", program, MAXR);
				exit (-1);
			}
			if (ptr = strchr (command, '#'))  *ptr = '\0';
			if (!strlen (command)) continue;		/* skip blank lines */
			if (ptr = strchr (command, '\n')) *ptr = '\0';	/* strip terminal nl */
			i = m = 0; while (m < MAXL && i < MAXL) {
				while (!isgraph ((int) command[i]) && command[i]) i++;
				if (!command[i]) break;
				srgv[m++] = command + i;
				while (isgraph ((int) command[i])) i++;
				if (!command[i]) break;
				command[i++] = '\0';
			}
			getroot (srgv[0], imgroot[nrun]);
			for (i = 1; i < m; i++) {
				if (!strncmp (srgv[i], "fmt=",  4)) strcpy (format[nrun], srgv[i] + 4);
				if (!strncmp (srgv[i], "ref=",  4)) ref_frame[nrun] = atoi (srgv[i] + 4) - 1;
				if (!strncmp (srgv[i], "mask=", 5)) getroot (srgv[i] + 5, mskroot[nrun]);
			}
			if (debug) printf ("imgroot=%s mskroot=%s format=%s ref_frame=%d\n",
				imgroot[nrun], mskroot[nrun], format[nrun], ref_frame[nrun]);
			nrun++;
		}
		fclose (imgfp);
	}

	if (!nrun) {
		fprintf (stderr, "Usage:\t%s -l4dfp_list_file\n", program);
		fprintf (stderr, "   or:\t%s <run1_4dfp> <run2_4dfp> ...\n", program);
		fprintf (stderr, " e.g.:\t%s run1_4dfp run2_4dfp run3_4dfp\n", program);
		fprintf (stderr, "	%s -sqwv -lruns_4dfp.lst\n", program);
		fprintf (stderr, "	%s -pwqsf -n3 -lruns_4dfp.lst\n", program);
		fprintf (stderr, "	option\n");
		fprintf (stderr, "	-d	debug mode\n");
		fprintf (stderr, "	-@<b|l>\toutput big or little endian (default CPU endian)\n");
		fprintf (stderr, "	-f	force recomputing even if output files exist\n");
		fprintf (stderr, "	-g	enable linear intensity gradient compensation\n");
		fprintf (stderr, "	-c	use cross-modal registration always\n");
		fprintf (stderr, "	-l<str>	specify list file of 4dfp filenames\n");
		fprintf (stderr, "	-m<str>	specify 4dfp mask to be applied to all runs (default compute)\n");
		fprintf (stderr, "	-n<int>	specify number of pre-functional frames\n");
		fprintf (stderr, "	-b<flt>	specify pre-blur in reciprocal mm (default=%.2f)\n", FBLUR);
		fprintf (stderr, "	-p	2D (planar) realignment (default 3D)\n");
		fprintf (stderr, "	-q	minimize status reporting\n");
		fprintf (stderr, "	-r<int>	specify non-default reference frame\n");
		fprintf (stderr, "	-s	enable stretch\n");
		fprintf (stderr, "	-v[0|1]	disable/enable per frame intensity normalization (default enabled)\n");
		fprintf (stderr, "	-w	enable wrap addressing\n");
		fprintf (stderr, "	-Z	output undefined voxels as 0.0 (default 1.0e-37)\n");
		exit (1);
	}
	
/*****************************/
/* compute output endianness */
/*****************************/
	if (!control) control = (CPU_is_bigendian()) ? 'b' : 'l';

/**************************/
/* begin zeroth runs loop */
/**************************/
	for (i = 0;  i < nrun; i++) {
		sprintf (imgfile, "%s.4dfp.img", imgroot[i]);
		if (get_4dfp_dimoe (imgfile, imgdim, voxdim, orient + i, isbig + i)) errr (program, imgfile);
		if (!i) {
			xdim		= imgdim[0];
			ydim		= imgdim[1];
			zdim		= imgdim[2];
			mmppix[0]	= voxdim[0];
			mmppix[1]	= voxdim[1];
			mmppix[2]	= voxdim[2];
			volume_dim	= xdim * ydim * zdim;
		} else {
			ierr = (imgdim[0] != xdim || imgdim[1] != ydim || imgdim[2] != zdim);
			for (k = 0; k < 3; k++) ierr |= (fabs (mmppix[k] - voxdim[k]) > 1.e-4);
			ierr |= orient[i] - orient[0];
			if (ierr) errd (program, imgroot[0], imgroot[i]);
		}
		if (strlen (mskroot[i])) {
			sprintf (mskfile, "%s.4dfp.img", mskroot[i]);
			if (get_4dfp_dimoe (mskfile, mskdim, voxdim, &j, msbig + i)) errr (program, mskfile);
			ierr = (imgdim[0] != xdim || imgdim[1] != ydim || imgdim[2] != zdim);
			for (k = 0; k < 3; k++) ierr |= (fabs (mmppix[k] - voxdim[k]) > 1.e-4);
			ierr |= orient[i] - j;
			if (ierr) errd (program, mskroot[i], imgroot[i]);
		}
		nframtot += (nfram[i] = imgdim[3]);

		sprintf (matfile, "%s_xr3d.mat", imgroot[i]);
		sprintf (avgfile, "%s_r3d_avg.4dfp.img", imgroot[i]);
		matstat[i] = !access (matfile, R_OK) && !access (avgfile, R_OK);
		global_matstat &= matstat[i];
	}
/************************/
/* end zeroth runs loop */
/************************/
	if (!(intraR	= (float *) calloc (nframtot * 12, sizeof (float)))
	||  !(intraS	= (float *) calloc (nframtot * 4,  sizeof (float)))
	||  !(interR	= (float *) calloc (nrun * 12, sizeof (float)))
	||  !(interS	= (float *) calloc (nrun * 4,  sizeof (float)))
	||  !(imge	= (float *) calloc (nrun * volume_dim, sizeof (float)))
	||  !(scratch	= (float *) malloc (zdim * zdim * sizeof (float)))
	||  !(d2xi	= (float *) malloc (volume_dim * sizeof (float)))
	||  !(d2yi	= (float *) malloc (volume_dim * sizeof (float)))
	||  !(d2zi	= (float *) calloc (volume_dim,  sizeof (float)))
	||  !(imgf	= (float *) malloc (volume_dim * sizeof (float)))
	||  !(img2	= (float *) malloc (volume_dim * sizeof (float)))
	||  !(imgr	= (float *) malloc (volume_dim * sizeof (float)))
	||  !(imgs	= (float *) malloc (volume_dim * sizeof (float)))
	||  !(mdef	= (int *)   malloc (volume_dim * sizeof (int)))
	||  !(imgn	= (int *)   malloc (volume_dim * sizeof (int)))
	||  !(mask	= (int *)   malloc (volume_dim * sizeof (int)))) errm (program);

	param2warp_ (&mode, interR, interT);	/* reference run cross t4 */
	interS[0] = 1.0;			/* reference run scaled to unity */
	if (!(bigfmt = (char *) calloc (MAXF, sizeof (char)))) errm (program);

/************************/
/* begin main runs loop */
/************************/
	base_frame = 0;
	imga = imge;
	for (i = 0; i < nrun; i++) {
		sprintf (matfile, "%s_xr3d.mat", imgroot[i]);
		sprintf (imgfile, "%s.4dfp.img", imgroot[i]);
		fprintf (stdout, "Reading: %s\n", imgfile);
		if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
		if (global_matstat || (i && matstat[i])) goto RESAMPLE;
		if (!i && matstat[i]) {
			sprintf (avgfile, "%s_r3d_avg.4dfp.img", imgroot[i]);
			if (!get_4dfp_dimoe (avgfile, avgdim, voxdim, &j, &asbig)
			&& (avgfp = fopen (avgfile, "rb"))
			&& !gread ((char *) imga, sizeof (float), volume_dim, avgfp, asbig)
			&& !fclose (avgfp)) {
				filter (mode, imga, imgdim, mmppix, &fblur);
				goto RESAMPLE;
			}
		}

/*******************************/
/* prepare to self-align run i */
/*******************************/
		if (!strlen (format[i])) sprintf (format[i], "%dx%d+", nf_anat0, nfram[i] - nf_anat0);
		strcpy (bigfmt, format[i]);
		if (!(bigfmt = (char *) realloc (bigfmt, nfram[i]*sizeof(char) + 4))) errm (program);
		if (expandf (bigfmt, nfram[i] + 4)) exit (-1);
		nf_anat = k = 0;
		while (j = bigfmt[k++]) if (j == 'x') nf_anat++;
		if (!ref_frame[i]) ref_frame[i] = (ref_frame0) ? ref_frame0 - 1: (nfram[i] + nf_anat0) / 2;
		k = ref_frame[i];
		if (k < 0 || k >= nfram[i] || bigfmt[k] == 'x') {
			fprintf (stderr, "%s: illegal %s reference frame (%d)\n", program, imgroot[i], k + 1);
			exit (-1);
		}
		if (fseek (imgfp, (long) ref_frame[i] * volume_dim * sizeof (float), SEEK_SET)
		||  gread ((char *) imgf, sizeof (float), volume_dim, imgfp, isbig[i])) errr (program, imgfile);

/**********************/
/* prepare run i mask */
/**********************/
		k = -1; if (strlen (mskroot[i])) k = i; if (global_mask) k = 0;
		if (k >= 0) {
			sprintf (mskfile, "%s.4dfp.img", mskroot[k]);
			printf ("Reading: %s\n", mskfile);
			if (!(mskfp = fopen (mskfile, "rb"))
			|| gread ((char *) imgs, sizeof (float), volume_dim, mskfp, msbig[i])
			|| fclose (mskfp)) errr (program, mskfile);
			for (j = 0; j < volume_dim; j++) mask[j] = imgs[j];
		} else {
			img2lmask (&xdim, &ydim, &zdim, imgf, mask, mmppix, &fhalf, &crit);
		}

/**************************/
/* filter reference frame */
/**************************/
		filter (mode, imgf, imgdim, mmppix, &fblur);

/***************************/
/* align functional frames */ 
/***************************/
		printf ("Computing: %s to frame %d intra-run alignment\n", imgroot[i], ref_frame[i] + 1);
		for (k = 0; k < volume_dim; k++) imga[k] = imgn[k] = 0;
		for (j = 0; j < nfram[i]; j++) {
			if (bigfmt[j] == 'x') continue;
			printf ("frame %4d ", j + 1); fflush (stdout);
			rp	= intraR + (base_frame + j) * 12;
			sp	= intraS + (base_frame + j) * 4;
			if (fseek (imgfp, (long) j * volume_dim * sizeof (float), SEEK_SET)
			||  gread ((char *) img2, sizeof (float), volume_dim, imgfp, isbig[i])) errr (program, imgfile);

			for (k = 0; k < volume_dim; k++) imgs[k] = img2[k];
			filter (mode, imgs, imgdim, mmppix, &fblur);
			if (mode & ALIGN3D) splinez_ (imgs, &xdim, &ydim, &zdim, d2zi, scratch);
			alignmiz_ (&xdim, &ydim, &zdim, imgf, imgs, d2zi, mask, rp, sp, mmppix, &mode, &err);
			for (k = 0; k < 12; k++) if (isnan (rp[k])) goto ERRC;

			if (spline (mode, imgdim, img2, d2xi, d2yi, d2zi)) errm (program);
			for (k = 0; k < volume_dim; k++) mdef[k] = 1;
			alignmrixyz_ (&xdim, &ydim, &zdim, img2, d2xi, d2yi, d2zi, imgr, mdef, rp, sp, mmppix, &mode);
/***********************************************/
/* add realigned frame to functional frame sum */
/***********************************************/
			for (k = 0; k < volume_dim; k++) if (mdef[k]) {
				imga[k] += imgr[k];
				imgn[k]++;
			}
		}
		for (k = 0; k < volume_dim; k++) if (imgn[k]) imga[k] /= imgn[k];

/******************************************/
/* save self-aligned functional frame sum */
/******************************************/
		sprintf (outfile, "%s_r3d_avg.4dfp.img", imgroot[i]);
		fprintf (stdout, "Writing: %s\n", outfile);
		if (!(outfp = fopen (outfile, "wb"))
		|| gwrite ((char *) imga, sizeof (float), volume_dim, outfp, control)
		|| fclose (outfp)) errw (program, outfile);
		imgdim[3] = 1;
		if (writeifhe (program, outfile, imgdim, voxdim, orient[i], control)) errw (program, outfile);
		startrece (outfile, argc, argv, rcsid, control);
		
		sprintf (command, "Paradigm format: %s\n", format[i]); printrec (command);
		sprintf (command, "Intra-run reference frame = %d\n", ref_frame[i] + 1); printrec (command);
		sprintf (command, "Gaussian pre-blur fhalf = %.4f (1/mm)\n", fblur); printrec (command);
		k = -1; if (strlen (mskroot[i])) k = i; if (global_mask) k = 0;
		if (k >= 0) {
			sprintf (command, "Externally supplied mask: %s\n", mskroot[k]); printrec (command);
		}
		catrec (imgfile);
		endrec ();

/************************************/
/* self-align non-functional frames */
/************************************/
		modef = mode | CROSSMODAL;
		printf ("Computing: %s to functional frame avg intra-run alignment\n", imgroot[i]);
		for (j = 0; j < nfram[i]; j++) { 
			if (bigfmt[j] != 'x') continue;
			printf ("frame %4d ", j + 1); fflush (stdout);
			rp	= intraR + (base_frame + j) * 12;
			sp	= intraS + (base_frame + j) * 4;
			if (fseek (imgfp, (long) j * volume_dim * sizeof (float), SEEK_SET)
			||  gread ((char *) imgs, sizeof (float), volume_dim, imgfp, isbig[i])) errr (program, imgfile);

			if (0) filter (mode, imgs, imgdim, mmppix, &fblur);	/* works better without */
			if (mode & ALIGN3D) splinez_ (imgs, &xdim, &ydim, &zdim, d2zi, scratch);
			alignmiz_ (&xdim, &ydim, &zdim, imga, imgs, d2zi, mask, rp, sp, mmppix, &modef, &err);
			for (k = 0; k < 12; k++) if (isnan (rp[k])) goto ERRC;
		}

/***********************************/
/* prepare for cross-run alignment */
/***********************************/
		if (strlen (mskroot[0])) {
			sprintf (mskfile, "%s.4dfp.img", mskroot[0]);
			printf ("Reading: %s\n", mskfile);
			if (!(mskfp = fopen (mskfile, "rb"))
			|| gread ((char *) imgs, sizeof (float), volume_dim, mskfp, msbig[0])
			|| fclose (mskfp)) errr (program, mskfile);
			for (j = 0; j < volume_dim; j++) mask[j] = imgs[j];
		} else {
			img2lmask (&xdim, &ydim, &zdim, imga, mask, mmppix, &fhalf, &crit);
		}
/************************************************************/
/* prepare functional frame average for cross-run alignment */
/************************************************************/
		filter (mode, imga, imgdim, mmppix, &fblur);

/*******************************/
/* compute cross run alignment */
/*******************************/
		if (i) {
			rp	= interR + i * 12;
			sp	= interS + i * 4;
			printf ("Aligning: %s to %s\n", imgroot[i], imgroot[0]);
			if (mode & ALIGN3D) splinez_ (imga, &xdim, &ydim, &zdim, d2zi, scratch);
			alignmiz_ (&xdim, &ydim, &zdim, imge, imga, d2zi, mask, rp, sp, mmppix, &mode, &err);
			for (k = 0; k < 12; k++) if (isnan (rp[k])) goto ERRC;
			param2warp_ (&mode, rp, interT);
			WinterT (stdout, imgroot[i], interT, sp);
/**********************/
/* compose transforms */
/**********************/
			for (j = 0; j < nfram[i]; j++) {
				rp = intraR + (base_frame + j) * 12;
				sp = intraS + (base_frame + j) * 4;
				param2warp_ (&mode, rp, intraT);
				matmul_ (intraT, interT, totalT, &four);
				warp2param_ (&mode, totalT, rp);
				sp[0] *= interS[4 * i];
				for (k = 1; k < 4; k++) sp[k] += interS[4 * i + k];
			}
		}
/*****************/
/* save mat file */
/*****************/
		if (!(matfp = fopen (matfile, "w"))) errw (program, matfile);
		fprintf (stdout, "Writing: %s\n", matfile);
		for (j = 0; j < nfram[i]; j++) {
			rp = intraR + (base_frame + j) * 12;
			sp = intraS + (base_frame + j) * 4;
			param2warp_ (&mode, rp, totalT);
			WintraT (matfp, j, totalT, sp);
		}
		fclose (matfp);

RESAMPLE:	sprintf (outfile, "%s_xr3d.4dfp.img", imgroot[i]);
		sprintf (ifhfile, "%s_xr3d.4dfp.ifh", imgroot[i]);
		if (!force_flag && matstat[i] && !access (outfile, R_OK) && !access (ifhfile, R_OK)) goto LOOP;
/*****************/
/* read mat file */
/*****************/
		if (!(matfp = fopen (matfile, "r"))) errr (program, matfile);
		fprintf (stdout, "Reading: %s\n", matfile);
		for (j = 0; j < nfram[i]; j++) {
			rp	= intraR + (base_frame + j) * 12;
			sp	= intraS + (base_frame + j) * 4;
			fscanf (matfp, "%*s %*s %*s");
			for (k = 0; k < 4; k++) fscanf (matfp, "%f %f %f %f",
				totalT + 0 + k, totalT + 4 + k, totalT + 8 + k, totalT + 12 + k);
			warp2param_ (&mode, totalT, rp);
			fscanf (matfp, "%*s %*s %*s");
			fscanf (matfp, "%f %f %f %f", sp + 0, sp + 1, sp + 2, sp + 3);
		}
		fclose (matfp);

		fprintf (stdout, "Writing: %s\n", outfile);
		if (!(outfp = fopen (outfile, "w+b"))) errw (program, outfile);
		for (j = 0; j < nfram[i]; j++) {
			if (gwrite ((char *) imgr, sizeof (float), volume_dim, outfp, control)) errw (program, outfile);
		}
		printf ("Resampling: %s\n", imgroot[i]);
		for (j = 0; j < nfram[i]; j++) {
			printf ("frame %4d ", j + 1); fflush (stdout);
			rp	= intraR + (base_frame + j) * 12;
			sp	= intraS + (base_frame + j) * 4;
			param2warp_ (&mode, rp, totalT);
			tp = (norm_flag) ? sp : up;	/* unit scaling factors if norm disabled */
			if (debug) WintraT (stdout, j, totalT, tp);

			if (fseek (imgfp, (long) j * volume_dim * sizeof (float), SEEK_SET)
			||  gread ((char *) img2, sizeof (float), volume_dim, imgfp, isbig[i])) errr (program, imgfile);
			if (0) write_min_max (img2, volume_dim);
			if (spline (mode, imgdim, img2, d2xi, d2yi, d2zi)) errm (program);
			for (k = 0; k < volume_dim; k++) mdef[k] = 1;
			alignmrixyz_ (&xdim, &ydim, &zdim, img2, d2xi, d2yi, d2zi, imgr, mdef, rp, tp, mmppix, &mode);
			if (!e37_flag) for (k = 0; k < volume_dim; k++) if (!mdef[k]) imgr[k] = 0.0;
			if (fseek (outfp, (long) volume_dim * j * sizeof (float), SEEK_SET)
			|| gwrite ((char *) imgr, sizeof (float), volume_dim, outfp, control)) errw (program, outfile);
		}
		fclose (imgfp);
		fclose (outfp);

		imgdim[3] = nfram[i];
		if (writeifhe (program, outfile, imgdim, voxdim, orient[i], control)) errw (program, outfile);
		startrece (outfile, argc, argv, rcsid, control);
		if (matstat[i]) {
			sprintf (command, "Cross-realigned xr3d stack resampled using existing mat file\n");
			printrec (command);
		} else {
			sprintf (command, "Cross-run reference = %s\n", imgroot[0]); printrec (command);
			sprintf (outfile, "%s+", imgroot[i]);
			if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
			WinterT (outfp, imgroot[i], interT, interS + i * 4); fclose (outfp);
			catrec (outfile); remove (outfile);
		}
		catrec (imgfile);
		endrec ();

LOOP:		base_frame	+= nfram[i];
		imga		+= volume_dim;
	}
/***********************/
/* end main runs loops */
/***********************/

	free (interR); free (intraR); free (interS); free (intraS);
	free (scratch); free (d2xi); free (d2yi); free (d2zi);
	free (imgs); free (img2); free (imgr); free (imgf); free (imge);
	free (mask); free (imgn); free (mdef);
	free (bigfmt);

	exit (0);
ERRC:	fprintf (stderr, "%s: convergence failed\n", program);		exit (-1);
}

void filter (int mode, float *imgs, int *imgdim, float *mmppix, float *pfblur) {
	extern void	gauss3diz (float *imag, int *pnx, int *pny, int *pnz, float *mmppix, float *pfhalf);	 /* gauss3diz.c */
	extern void	gauss2diz (float *imag, int *pnx, int *pny, int *pnz, float *mmppix, float *pfhalf);	 /* gauss2diz.c */

	if (mode & ALIGN3D) {
		gauss3diz (imgs, imgdim + 0, imgdim + 1, imgdim + 2, mmppix, pfblur);
	} else {
		gauss2diz (imgs, imgdim + 0, imgdim + 1, imgdim + 2, mmppix, pfblur);
	}
}

int spline (int mode, int *imgdim, float *img2, float *d2xi, float *d2yi, float *d2zi) {
	float		*scratch;
	int		k;

	splinex_ (img2, imgdim + 0, imgdim + 1, imgdim + 2, d2xi);
	spliney_ (img2, imgdim + 0, imgdim + 1, imgdim + 2, d2yi);
	if (mode & ALIGN3D) {
		if (!(scratch = (float *) malloc (imgdim[2] * imgdim[2] * sizeof (float)))) return -1;
		splinez_ (img2, imgdim + 0, imgdim + 1, imgdim + 2, d2zi, scratch);
		free (scratch);
	} else {
		for (k = 0; k < imgdim[0] * imgdim[1] * imgdim[2]; k++) d2zi[k] = 0.0;
	}
	return 0;
}

void WintraT (FILE *fp, int frame, float *totalT, float* tp) {
	int		k;

	fprintf (fp, "t4 frame %d\n", frame + 1);
	for (k = 0; k < 4; k++) fprintf (fp, "%10.6f%10.6f%10.6f%10.4f\n",
		totalT[0 + k], totalT[4 + k], totalT[8 + k], totalT[12 + k]);
	fprintf (fp, "s4 frame %d\n", frame + 1);
	fprintf (fp, "%10.6f%10.6f%10.6f%10.6f\n", tp[0], tp[1], tp[2], tp[3]);
}

void WinterT (FILE *fp, char* imgroot, float *totalT, float* tp) {
	int		k;

	fprintf (fp, "t4 run %s\n", imgroot);
	for (k = 0; k < 4; k++) fprintf (fp, "%10.6f%10.6f%10.6f%10.4f\n",
		totalT[0 + k], totalT[4 + k], totalT[8 + k], totalT[12 + k]);
	fprintf (fp, "s4 run %s\n", imgroot);
	fprintf (fp, "%10.6f%10.6f%10.6f%10.6f\n", tp[0], tp[1], tp[2], tp[3]);
}
