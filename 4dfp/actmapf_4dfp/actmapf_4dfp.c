/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/actmapf_4dfp.c,v 1.33 2013/01/03 07:14:30 avi Exp $*/
/*$Log: actmapf_4dfp.c,v $
 * Revision 1.33  2013/01/03  07:14:30  avi
 * correct tmpfp fopen error
 *
 * Revision 1.32  2012/08/04  01:22:08  avi
 * trap tmpfile write errors
 *
 * Revision 1.31  2012/03/06  23:10:49  avi
 * trap nnez < 1
 *
 * Revision 1.30  2006/09/24  00:25:39  avi
 * correct control initialization
 *
 * Revision 1.29  2006/09/23  23:03:44  avi
 * Solaris 10
 *
 * Revision 1.28  2006/09/21  19:52:02  avi
 * eliminate istart and istop variables
 * call #includes now use <> convention
 *
 * Revision 1.27  2005/12/02  06:35:02  avi
 * correct usage
 *
 * Revision 1.26  2005/09/05  01:07:29  justinv
 * MAXF increased to 16384
 *
 * Revision 1.25  2005/07/22  05:09:45  avi
 * accept conc file input
 * save mmppix and center in output map
 * MAXF 1024->4096
 *
 * Revision 1.24  2005/01/21  18:28:29  avi
 * allow comments in input profiles
 *
 * Revision 1.23  2004/05/23  04:05:30  avi
 * remove -b (baseline correct) option
 * add -R (relative modulation) option
 * run ifh2hdr
 *
 * Revision 1.22  2004/05/16  01:52:33  avi
 * -w option
 *
 * Revision 1.21  2000/09/04  00:33:20  avi
 * minor fix in input filename extension stripping
 *
 * Revision 1.20  1999/03/11  09:18:09  avi
 * compute sin images with correct (negative) sine
 * newer usage and utility programs
 *
 * Revision 1.19  1999/01/04  08:21:20  avi
 * remove #include <mri/mri.h>
 *
 * Revision 1.18  1999/01/04  08:14:45  avi
 * correct core dumps on ERRR and ERRW
 *
 * Revision 1.17  1998/11/16  20:54:17  avi
 * MAXF 520 -> 1024
 *
 * Revision 1.16  1998/10/09  03:49:51  avi
 * read one frame at a time
 *
 * Revision 1.15  1998/10/08  23:51:08  avi
 * MAXF -> 520
 * new rec calls
 *
 * Revision 1.14  1998/10/08  23:11:06  avi
 * before increasing MAXF
 *
 * Revision 1.13  1997/12/05  04:49:07  avi
 * compute cos and sin weights
 *
 * Revision 1.12  1997/05/23  02:37:50  yang
 * new rec macros
 *
 * Revision 1.11  1997/04/28  21:00:55  yang
 * Working Solaris version.
 *
 * Revision 1.10  1996/07/12  07:17:31  avi
 * rec file modifications
 *
 * Revision 1.9  1996/06/22  21:24:31  avi
 * Redo switch processing
 *
 * Revision 1.8  1996/06/18  03:49:54  avi
 * compute weights variance with respect to zero unless -z switch is set
 *
 * Revision 1.7  1996/05/22  06:22:08  avi
 * -a switch and new variable trailer
 *
 * Revision 1.6  1996/05/20  07:51:30  avi
 * exit if expandf returns nonzero
 *
 * Revision 1.5  1996/05/20  07:38:52  avi
 * Revision 1.4  1996/05/20  07:35:50  avi
 * Revision 1.3  1996/05/20  07:18:36  avi
 * cscale (activation map intensity scaling factor)
 *
 * Revision 1.2  1996/05/19  06:49:18  avi
 * correct frames vs. npts check
 *
 * Revision 1.1  1996/05/19  03:09:16  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <unistd.h>			/* getpid () */
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>
#include <conc.h>

#define MAXL		256
#define MAXF		16384		/* maximum npts coded in format */
#define MAX(a,b)	(a>b? a:b)
#define MIN(a,b)	(a<b? a:b)

int	expandf (char *string, int len);				/* expandf.c */
float	fimg_mode (float* fimg, int nval);				/* fimg_mode.c */

/********************/
/* global variables */
/********************/
static char rcsid[] = "$Id: actmapf_4dfp.c,v 1.33 2013/01/03 07:14:30 avi Exp $";

void usage (char *program) {
	fprintf (stderr, "Usage:\t%s <format> <4dfp|conc input>\n", program);
	fprintf (stderr, " e.g.,\t%s -zu \"3x3(11+4x15-)\" b1_rmsp_dbnd_xr3d_norm\n", program);
	fprintf (stderr, " e.g.,\t%s -aanatomy -c10 -u \"+\" ball_dbnd_xr3d.conc\n", program);
	fprintf (stderr, " e.g.,\t%s -zu \"4x124+\" b1_rmsp_dbnd_xr3d -wweights.txt\n", program);
	fprintf (stderr, "\toption\n");
	fprintf (stderr, "\t-a<str>\tspecify 4dfp output root trailer (default = \"actmap\")\n");
	fprintf (stderr, "\t-c<flt>\tscale output by specified factor\n");
	fprintf (stderr, "\t-u\tscale weights to unit variance\n");
	fprintf (stderr, "\t-z\tadjust weights to zero sum\n");
	fprintf (stderr, "\t-R\tcompute relative modulation (default absolute)\n");
	fprintf (stderr, "\t-w<weight file>\tread (text) weights from specified filename\n");
	fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default input endian)\n");
	fprintf (stderr, "N.B.:\tconc files must have extension \"conc\"\n");
	fprintf (stderr, "N.B.:\twhen using weight files 'x' frames in format are not counted\n");
	fprintf (stderr, "N.B.:\trelative modulation images are zeroed where mean intensity < 0.5*whole_image_mode\n");
	exit (1);
}

int main (int argc, char *argv[]) {
/**********************/
/* filename variables */
/**********************/
	CONC_BLOCK	conc_block;			/* conc i/o control block */
	FILE		*tmpfp, *imgfp, *weifp;
	char		imgroot[MAXL], imgfile[MAXL];	/* input 4dfp stack filename */
	char		outfile[MAXL], tmpfile[MAXL], weifile[MAXL] = "";
	char		trailer[MAXL] = "actmap";	/* appended to outfile name */

/**********************/
/* image dimensioning */
/**********************/
	IFH		ifh;
	float		*imgt;			/* one volume */
	float		*imgs;			/* weighted sum */
	float		*imga;			/* simple sum */
	int		index, imgdim[4], vdim;	/* image dimensions */
	int		isbig, osbig;
	char		control = '\0';

/*************************/
/* timeseries processing */
/*************************/
	char		*str, format[MAXF];
	float		theta, twopi;
	float		weight[MAXF];
	float		fmin, fmax, fmode;
	float 		cscale = 1.0;
	float		weight_sd, wt_mean, wt_var;
	int		npts, npos, nneg, nzer, nnez, nsin;

/***********/
/* utility */
/***********/
	char		command[MAXL], *ptr, program[MAXL], *srgv[MAXL];
	int		c, i, j, k, m;

/*********/
/* flags */
/*********/
	int		conc_flag = 0;
	int		unit_var = 0;
	int		zero_mean = 0;
	int		scale_rel = 0;		/* normalize weighted sum by voxelwise mean intensity */
	int		status = 0;
	int		debug = 0;

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
				case 'd': debug++;		break;
				case 'z': zero_mean++;		break;
				case 'u': unit_var++;		break;
				case 'R': scale_rel++;		break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
				case 'c': cscale = atof (ptr);		*ptr = '\0'; break;
				case 'w': strcpy (weifile, ptr);	*ptr = '\0'; break;
				case 'a': strcpy (trailer, ptr);	*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	strcpy (format,  argv[i]);	k++; break;
			case 1:	getroot (argv[i], imgroot);
				conc_flag = (strstr (argv[i], ".conc") == argv[i] + strlen (imgroot));
								k++; break;
		}	
	}
	if (k < 2) usage (program);

/******************************/
/* open temp process log file */
/******************************/
	sprintf (tmpfile, "temp%d", getpid ());
	if (!(tmpfp = fopen (tmpfile, "w"))) errw (program, tmpfile);

/*******************************************/
/* execute preliminary timeseries analysis */
/*******************************************/
	if (k = expandf (format, MAXF)) exit (k);
	printf ("%s\n", format);
	npts = strlen (format);
	printf ("%s: time series defined for %d frames\n", program, npts);
	npos = nneg = nzer = nsin = 0;
	twopi = 8.*atan (1.);
	for (k = 0; k < npts; k++) {
		switch (format[k]) {
			case 'x': weight[k] =  0.; nzer++; break;
			case '+': weight[k] =  1.; npos++; break;
			case '-': weight[k] = -1.; nneg++; break;
			case 'C': case 'c': case 'S': case 's': {
				str = strchr (format + k, '~');
				j = str - format - k + 1;
				for (i = 0; i < j; i++) {
					theta = twopi * (float) i / (float) j;
					switch (format[k]) {
						case 'C': case 'c': weight[k + i] =  cos (theta); break;
						case 'S': case 's': weight[k + i] = -sin (theta); break;
					}
				}
				nsin += j;
				k += j - 1;
				break;
			}
		}
	}
	printf ("%s: number positive=%d  negative=%d  sinusoidal=%d  zero=%d\n",
		program, npos, nneg, nsin, nzer);
	nnez = npts - nzer;
	fprintf (tmpfp, "Timepoint counts: positive=%d  negative=%d  sinusoidal=%d  zero=%d\n",
		npos, nneg, nsin, nzer);
	fprintf (tmpfp, "%s\n", format);
	if (nnez < 1) {
		fprintf (stderr, "%s: invalid format\n", program);
		exit (-1);
	}

	if (strlen (weifile)) {
/************************************/
/* read externally supplied weights */
/************************************/
		printf ("Reading: %s\n", weifile);
		if (!(weifp =  fopen (weifile, "r"))) errr (program, weifile);
		k = 0;	while (fgets (command, MAXF, weifp)) {
			if (!(m = split (command, srgv, MAXL))) continue;
			if (command[0] == '#') {
				printf ("%s", command);
			} else {
				weight[k++] = atof (srgv[0]);
			}
		}
		fclose (weifp);
		if (k < npts) {
			fprintf (stderr, "%s: %s lines (%d) less than format npts (%d)\n",
			program, weifile, k, npts);
			exit (-1);
		}
	}

/********************************************/
/* compute mean and variance of the weights */
/********************************************/
	for (k = 0; k < 2; k++) {
		wt_mean = wt_var = 0.;
		for (i = 0; i < npts; i++) if (format[i] != 'x') {
			wt_mean += weight[i];
			wt_var  += weight[i] * weight[i];
		}
		wt_mean /= (float) nnez;
		if (zero_mean) wt_var -= (float) nnez * wt_mean * wt_mean;
		wt_var /= (float) nnez; 
		weight_sd = sqrt (wt_var);

		if (zero_mean) {
			for (i = 0; i < npts; i++) if (format[i] != 'x') weight[i] -= wt_mean;
			zero_mean = 0;
		}

		if (unit_var) {
			for (i = 0; i < npts; i++) if (format[i] != 'x') weight[i] /= weight_sd;
			unit_var = 0;
		}
	}
	if (debug) {
		for (k = 0; k < npts; k++) {
			printf ("%4d %c %10.6f\n", k + 1, format[k], weight[k]);
		}
	}
	printf ("%s: weight_mean=%10.6f  weight_sd=%10.6f\n", program, wt_mean, weight_sd);
	fprintf (tmpfp, "Time series weights:  mean=%10.6f  sd=%10.6f\n", wt_mean, weight_sd);

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
		if (Getifh (imgfile, &ifh)) errr (program, imgfile);
		for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];
		if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
		printf ("Reading: %s\n", imgfile);
		isbig = strcmp (ifh.imagedata_byte_order, "littleendian");
	}
	if (!control) control = (isbig) ? 'b' : 'l';
	vdim = imgdim[0] * imgdim[1] * imgdim[2];
	if (imgdim[3] < npts) {
		fprintf (stderr, "%s: more defined npts (%d) than frames (%d)\n", program, npts, imgdim[3]);
		exit (-1);
	}
	imgs =	(float *) calloc (vdim, sizeof (float));
	imga =	(float *) calloc (vdim, sizeof (float));
	imgt =	(float *) calloc (vdim, sizeof (float));
	if (!imgs || !imgt || !imga) errm (program);

	printf ("computing weighted sum image scaled by %.4f\n", cscale);
	fprintf (tmpfp, "Weighted sum image scaled by %.4f\n", cscale);
	for (i = 0; i < npts; i++) {
		if (conc_flag) {
			conc_read_vol (&conc_block, imgt);
		} else {
			if (eread (imgt, vdim, isbig, imgfp)) errr (program, imgfile);
		}
		if (format[i] == 'x') continue;
		for (j = 0; j < vdim; j++) {
			imga[j] += imgt[j];
			imgs[j] += imgt[j]*weight[i];
		}
	}
	if (conc_flag) {
		conc_free (&conc_block);
	} else {
		if (fclose (imgfp)) errr (program, imgfile);
	}

	fmode = fimg_mode (imga, vdim);
/******************************/
/* compute final output image */
/******************************/
	fmin = FLT_MAX; fmax = -fmin;
	for (j = 0; j < vdim; j++) {
		if (scale_rel) {
			if (imga[j] > 0.5*fmode) {
				imgs[j] *= cscale / imga[j];
			} else {
				imgs[j] = 0.;
			}
		} else {
			imgs[j] *= cscale / (float) nnez;
		}
		fmin = MIN (fmin, imgs[j]);
		fmax = MAX (fmax, imgs[j]);
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
	if (!(imgfp = fopen (outfile, "wb")) || ewrite (imgs, vdim, control, imgfp)
	|| fclose (imgfp)) errw (program, outfile);

/*******/
/* ifh */
/*******/
	imgdim[3] = 1;
	if (conc_flag) {
		writeifhmce (program, outfile, imgdim,
			conc_block.voxdim, conc_block.orient, conc_block.mmppix, conc_block.center, control);
	} else {
		writeifhmce (program, outfile, imgdim,
			ifh.scaling_factor, ifh.orientation, ifh.mmppix, ifh.center, control);
	}

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr -r%dto%d %s", (int) (fmin-0.5), (int) (fmax+0.5), outfile);
	printf ("%s\n", command); status |= system (command);

/*******/
/* rec */
/*******/
	startrecle (outfile, argc, argv, rcsid, control);
	if (fclose (tmpfp)) errw (program, tmpfile);
	catrec (tmpfile); remove (tmpfile);
	catrec (imgfile);
	endrec ();

	free (imgt); free (imgs); free (imga);
	exit (status);
}
