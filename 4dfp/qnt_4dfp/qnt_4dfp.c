/*$Header: /data/petsun4/data1/src_solaris/qnt_4dfp/RCS/qnt_4dfp.c,v 1.19 2010/05/07 03:12:23 avi Exp $*/
/*$Log: qnt_4dfp.c,v $
 * Revision 1.19  2010/05/07  03:12:23  avi
 * option -A
 *
 * Revision 1.18  2008/02/07  04:49:46  avi
 * enable command buffer to hold long format strings
 *
 * Revision 1.17  2007/08/14  02:25:09  avi
 * format length now data determined
 *
 * Revision 1.16  2007/08/09  03:37:05  avi
 * linux compliant; use isnormal()
 *
 * Revision 1.15  2007/03/17  01:03:50  avi
 * trap format vs. data frame count mismatch
 *
 * Revision 1.14  2006/09/25  03:22:40  avi
 * Solaris 10
 *
 * Revision 1.13  2006/08/07  03:30:35  avi
 * safe 1.e-37 test
 *
 * Revision 1.12  2006/05/04  02:17:19  avi
 * -W does not imply -D
 *
 * Revision 1.11  2006/05/04  01:46:49  avi
 * option -W (interpret mask a weights)
 *
 * Revision 1.10  2005/09/19  23:22:25  avi
 * options -D -d -f
 *
 * Revision 1.9  2005/09/16  03:43:29  avi
 * conc capability
 *
 * Revision 1.8  2005/09/16  01:37:33  avi
 * separate image and mask operations in preparation for conc capability
 *
 * Revision 1.7  2005/08/01  19:46:07  jon
 * Change total to type double to avoid precision problems.
 *
 * Revision 1.6  2004/09/22  05:51:09  avi
 * use get_4dfp_dimo_quiet ()
 *
 * Revision 1.5  2004/09/21  21:24:24  rsachs
 * Installed 'errm','errr','errw','setprog'. Replaced 'Get4dfpDimN' with 'get_4dfp_dimo'.
 *
 * Revision 1.4  2004/05/26  23:56:01  avi
 * -v option
 *
 * Revision 1.3  2000/04/22  03:18:08  avi
 * time series mode
 *
 * Revision 1.2  2000/04/07  15:18:40  jon
 * Updated to correct file reading.
 *
 * Revision 1.1  2000/03/30  17:34:01  jon
 * Initial revision
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef _FEATURES_H		/* defined by math.h in Linux */
	#include <sunmath.h>	/* isnormal() */
	#include <ieeefp.h>
#endif
#include <endianio.h>
#include <conc.h>

/*************/
/* externals */
/*************/
extern int expandf (char *format, int bufsiz);		/* expandf.c */

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void getrange (char *string, float *minval, float *maxval) {
        char	*str;

	str = strstr (string, "to");
	if (str) {
		*str = '\0';
		*minval = atof (string);
		*maxval = atof (str + 2);
	} else {
		*minval = 0.0;
		*maxval = atof (string);
	}
}

/********************/
/* global variables */
/********************/
char	program[MAXL];
static char rcsid[]= "$Id: qnt_4dfp.c,v 1.19 2010/05/07 03:12:23 avi Exp $";

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
}

static void usage (char *program) {
	printf ("%s\n", rcsid);
	printf ("Usage:\t%s <(4dfp)|(conc) image> <(4dfp) mask>\n", program);
	printf (" e.g.:\t%s -t23.2 va1234_mpr mask\n", program);
	printf ("\toption\n");
	printf ("\t-s\ttime series mode\n");
	printf ("\t-d\tinclude backwards differences (differentiated signal) in output\n");
	printf ("\t-D\tcount only defined <image> voxels (see note below)\n");
	printf ("\t-A\tapply threshold test to absolute value of <mask>\n");
	printf ("\t-W\tinterpret <mask> as spatial weights (negative values allowed)\n");
	printf ("\t-v<flt>[to<flt>] count only <image> voxels within specified range\n");
	printf ("\t-f<str>\tspecify frames to count format, e.g., \"4x120+4x76+\"\n");
	printf ("\t-p<flt>\tspecify mask threshold as percent of <mask> max\n");
	printf ("\t-t<flt>\tspecify absolute <mask> threshold (default = 0.0)\n");
	printf ("\t-c<flt>\tscale output mean values by specified constant (default = 1.0)\n");
	printf ("N.B.:\tonly the first frame of <mask> is used\n");
	printf ("N.B.:\t<image> and <mask> may be the same\n");
	printf ("N.B.:\tconc files must have extension \"conc\"\n");
	printf ("N.B.:\tdefined means not 0.0 and not NaN and not 1.e-37 and finite\n");
	printf ("N.B.:\t-d requires -f and implies -s\n");
	printf ("N.B.:\t-W disables mask threshold testing\n");
	exit (1);
}

int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	CONC_BLOCK	conc_block;			/* conc i/o control block */
	FILE            *imgfp, *mskfp;
	char            imgroot[MAXL], imgfile[MAXL], mskroot[MAXL], mskfile[MAXL];

/**************/
/* processing */
/**************/
	char		*format;			/* pattern of frames to count */
	int             imgdim[4], mskdim[4], dimension, nvox, imgori, mskori, isbig, isbigm;
	int		nf_total = 0, nfile, ifile, iframe;
	float           imgvox[3], mskvox[3];
        float           *imgt, *imgm;
	float		threshold = 0.0, pct, imgmax, vneg, vpos;
	double		total, scale = 1.0;

/***********/
/* utility */
/***********/
	char            *ptr, command[1024];
	int             c, i, j, k;
	double		q, p, w, diff;

/*********/
/* flags */
/*********/
	int		wei_flag = 0;
	int		conc_flag = 0;
	int		deriv_flag = 0;
	int             pctflag = 0;
	int             absflag = 0;
	int		series = 0;
	int             vrange = 0;
	int             status = 0;
	int		defined;	/* test on imgt voxel value */
	int		D_flag = 0;	/* gates defined test effect on ROI value */
	int             pdefined;	/* previous frame is defined */

	setprog (program, argv);
/*	q = 1./0.; printf ("isnormal(%f)=%d\n", q, isnormal(q));	*/

	if (!(format = (char *) calloc (1024, sizeof (char)))) errm (program);
/************************/
/* process command line */
/************************/
	for (j = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'd': deriv_flag++;
				case 's': series++;					break;
				case 'D': D_flag++;					break;
				case 'W': wei_flag++;					break;
				case 'A': absflag++;					break;
				case 'f': strncpy (format, ptr, 1023);			*ptr = '\0'; break;
				case 'c': scale = atof (ptr);				*ptr = '\0'; break;
				case 'p': pctflag++; pct = atof (ptr);			*ptr = '\0'; break;
				case 't': threshold = atof (ptr);			*ptr = '\0'; break;
				case 'v': getrange (ptr, &vneg, &vpos);	vrange++;	*ptr = '\0'; break;
				default:						break;
			}
		} else switch (j) {
			case 0:	getroot (argv[i], imgroot);
				conc_flag = (strstr (argv[i], ".conc") == argv[i] + strlen (imgroot));
                						j++; break;
                	case 1: getroot (argv[i], mskroot);	j++; break;
			default:				break;
		}
	}
	if (j < 2) usage (program);
	write_command_line (stdout, argc, argv);

/********************************************/
/* get input 4dfp dimensions and open files */
/********************************************/
	if (conc_flag) {
		conc_init_quiet (&conc_block, program);
		conc_open_quiet (&conc_block, imgroot);
		strcpy (imgfile, conc_block.lstfile);
		for (k = 0; k < 4; k++) imgdim[k] = conc_block.imgdim[k];
		for (k = 0; k < 3; k++) imgvox[k] = conc_block.voxdim[k];
		imgori = conc_block.orient;
		nfile  = conc_block.rnfile;
		isbig = conc_block.isbig;
	} else {
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		if (get_4dfp_dimoe_quiet (imgfile, imgdim, imgvox, &imgori, &isbig) < 0) errr (program, imgfile);
		nfile = 1;
		if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
	}

/***********************************/
/* set up frames to count (format) */
/***********************************/
	if (!(format = (char *) realloc (format, (imgdim[3] + 1)*sizeof (char)))) errm (program);
	if (strlen (format)) {
		if (k = expandf (format, imgdim[3] + 1)) exit (k);
		if ((nf_total = strlen (format)) != imgdim[3]) {
			fprintf (stderr, "format codes for %s frames (%d) than data (%d)\n",
				(nf_total < imgdim[3]) ? "fewer" : "more", nf_total, imgdim[3]);
			exit (-1);
		}
	} else {
		for (k = 0; k < imgdim[3]; k++) format[k] = '+';
		format[imgdim[3]] = '\0';
	}
	if (deriv_flag && !nf_total) usage (program);

/*****************/
/* alloc buffers */
/*****************/
	dimension  = imgdim[0]*imgdim[1]*imgdim[2];
	if (!(imgt = (float *) malloc (dimension * sizeof (float)))) errm (program);
	if (!(imgm = (float *) malloc (dimension * sizeof (float)))) errm (program);

/****************/
/* prepare mask */
/****************/
	sprintf (mskfile, "%s.4dfp.img", mskroot);
	if (get_4dfp_dimoe_quiet (mskfile, mskdim, mskvox, &mskori, &isbigm) < 0) errr (program, mskfile);
	status = mskori - imgori;
	for (k = 0; k < 3; k++) status |= (imgdim[k] != mskdim[k]);
	for (k = 0; k < 3; k++) status |= (fabs (imgvox[k] - mskvox[k]) > 0.0001);
	if (status) {
		fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot, mskroot);
		exit (-1);
	}
	if (!(mskfp = fopen (mskfile, "rb"))) errr (program, mskfile);
	fprintf (stdout, "#Reading: %s\n", mskfile);
	if (eread (imgm, dimension, isbigm, mskfp)) errr (program, mskfile);
	if (absflag) for (i = 0; i < dimension; i++) if (imgm[i] < 0.) imgm[i] = -imgm[i];
	if (pctflag) {
		imgmax = -FLT_MAX;
		for (i = 0; i < dimension; i++) if (imgm[i] > imgmax) imgmax = imgm[i];
		printf ("#Maximum= %f\n#Percent= %f\n", imgmax, pct);
		threshold = pct * imgmax / 100.;
		rewind (mskfp);
	}
	printf ("#Threshold= %-.4f\n", threshold);

/***********/
/* process */
/***********/
	printf ("#Reading: %s\n", imgfile);
	if (series) {
		if (deriv_flag) {
			printf ("#Frame      Mean     Deriv      Nvox\n");
		} else {
			printf ("#Frame      Mean      Nvox\n");
		}
	}
	for (iframe = pdefined = ifile = 0; ifile < nfile; ifile++) {
		for (k = 0; k < ((conc_flag) ? conc_block.nvol[ifile] : imgdim[3]); k++, iframe++) {
			if (conc_flag) {
				conc_read_vol (&conc_block, imgt);
			} else {
				if (eread (imgt, dimension, isbig, imgfp)) errr (program, imgfile);
			}
			total = nvox = 0;
			for (i = 0; i < dimension; i++) {
				j = (vrange) ? (imgt[i] > vneg && imgt[i] < vpos) : 1;
				defined = isnormal (imgt[i]) && imgt[i] != (float) 1.e-37;
				if (D_flag) j &= defined;
				if (j && (wei_flag || (imgm[i] > threshold))) {
					w = (wei_flag) ? imgm[i] : 1.0;
					total += imgt[i]*w;
					nvox++;
				}
			}
			p = q;
			q = scale*total/nvox;
/*****************************************************/
/* print frame, mean, and voxel count for each frame */
/*****************************************************/
			if (series) {
				pdefined = k && format[iframe - 1] != 'x';
				if (deriv_flag) {
					diff = (pdefined) ? q - p : 0.0;
					printf ("%6d%10.4f%10.4f%10d\n", iframe + 1, q, diff, nvox);
				} else {
					printf ("%6d%10.4f%10d\n", iframe + 1, q, nvox);
				}
			} else {
				printf ("Total= %f\nMean= %f\nVoxels= %d\n", total, q, nvox);
			}
		}
	}

/*********************/
/* clean up and exit */
/*********************/
	if (conc_flag) {
		conc_free (&conc_block);
	} else {
		if (fclose (imgfp)) errr (program, imgfile);
	}
	fclose (mskfp);
	free (imgt); free (imgm);
	free (format);
	exit (status);
}

