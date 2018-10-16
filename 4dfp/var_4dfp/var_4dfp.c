/*$Header: /data/petsun4/data1/src_solaris/var_4dfp/RCS/var_4dfp.c,v 1.20 2010/12/23 03:36:23 avi Exp $*/
/*$Log: var_4dfp.c,v $
 * Revision 1.20  2010/12/23  03:36:23  avi
 * long -> int in rnan()
 *
 * Revision 1.19  2010/10/01  01:38:49  avi
 * correct serious error in coding of format pointer
 * option -G
 *
 * Revision 1.18  2010/07/08  19:34:39  avi
 * MAXF inreased to 4096
 *
 * Revision 1.17  2008/02/06  04:05:36  avi
 * correct crash on attempted read of long format string from command line
 *
 * Revision 1.16  2007/08/14  03:24:36  avi
 * format length now input determined
 *
 * Revision 1.15  2007/08/14  01:36:30  avi
 * many more controls on format in relation to total data length
 *
 * Revision 1.14  2007/08/05  01:22:29  avi
 * linux compliant (code rnan())
 *
 * Revision 1.13  2006/09/24  02:20:47  avi
 * Solaris 10
 *
 * Revision 1.12  2006/08/07  03:35:20  avi
 * safe 1.e-37 test
 *
 * Revision 1.11  2005/08/21  02:10:47  avi
 * output unbiased (1/(n-1)) variance estimate
 * mark_undefined and options -E -N -Z
 *
 * Revision 1.10  2005/08/20  03:04:33  avi
 * read/write conc
 *
 * Revision 1.9  2004/11/03  05:33:38  avi
 * eliminate dependence of lmri
 * frames to count format
 * eliminate support for pre-ifh 4dfp
 *
 * Revision 1.8  2003/06/28  02:53:40  avi
 * -z option
 *
 * Revision 1.7  2000/04/25  06:17:43  avi
 * -m (remove mean volume) mode
 *
 * Revision 1.6  2000/04/25  05:09:31  avi
 * eliminate FORTRAN
 * allocate one frame instead of whole stack
 *
 * Revision 1.5  1997/05/23  03:56:17  yang
 * new rec macros
 * Revision 1.4  1997/04/28  21:10:15  yang
 * Working Solaris version.
 * Revision 1.3  1996/10/25  02:41:40  avi
 * add -c switch
 * Revision 1.2  1996/08/05  22:29:50  avi
 * correct output in old (pre ifh) mode
 * Revision 1.1  1996/08/05  02:13:36  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>		/* defines isnan() in Linux */
#include <float.h>
#ifndef _FEATURES_H		/* defined by math.h in Linux */
	#include <ieeefp.h>	/* isnan() */
#endif
#include <assert.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>
#include <conc.h>		/* /data/petsun4/data1/src_solaris/actmapf_4dfp */

#define MAXL	256
#define MAXF	4096		/* format characters as input from command line */

/*************/
/* externals */
/*************/
int	expandf (char *string, int len);		/* expandf.c */

/********************/
/* global variables */
/********************/
static char rcsid[] = "$Id: var_4dfp.c,v 1.20 2010/12/23 03:36:23 avi Exp $";
static char	program[MAXL];
static int	debug = 0;

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

float rnan (void) {
	union {
		float		r;
		unsigned int	j;
        } word;
	word.j = 0x7fffffff;
	return word.r;
}

void mark_undefined (float *imgt, int dimension, float *imgs, int NaN_flag) {
	int		i;

	for (i = 0; i < dimension; i++) if (!imgs[i]) switch (NaN_flag) {
		case 'Z': imgt[i] = 0.0;	break;
		case 'E': imgt[i] = 1.e-37;	break;
		case 'N': imgt[i] = rnan ();	break;
		default:			break;
	}
}

int main (int argc, char *argv[]) {
	FILE		*imgfp, *outfp;
	char		imgroot[MAXL], imgfile[MAXL];
	char		outroot[MAXL], outfile[MAXL], ifhfile[MAXL];
	CONC_BLOCK	conc_block;		/* conc i/o control block */
	IFH		ifh;

/*********/
/* flags */
/*********/
	int		opr = 'v';		/* 'v' : variance; 's' : sd1; 'm' : remove mean */
	int		conc_flag = 0;
	int		global = 0;		/* when set compute voxelwise mean ignoring run boundaries */
	int		NaN_flag = 'E';		/* 'E' 1.e-37; 'Z' 0.0; 'N' NaN; */
	int		status = 0;

/***********/
/* utility */
/***********/
	double		q;
	int		c, i, j, k, kk;
	char		*ptr, command[1024];		/* capacity must accommodate long format string */

/*******************/
/* imge processing */
/*******************/
	char		*format;			/* pattern of frames to count */
	float		*fptr;				/* general float pointer */
	float		*imgt;				/* input 4dfp volume */
	float		*imgs;				/* mask of nonzero (sampled) voxels */
	float		*imgv;				/* variance over all frames */
	float		*imga;				/* average over all frames */
	float		**imgu;				/* multivolume average over all frames */
	float		factor = 1.0;			/* multiplier for output image values */
	int		imgdim[4], dimension, nf_total, nf_func, nf_anat = 0;
	int		nfile, ifile, iframe;
	int		isbig;
	char		control = '\0';

	fprintf (stdout, "%s\n", rcsid);
	if (ptr = strrchr (argv[0], '/')) ptr++; else ptr = argv[0];
	strcpy (program, ptr);

	if (!(format = (char *) calloc (MAXF + 1, sizeof (char)))) errm (program);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'd': debug++;				break;
				case 'G': global++;				break;
				case 's': case 'm': case 'z': opr = c;		break;
				case 'N': case 'Z': case 'E': NaN_flag = c;	break;
				case '@': control = *ptr++;			*ptr = '\0'; break;
				case 'f': strncpy (format, ptr, MAXF);		*ptr = '\0'; break;
				case 'n': nf_anat = atoi (ptr);			*ptr = '\0'; break;
				case 'c': factor = atof (ptr);			*ptr = '\0'; break;
				default:					break;
			}
		} else switch (k) {
		 	case 0: getroot (argv[i], imgroot);
				conc_flag = (strstr (argv[i], ".conc") == argv[i] + strlen (imgroot));
										k++; break;
			default:						break;
		}
	}
	if (k < 1) {
		printf ("Usage:\t%s <(4dfp|conc) input>\n", program);
		printf (" e.g.,\t%s -sn3 -c10 test_b1_rmsp_dbnd\n", program);
		printf ("\toption\n");
		printf ("\t-d\tdebug mode\n");
		printf ("\t-m\tremove mean volume from stack\n");
		printf ("\t-s\tcompute s.d. about mean\n");
		printf ("\t-G\tcompute mean ignoring run boundaries (default within runs)\n");
		printf ("\t-v\tcompute variance about mean (default operation)\n");
		printf ("\t-z\toutput logical and of all nonzero defined voxels\n");
		printf ("\t-n<int>\tspecify number of pre-functional frames per run (default = 0)\n");
		printf ("\t-f<str>\tspecify frames to count format, e.g., \"4x120+4x76+\"\n");
		printf ("\t-c<flt>\tscale output image values by specified factor\n");
		printf ("\t-N\toutput undefined voxels as NaN\n");
		printf ("\t-Z\toutput undefined voxels as 0\n");
		printf ("\t-E\toutput undefined voxels as 1.e-37 (default)\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:\tinput conc files must have extension \"conc\"\n");
		printf ("N.B.:\tvoxelwise mean is individually computed over each run in conc\n");
		printf ("N.B.:\t-f option overrides -n\n");
		exit (1);
	}

/***********************************/
/* get input 4dfp image dimensions */
/***********************************/
	if (conc_flag) {
		conc_init (&conc_block, program);
		conc_open (&conc_block, imgroot);
		strcpy (imgfile, conc_block.lstfile);
		for (k = 0; k < 4; k++) imgdim[k] = conc_block.imgdim[k];
		dimension = conc_block.vdim;
		nfile = conc_block.rnfile;
		isbig = conc_block.isbig;
	} else {
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		if (Getifh (imgfile, &ifh)) errr (program, imgfile);
		for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];
		dimension = imgdim[0] * imgdim[1] * imgdim[2];
		nfile = 1;
		if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
		isbig = strcmp (ifh.imagedata_byte_order, "littleendian");
	}
	if (!control) control = (isbig) ? 'b' : 'l';

/***********************************/
/* set up frames to count (format) */
/***********************************/
	if (!(format = (char *) realloc (format, (imgdim[3] + 1)*sizeof (char)))) errm (program);
	if (strlen (format)) {
		if (k = expandf (format, imgdim[3] + 1)) exit (k);
		if ((nf_total = strlen (format)) != imgdim[3]) {
			fprintf (stderr, "format codes frame count (%d) not equal to data length (%d)\n", nf_total, imgdim[3]);
			exit (-1);
		}
	} else {
		nf_total = imgdim[3];
		printf ("var_4dfp: nf_anat=%d\n", nf_anat);
		for (k = ifile = 0; ifile < nfile; ifile++) {
			nf_func = ((conc_flag) ? conc_block.nvol[ifile] : imgdim[3]) - nf_anat;
			if (nf_func < 0) {
				fprintf (stderr, "%s: %s has more skip than total frames\n", program,
				(conc_flag) ? conc_block.imgfile0[ifile] : imgfile);
			}
			for (j = 0; j < nf_anat; j++) format[k++] = 'x';
			for (j = 0; j < nf_func; j++) format[k++] = '+';
		}
		format[k] = '\0';
		if (k != imgdim[3]) {
			fprintf (stderr, "k != imgdim[3] (should never happen)\n");
			exit (-1);
		}
	}
	printf ("%s\n", format);
	for (nf_func = j = 0; j < nf_total; j++) if (format[j] == '+') nf_func++;
	printf ("frames total=%d counted=%d skipped=%d\n", nf_total, nf_func, nf_total - nf_func);

	if (!(imgt = (float *) malloc (dimension * sizeof (float)))
	||  !(imgs = (float *) malloc (dimension * sizeof (float)))
	||  !(imga = (float *) calloc (dimension,  sizeof (float)))
	||  !(imgv = (float *) calloc (dimension,  sizeof (float)))) errm (program);
	imgu = calloc_float2 (nfile, dimension);

/**************/
/* initialize */
/**************/
	for (i = 0; i < dimension; i++) imgs[i] = 1.0;

/********************************/
/* compute mean and logical and */
/********************************/
	for (iframe = kk = ifile = 0; ifile < nfile; ifile++) {
		printf ("Reading: %s", (conc_flag) ? conc_block.imgfile0[ifile] : imgfile);
		printf ("\t%d frames\n", (conc_flag) ? conc_block.nvol[ifile] : imgdim[3]);
		for (k = j = 0; j < ((conc_flag) ? conc_block.nvol[ifile] : imgdim[3]); j++) {
			if (conc_flag) {
				conc_read_vol (&conc_block, imgt);
			} else {
				if (eread (imgt, dimension, isbig, imgfp)) errr (program, imgfile);
			}
			if (format[iframe++] != '+') continue;
			k++;		/* single run functional frame count */
			kk++;		/* global     functional frame count */
			for (i = 0; i < dimension; i++) {
				if (imgt[i] == (float) 1.e-37 || isnan (imgt[i]) || imgt[i] == 0.0) imgs[i] = 0.;
				imgu[ifile][i]	+= imgt[i];
				imga[i]		+= imgt[i];
			}
		}
		for (i = 0; i < dimension; i++) imgu[ifile][i] /= k;
	}
	for (i = 0; i < dimension; i++) imga[i] /= kk;
	assert (iframe == imgdim[3]);

/********************/
/* compute variance */
/********************/
	if (!conc_flag) rewind (imgfp);
	printf ("Reading files ");
	for (iframe = ifile = 0; ifile < nfile; ifile++) {printf (" %d", ifile + 1); fflush (stdout);
		for (j = 0; j < ((conc_flag) ? conc_block.nvol[ifile] : imgdim[3]); j++) {
			if (conc_flag) {
				conc_read_vol (&conc_block, imgt);
			} else {
				if (eread (imgt, dimension, isbig, imgfp)) errr (program, imgfile);
			}
			if (format[iframe++] != '+') continue;
			for (i = 0; i < dimension; i++) {
				if (global) {
					q = imgt[i] - imga[i];
				} else {
					q = imgt[i] - imgu[ifile][i];
				}
				imgv[i] += q*q;
			}
		}
	}
	printf ("\n");
	for (i = 0; i < dimension; i++) imgv[i] /= (global) ? kk : (nf_func - nfile);

	switch (opr) {
	case 'z':
		sprintf (outroot, "%s_sam", imgroot);
		fptr = imgs;
		break;
	case 'm':
		sprintf (outroot, "%s_uout", imgroot);
		break;
	case 's':
		sprintf (outroot, "%s_sd1", imgroot);
		for (k = 0; k < dimension; k++) imgv[k] = factor * sqrt (imgv[k]);
		mark_undefined (imgv, dimension, imgs, NaN_flag);
		fptr = imgv;
		break;
	case 'v': default:
		sprintf (outroot, "%s_var", imgroot);
		for (k = 0; k < dimension; k++) imgv[k] *= factor;
		mark_undefined (imgv, dimension, imgs, NaN_flag);
		fptr = imgv;
		break;
	}
				
/****************/
/* write output */
/****************/
	switch (opr) {
	case 'm':
		if (conc_flag) {
			conc_newe (&conc_block, "uout", control);
			sprintf (outfile, "%s.conc", outroot);
			printf ("Writing: %s\n", outfile);
		} else {
			sprintf (outfile, "%s.4dfp.img", outroot);
			if (!(outfp = fopen (outfile, "wb"))) errw (program, outfile);
			rewind (imgfp);
		}
		for (ifile = 0; ifile < nfile; ifile++) {
			printf ("Writing: %s", (conc_flag) ? conc_block.imgfile1[ifile] : outfile);
			printf ("\t%d frames\n", (conc_flag) ? conc_block.nvol[ifile] : imgdim[3]);
			for (j = 0; j < ((conc_flag) ? conc_block.nvol[ifile] : imgdim[3]); j++) {
				if (conc_flag) {
					conc_read_vol (&conc_block, imgt);
				} else {
					if (eread (imgt, dimension, isbig, imgfp)) errr (program, imgfile);
				}
				for (i = 0; i < dimension; i++) {imgt[i] -= imgu[ifile][i]; imgt[i] *= factor;}
				mark_undefined (imgt, dimension, imgs, NaN_flag);
				if (conc_flag) {
					conc_write_vol (&conc_block, imgt);
				} else {
					if (ewrite (imgt, dimension, control, outfp)) errw (program, outfile);
				}
			}
		}
		if (!conc_flag) if (fclose (outfp)) errw (program, outfile);
		break;
	case 'z': case 'v': case 's':
		imgdim[3] = 1;
		sprintf (outfile, "%s.4dfp.img", outroot);
		printf ("Writing: %s\n", outfile);
		if (!(outfp = fopen (outfile, "wb"))) errw (program, outfile);
		if (ewrite (fptr, dimension, control, outfp)) errw (program, outfile);
		if (fclose (outfp)) errw (program, outfile);
		break;
	}

/***************/
/* ifh and hdr */
/***************/
	if (opr == 'm' && conc_flag) {
		conc_ifh_hdr_rec (&conc_block, argc, argv, rcsid);
	} else {
		if (conc_flag) {
			writeifhmce (program, outfile, imgdim,
				conc_block.voxdim, conc_block.orient, conc_block.mmppix, conc_block.center, control);
		} else {
			writeifhmce (program, outfile, imgdim,
				ifh.scaling_factor, ifh.orientation, ifh.mmppix, ifh.center, control);
		}
		sprintf (command, "ifh2hdr %s", outroot); status |= system (command);
	}

/*******/
/* rec */
/*******/
	startrece (outfile, argc, argv, rcsid, control);
	sprintf (command, "Output scaled by %.4f\n", factor); printrec (command);
	if (conc_flag) {
		for (ifile = 0; ifile < nfile; ifile++) catrec (conc_block.imgfile0[ifile]);
	} else {
		catrec (imgfile);
	}
	endrec ();

/*************/
/* close i/o */
/*************/
	if (conc_flag) {
		conc_free (&conc_block);
	} else {
		if (fclose (imgfp)) errr (program, imgfile);
	}

	free (imgt); free (imgv); free (imgs); free (imga);
	free_float2 (imgu);
	free (format);
	exit (status);
}
