/*$Header: /data/petsun4/data1/src_solaris/qnt_4dfp/RCS/qntm_4dfp.c,v 1.8 2012/06/22 03:59:06 avi Exp $*/
/*$Log: qntm_4dfp.c,v $
 * Revision 1.8  2012/06/22  03:59:06  avi
 * option -h
 *
 * Revision 1.7  2010/11/30  04:46:00  avi
 * trap empty ROI files
 *
 * Revision 1.6  2010/10/19  02:48:44  avi
 * integerize code-by-value logic
 *
 * Revision 1.5  2010/08/28  04:38:30  avi
 * option -Z (count zero timeseries voxels)
 *
 * Revision 1.4  2010/05/20  01:22:37  avi
 * option -V
 *
 * Revision 1.3  2010/05/19  22:33:30  avi
 * minor safety code
 *
 * Revision 1.2  2010/05/19  02:36:10  avi
 * change ROI index memory strategy from imgn based to voxdat based
 *
 * Revision 1.1  2010/05/19  00:49:21  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <float.h>
#include <assert.h>
#ifndef _FEATURES_H		/* defined by math.h in Linux */
	#include <sunmath.h>	/* isnormal() */
	#include <ieeefp.h>
#endif
#include <endianio.h>
#include <conc.h>		/* includes definition of MAXL */

#define	UCHAR	unsigned char
#define	USHORT	unsigned short
#define INCR	4		/* ROI index memory increment size in items */
				/* must be a divisor of sizeof(UCHAR) */

typedef struct {
	USHORT		*p;	/* pointer to ROI index list */
	UCHAR		m;	/* memory allocated for this many items */
	UCHAR		n;	/* items entered into list */
} VOXDAT;

/********************/
/* global variables */
/********************/
char	program[MAXL];
static char rcsid[]= "$Id: qntm_4dfp.c,v 1.8 2012/06/22 03:59:06 avi Exp $";

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
	printf ("\t-Z\tcount zero voxels in <image> as defined\n");
	printf ("\t-V\tforce code_by_volume even if the number of volumes is 1\n");
	printf ("\t-h\tsuppress printing output header\n");
	printf ("\t-o<str>\twrite output to specified text file (default stdout)\n");
	printf ("N.B.:\tconc files must have extension \"conc\"\n");
	printf ("N.B.:\tonly defined voxels (not 0.0 and not NaN and not 1.e-37 and finite) are counted\n");
	printf ("N.B.:\t<(4dfp) ROI> may either a value-coded single volume ROI image or a multi-volume mask\n");
	printf ("N.B.:\t<(4dfp) ROI> coded values are integerized\n");
	printf ("N.B.:\t%s ignores <(4dfp) ROI> ifh center and mmppix fields\n", program);
	exit (1);
}

int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	CONC_BLOCK	conc_block;			/* conc i/o control block */
	FILE            *imgfp, *ROIfp, *outfp;
	char            imgroot[MAXL], imgfile[MAXL], ROIroot[MAXL], ROIfile[MAXL], outfile[MAXL] = "";

/**************/
/* processing */
/**************/
	int             imgdim[4], ROIdim[4], dimension, imgori, ROIori, isbig, isbigm;
	int		nROI, *nvox;
	float           imgvox[3], ROIvox[3];
        float           *imgt, *imgr;
	int		jmin, jmax;
	double		*total;
	VOXDAT		*voxdat;

/***********/
/* utility */
/***********/
	char            *ptr, command[2*MAXL];
	int             c, i, j, k, l;

/*********/
/* flags */
/*********/
	int		print_header = 1;
	int		conc_flag = 0;
	int             status = 0;
	int		count_zero = 0;
	int		code_by_volume = 0;
	int		defined, definedm;	/* test on imgt and imgr voxel values */

	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (j = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'Z': count_zero++;			break; 
				case 'V': code_by_volume++;		break; 
				case 'h': print_header = 0;		break; 
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

	if (strlen (outfile)) {
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	} else {
		outfp = stdout;
	}
	if (print_header) write_command_line (outfp, argc, argv);

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
		isbig  = conc_block.isbig;
	} else {
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		if (get_4dfp_dimoe_quiet (imgfile, imgdim, imgvox, &imgori, &isbig) < 0) errr (program, imgfile);
		if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
	}

/*****************/
/* alloc buffers */
/*****************/
	dimension  = imgdim[0]*imgdim[1]*imgdim[2];
	if (!(imgt = (float *) malloc (dimension*sizeof (float)))) errm (program);
	if (!(imgr = (float *) malloc (dimension*sizeof (float)))) errm (program);
	if (!(voxdat = (VOXDAT *) calloc (dimension, sizeof (VOXDAT)))) errm (program);
	for (i = 0; i < dimension; i++) {
		if (!(voxdat[i].p = (USHORT *) calloc (INCR, sizeof (USHORT)))) errm (program);
		voxdat[i].m = INCR;
	}

/****************/
/* prepare ROIs */
/****************/
	sprintf (ROIfile, "%s.4dfp.img", ROIroot);
	if (get_4dfp_dimoe_quiet (ROIfile, ROIdim, ROIvox, &ROIori, &isbigm) < 0) errr (program, ROIfile);
	status = ROIori - imgori;
	for (k = 0; k < 3; k++) status |= (imgdim[k] != ROIdim[k]);
	for (k = 0; k < 3; k++) status |= (fabs (imgvox[k] - ROIvox[k]) > 0.0001);
	if (status) {
		fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot, ROIroot);
		exit (-1);
	}

	if (ROIdim[3] > 1) code_by_volume++;
	fprintf (stdout, "Reading: %s\n", ROIfile);
	if (!(ROIfp = fopen (ROIfile, "rb"))) errr (program, ROIfile);
	if (code_by_volume) {
		if (print_header) fprintf (outfp, "#code ROI by volume\n");
		nROI = ROIdim[3];
		for (j = 0; j < nROI; j++) {
			if (eread (imgr, dimension, isbigm, ROIfp)) errr (program, ROIfile);
			for (i = 0; i < dimension; i++) {
				definedm = isnormal (imgr[i]) && imgr[i] != (float) 1.e-37 && imgr[i] != 0.;
				if (definedm) {
					if (voxdat[i].n == voxdat[i].m) {
						if (!(voxdat[i].m += INCR)) {
							fprintf (stderr, "%s: ROI/voxel limit (256) exceeded\n", program);
							exit (-1);
						}
						voxdat[i].p = realloc (voxdat[i].p, voxdat[i].m*sizeof (USHORT));
						if (!voxdat[i].p) errm (program);
					}
					voxdat[i].p[voxdat[i].n++] = j;
				}
			}
		}
	} else {	/* code by value */
		nROI = 0;
		if (print_header) fprintf (outfp, "#code ROI by value\n");
		if (eread (imgr, dimension, isbigm, ROIfp)) errr (program, ROIfile);
		jmin = 1000000; jmax = -jmin;
		for (i = 0; i < dimension; i++) {
			definedm = isnormal (imgr[i]) && imgr[i] != (float) 1.e-37;
			if (!definedm) continue;
			j = (int) imgr[i] + FLT_EPSILON; if (!j) continue;
			if (j < jmin) jmin = j;
			if (j > jmax) jmax = j;
			nROI = jmax - jmin + 1;
		}
		if (!nROI) {
			fprintf (stderr, "%s: no ROIs in %s\n", program, ROIfile);
			exit (-1);
		}
		if (print_header) fprintf (outfp, "#%s range %d to %d; nROI = %d\n", ROIfile, jmin, jmax, nROI);
		
		for (i = 0; i < dimension; i++) {
			definedm = isnormal (imgr[i]) && imgr[i] != (float) 1.e-37;
			if (!definedm) continue;
			j =  (int) (imgr[i] + FLT_EPSILON); if (!j) continue;
			j -= jmin;
			assert (j >= 0);
			voxdat[i].p[voxdat[i].n++] = j;
		}
	}
	if (fclose (ROIfp)) errr (program, ROIfile);

	if (!(total = (double *) malloc (nROI*sizeof (double)))) errm (program);
	if (!(nvox  =    (int *) malloc (nROI*sizeof (int))))    errm (program);
/***********/
/* process */
/***********/
	printf ("Reading: %s\n", imgfile);
	for (k = 0; k < imgdim[3]; k++) {
		if (conc_flag) {
			conc_read_vol (&conc_block, imgt);
		} else {
			if (eread (imgt, dimension, isbig, imgfp)) errr (program, imgfile);
		}
		for (j = 0; j < nROI; j++) total[j] = nvox[j] = 0;
		for (i = 0; i < dimension; i++) {
			defined =  isnormal (imgt[i]) && imgt[i] != (float) 1.e-37;
			if (count_zero && imgt[i] == 0.) defined = 1;
			if (!defined) continue;
			for (l = 0; l < voxdat[i].n; l++) {
				j = voxdat[i].p[l];
				total[j] += imgt[i];
				nvox[j]++;
			}
		}
		for (j = 0; j < nROI; j++) {
			if (!nvox[j]) {
				fprintf (stderr, "%s: zero size ROI detected\n", program);
				status = 1;
			} else {
				fprintf (outfp, " %9.4f", total[j]/nvox[j]);
			}
		}
		fprintf (outfp, "\n");
	}
	if (strlen (outfile)) if (fclose (outfp)) errw (program, outfile);

	if (code_by_volume) {
/***************************************/
/* report max number of ROIs per voxel */
/***************************************/
		for (k = i = 0; i < dimension; i++) if (voxdat[i].n > k) k = voxdat[i].n;
		printf ("%s: maximum number of ROIs per voxel = %d\n", program, k);
	}

/*********************/
/* clean up and exit */
/*********************/
	if (conc_flag) {
		conc_free (&conc_block);
	} else {
		if (fclose (imgfp)) errr (program, imgfile);
	}
	for (i = 0; i < dimension; i++) free (voxdat[i].p);
	free (voxdat);
	free (total); free (nvox);
	free (imgt); free (imgr);
	exit (status);
}
