/*$Header: /data/petsun4/data1/src_solaris/qnt_4dfp/RCS/qntw_4dfp.c,v 1.1 2011/04/18 07:25:23 avi Exp $*/
/*$Log: qntw_4dfp.c,v $
 * Revision 1.1  2011/04/18  07:25:23  avi
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
#include <rec.h>
#include <endianio.h>
#include <conc.h>		/* includes definition of MAXL */

/********************/
/* global variables */
/********************/
char	program[MAXL];
static char rcsid[]= "$Id: qntw_4dfp.c,v 1.1 2011/04/18 07:25:23 avi Exp $";

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

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

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;
	char            string[MAXL];

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
        get_time_usr (string);
        fprintf (outfp, "#%s", string);
        if (get_machine_info (string)) {
                fprintf (outfp, "@%s\n", string);
        } else {
                fprintf (outfp, "\n");
        }
}

static void usage (char *program) {
	printf ("%s\n", rcsid);
	printf ("Usage:\t%s <(4dfp)|(conc) image> <(4dfp) ROI>\n", program);
	printf (" e.g.:\t%s TC30274_rmsp_faln_dbnd_xr3d_atl.conc iter10_roi_-02_-37_+27m_ROI\n", program);
	printf ("\toption\n");
	printf ("\t-L<int>\tspecify ROI weight L-norm (default = 0)\n");
	printf ("\t-o<str>\twrite output to specified text file (default stdout)\n");
	printf ("\t-Z\tcount zero voxels in <image> as defined\n");
	printf ("\t-H\tinclude heaer info in output\n");
	printf ("N.B.:\tconc files must have extension \"conc\"\n");
	printf ("N.B.:\t<(4dfp) ROI> is interpreted as a multi-volume voxel-wise set of weights\n");
	printf ("N.B.:\tonly defined voxels (not 0.0 and not NaN and not 1.e-37 and finite) are counted\n");
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
	int		nROI, norm = 0;
	float           imgvox[3], ROIvox[3];
        float           *imgt;
	double		*total, **imgr;

/***********/
/* utility */
/***********/
	char            *ptr, command[2*MAXL];
	int             c, i, j, k, n;
	double		q, t;

/*********/
/* flags */
/*********/
	int		H_flag = 0;
	int		conc_flag = 0;
	int             status = 0;
	int		count_zero = 0;
	int		defined;		/* test on imgt and imgr voxel values */

	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (j = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'H': H_flag++;			break; 
				case 'Z': count_zero++;			break; 
				case 'o': strcpy (outfile, ptr);	*ptr = '\0'; break; 
				case 'L': norm = atoi (ptr);		*ptr = '\0'; break; 
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
	if (norm < 0 || norm > 2) {
		fprintf (stderr, "option -L aguments must be in the range 0-2\n");
		usage (program);
	}
	if (j < 2) usage (program);

	if (strlen (outfile)) {
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	} else {
		outfp = stdout;
	}
 	if (H_flag) write_command_line (outfp, argc, argv);

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
	nROI = ROIdim[3];
	imgr = calloc_double2 (dimension, nROI);
	fprintf (stdout, "Reading: %s\n", ROIfile);
	if (!(ROIfp = fopen (ROIfile, "rb"))) errr (program, ROIfile);
	for (j = 0; j < nROI; j++) {
		if (eread (imgt, dimension, isbigm, ROIfp)) errr (program, ROIfile);
		for (t = n = i = 0; i < dimension; i++) {
			defined = isnormal (imgt[i]) && imgt[i] != (float) 1.e-37;
			imgr[i][j] = q = (defined) ? imgt[i] : 0.;
			if (defined) switch (norm) {
				case 0: t += imgr[i][j] = 1.0; 	break;
				case 1: t += fabs (q);		break;
				case 2: t += q*q; n++;		break;
			}
		}
		if (norm == 2) t = sqrt(n*t);
		for (i = 0; i < dimension; i++) imgr[i][j] /= t;
	}
	if (fclose (ROIfp)) errr (program, ROIfile);

/***********/
/* process */
/***********/
	if (strlen (outfile)) fprintf (stdout, "Writing: %s\n", outfile);
	if (!(total = (double *) malloc (nROI*sizeof (double)))) errm (program);
	printf ("Reading: %s\n", imgfile);
	for (k = 0; k < imgdim[3]; k++) {
		if (conc_flag) {
			conc_read_vol (&conc_block, imgt);
		} else {
			if (eread (imgt, dimension, isbig, imgfp)) errr (program, imgfile);
		}
		for (j = 0; j < nROI; j++) total[j] = 0.0;
		for (i = 0; i < dimension; i++) {
			defined = isnormal (imgt[i]) && imgt[i] != (float) 1.e-37;
			if (count_zero && imgt[i] == 0.) defined = 1;
			if (!defined) continue;
			for (j = 0; j < nROI; j++) total[j] += imgt[i]*imgr[i][j];
		}
		for (j = 0; j < nROI; j++) fprintf (outfp, " %9.4f", total[j]);
		fprintf (outfp, "\n");
	}
	if (strlen (outfile)) if (fclose (outfp)) errw (program, outfile);

/*********************/
/* clean up and exit */
/*********************/
	if (conc_flag) {
		conc_free (&conc_block);
	} else {
		if (fclose (imgfp)) errr (program, imgfile);
	}
	free (total);
	free (imgt); free_double2 (imgr);
	exit (status);
}
