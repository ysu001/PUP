/*$Header: /data/petsun4/data1/src_solaris/sqrt_4dfp/RCS/rho2z_4dfp.c,v 1.6 2010/12/23 03:33:45 avi Exp $*/
/*$Log: rho2z_4dfp.c,v $
 * Revision 1.6  2010/12/23  03:33:45  avi
 * long -> int in rnan()
 *
 * Revision 1.5  2007/11/20  22:40:39  avi
 * linux compliant (locally code equivalent of quiet_nan())
 *
 * Revision 1.4  2006/09/24  23:37:24  avi
 * Solaris 10
 *
 * Revision 1.3  2006/08/07  02:31:51  avi
 * make 1.e-37 test safe
 *
 * Revision 1.2  2005/07/08  02:15:42  avi
 * -r (reverse operation) option
 *
 * Revision 1.1  2005/01/23  00:41:50  avi
 * Initial revision
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <unistd.h>
#ifndef _FEATURES_H		/* defined by math.h in linux */
	#include <ieeefp.h>	/* needed for isnan() in Solaris but not linux */
#endif
#include <rec.h>
#include <Getifh.h>
#include <endianio.h>

#define MAXL 256

float rnan (void) {
	union {
		float	r;
		int	l;
	} nan;
	nan.l = 0x7fffffff;
	return nan.r;
}

static double z (double rho) {
	return (fabs (rho) < 1.) ? 0.5*log ((1. + rho)/(1. - rho)) : rnan ();
}

static void numerical_test () {
	int		i;
	double		q;

	for (i = 0; i < 15; i++) {
		q = (double) random() / (double) 0x7fffffffL;
		q = 3.0*q - 1.0;
		printf ("%5d %10.6f %10.6f %10.6f %10.6f\n", i, q, z(q), atanh(q), tanh(z(q)));
	}
	exit (0);
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[]= "$Id: rho2z_4dfp.c,v 1.6 2010/12/23 03:33:45 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	FILE		*fp_img, *fp_out;
	IFH		ifh;
	char		imgfile[MAXL], outfile[MAXL], imgroot[MAXL] = "", outroot[MAXL] = "";
	char		*str, command[MAXL], program[MAXL];
	char		trailerz[] = "zfrm", trailerr[] = "corr";

/**************/
/* processing */
/**************/
	char		control = '\0';
	int		imgdim[4], orient, isbig;
	float		voxdim[3];
	float		*imgt;
	int		dimension, c, i, k;
	double		q;

/*********/
/* flags */
/*********/
	int		z2r_flag = 0;
	int		status = 0;
	int		E_flag = 0;
	int		debug = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);

/************************/
/* process command line */
/************************/
        for (k = 0, i = 1; i < argc; i++) {
                if (*argv[i] == '-') {
		strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case 'd': debug++;		break;
				case 'r': z2r_flag++;		break;
				case 'E': E_flag++;		break;
				case '@': control = *str++;	*str = '\0'; break;
			}
		} else switch (k) {
                        case 0: getroot (argv[i], imgroot);	k++; break;
                        case 1: getroot (argv[i], outroot);	k++; break;
                }
        }
	if (debug) numerical_test ();
        if (k < 1) {
		fprintf (stderr, "Usage:\t%s <(4dfp) image> [outroot]\n", program);
		fprintf (stderr, "e.g.,\t%s vce20_rho[.4dfp[.img]]\n", program);
		fprintf (stderr, "\toption\n");
		fprintf (stderr, "\t-r\treverse (convert z to r)\n");
		fprintf (stderr, "\t-E\toutput undefined voxels as 1.0e-37 (default 0.0)\n");
		fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default input-endian)\n");
		fprintf (stderr, "N.B.:\tdefault r to z output filename root = <image>_%s\t\n", trailerz);
		fprintf (stderr, "N.B.:\tdefault z to r output filename root = <image>_%s\t\n", trailerr);
		exit (1);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
/***************************************/
/* create output filename if not given */
/***************************************/
	if (!strlen (outroot)) sprintf (outroot, "%s_%s", imgroot, (z2r_flag) ? trailerr : trailerz);	
	sprintf (outfile, "%s.4dfp.img", outroot);

/*****************************/
/* get 4dfp input dimensions */
/*****************************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	if (Getifh (imgfile, &ifh)) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';
	dimension = imgdim[0] * imgdim[1] * imgdim[2];

/*****************/
/* alloc buffers */
/*****************/
	if (!(imgt = (float *) malloc (dimension * sizeof (float)))) errm (program);

/***********/
/* process */
/***********/
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	fprintf (stdout, "Reading: %s\n", imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);
	for (k = 0; k < imgdim[3]; k++) {
		if (eread  (imgt, dimension, isbig, fp_img)) errr (program, imgfile);
		for (i = 0; i < dimension; i++) {
			q = (z2r_flag) ? tanh ((double) imgt[i]) : atanh ((double) imgt[i]);
			if (E_flag) {
				if (imgt[i] == (float) 1.e-37) continue;
				imgt[i] = (isnan (q)) ? 1.e-37 : q;
			} else {
				imgt[i] = (isnan (q)) ?    0.0 : q;
			}
		}
		if (ewrite (imgt, dimension, control, fp_out)) errw (program, outfile);
	}
	fclose (fp_img);
	fclose (fp_out);

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, rcsid, control);
	catrec    (imgfile);
	endrec    ();

/*******/
/* ifh */
/*******/
	if (Writeifh (program, outfile, &ifh, control)) errw (program, outfile);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s", outroot);
	status |= system (command);

	free (imgt);
	exit (status);
}
