/*$Header: /data/petsun4/data1/src_solaris/sqrt_4dfp/RCS/t2z_4dfp.c,v 1.8 2009/01/14 03:34:32 avi Exp $*/
/*$Log: t2z_4dfp.c,v $
 * Revision 1.8  2009/01/14  03:34:32  avi
 * enable multi-volume operation
 *
 * Revision 1.7  2008/03/14  03:47:56  avi
 * linux compliant
 *
 * Revision 1.6  2006/09/25  00:42:25  avi
 * void usage ()
 *
 * Revision 1.5  2006/09/25  00:38:57  avi
 * Solaris 10
 *
 * Revision 1.4  2006/08/07  02:43:14  avi
 * safe 1.e-37 test
 *
 * Revision 1.3  2005/09/19  05:58:42  avi
 * output trailer for log image lnp -> log10p
 *
 * Revision 1.2  2005/09/19  05:43:11  avi
 * option -l (log10(prob(t))
 *
 * Revision 1.1  2005/09/19  00:46:23  avi
 * Initial revision
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <float.h>
#ifndef _FEATURES_H		/* defined by math.h in linux */
	#include <ieeefp.h>	/* needed for isnan() in Solaris but not linux */
#endif
#include <unistd.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MAXL 256

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

/********************************/
/* bertulani_betai.c prototypes */
/********************************/
double	betai (double a, double b, double x);
double	beta (double a, double b);
double	lbeta (double a, double b);

/************************/
/* statsub.c prototypes */
/************************/
double	lnp_z (double z);
double	p_z (double z);
double	dlnpdz (double z);
double	z_lnp (double lnp);
double	p_t (double t, double D);
double	lnp_t (double t, double D);
double	p_f (double f, double D1, double D2);
double	lnp_f (double f, double D1, double D2);

void usage (char *program) {
	fprintf (stderr, "Usage:\t%s <(4dfp) t-image>\n", program);
	fprintf (stderr, " e.g.,\t%s NP705_cond1_zfrm_RFX -nNP705_cond1_N\n", program);
	fprintf (stderr, "\toption\n");
	fprintf (stderr, "\t-l\toutput -log10(prob(t))\n");
	fprintf (stderr, "\t-N<flt>\tspecify global n\n");
	fprintf (stderr, "\t-n<str>\tspecify 4dfp n-image (up to two allowed)\n");
	fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default input endian)\n");
	fprintf (stderr, "N.B.:\tundefined (1.e-37, NaN) voxels in input are output as 1.e-37\n");
	fprintf (stderr, "N.B.:\toutput values are assigned the same sign as the input t value\n");
	fprintf (stderr, "N.B.:\tthe same n values apply to all volumes the input <t-image>\n");
	exit (1);
}

static char rcsid[]= "$Id: t2z_4dfp.c,v 1.8 2009/01/14 03:34:32 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	FILE		*fp_img, *fp_out;
	IFH		ifh;
	char		imgfile[MAXL], outfile[MAXL], imgroot[MAXL], outroot[MAXL];
	char		nmgroot[2][MAXL], nmgfile[2][MAXL];
	char		*str, command[MAXL], program[MAXL];

/**************/
/* processing */
/**************/
	char		control = '\0';
	int		imgdim[4], orient, imgdimn[4], orientn, isbig, isbign;
	float		voxdim[3], voxdimn[3];
	float		*imgt, *imgz, *imgn[2];
	int		dimension, c, i, k;
	int		nimgn = 0;		/* n-image count */
	double		an, nu, q;
	int		ix, iy, iz;

/*********/
/* flags */
/*********/
	int		status = 0;
	int		log_flag = 0;
	int		defined, sign;
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
				case 'd': debug++;				break;
				case 'l': log_flag++;				break;
				case '@': control = *str++;			*str = '\0'; break;
				case 'N': an = atof (str);			*str = '\0'; break;
				case 'n': if (nimgn == 2) usage (program);
					getroot (str, nmgroot[nimgn++]);	*str = '\0'; break;
			}
		} else switch (k) {
                        case 0: getroot (argv[i], imgroot);			k++; break;
                }
        }
        if (k < 1) usage (program);
	if (an < 2.0 && !nimgn) {
		fprintf (stderr, "%s: global N must be at least 2\n", program);
		usage (program);
	}

/*********************/
/* open 4dfp t-image */
/*********************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	if (Getifh (imgroot, &ifh)) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';
	dimension = imgdim[0] * imgdim[1] * imgdim[2];
	if (!(imgt = (float *) malloc (dimension * sizeof (float)))) errm (program);
	fprintf (stdout, "Reading: %s\n", imgfile);

/*********************/
/* get 4dfp n-images */
/*********************/
	for (k = 0; k < nimgn; k++) {
		sprintf (nmgfile[k], "%s.4dfp.img", nmgroot[k]);
		if (get_4dfp_dimoe (nmgfile[k], imgdimn, voxdimn, &orientn, &isbign) < 0)
				errr (program, nmgfile[k]);
		status = orient - orientn;
		for (i = 0; i < 3; i++) status |= (imgdimn[i] != imgdim[i]);
		for (i = 0; i < 3; i++) status |= (fabs (voxdimn[i] - voxdim[i]) > 0.0001);
		if (status) {
			fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot, nmgroot[k]);
			exit (-1);
		}
		fprintf (stdout, "Reading: %s\n", nmgfile[k]);
		if (!(imgn[k] = (float *) malloc (dimension * sizeof (float)))) errm (program);
		if (!(fp_img = fopen (nmgfile[k], "rb"))
		|| eread  (imgn[k], dimension, isbign, fp_img)
		|| fclose (fp_img)) errr (program, imgfile);
	}

/***********/
/* process */
/***********/
	sprintf (outroot, "%s_%s", imgroot, (log_flag) ? "log10p" : "z");	
	sprintf (outfile, "%s.4dfp.img", outroot);
	fprintf (stdout, "Writing: %s\n", outfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);

	if (!(imgz = (float *) malloc (dimension * sizeof (float)))) errm (program);
	for (k = 0; k < imgdim[3]; k++) {	/* loop on volumes */
		if (eread  (imgt, dimension, isbig, fp_img)) errr (program, imgfile);
		for (i = 0; i < dimension; i++) {
			defined = !isnan (imgt[i]) && imgt[i] != (float) 1.e-37;
			switch (nimgn) {
				case 0: nu = an - 1.0;				break;
				case 1: nu = imgn[0][i] - 1.0;			break;
				case 2: nu = imgn[0][i] + imgn[1][i] - 2.0;	break;
			}
			defined &= (nu >= 1.0);
			if (defined) {
				sign = (imgt[i] > 0.0) ? 1 : -1;
				if (log_flag) {
					imgz[i] = sign * (-1./M_LN10) * lnp_t (fabs ((double) imgt[i]), nu);
				} else {
					imgz[i] = sign * z_lnp (lnp_t (fabs ((double) imgt[i]), nu));
				}
			} else {
				imgz[i] = 1.e-37;
			}
		}
		if (ewrite (imgz, dimension, control, fp_out)) errw (program, outfile);
	}
	if (fclose (fp_img)) errr (program, imgfile);
	if (fclose (fp_out)) errw (program, outfile);

/*******/
/* ifh */
/*******/
	if (Writeifh (program, outfile, &ifh, control)) errw (program, outfile);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s", outroot);
	status |= system (command);

/*******/
/* rec */
/*******/
	startrece (outfile, argc, argv, rcsid, control);
	catrec (imgfile);
	for (k = 0; k < nimgn; k++) catrec (nmgfile[k]);
	endrec ();

	for (k = 0; k < nimgn; k++) free (imgn[k]);
	free (imgt); free (imgz);
	exit (status);
}
