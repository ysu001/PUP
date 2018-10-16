/*$Header: /data/petsun4/data1/src_solaris/sqrt_4dfp/RCS/z2logp_4dfp.c,v 1.6 2010/01/31 07:48:35 avi Exp $*/
/*$Log: z2logp_4dfp.c,v $
 * Revision 1.6  2010/01/31  07:48:35  avi
 * option -p
 * correct error in compuation of two-sided -log10(p) that inflated results by 0.994177
 *
 * Revision 1.5  2009/08/03  23:34:36  avi
 * correct several bugs including failure to parse command line '-2'
 * accept multi-volume input
 *
 * Revision 1.4  2008/03/14  04:00:54  avi
 * linux compliant
 *
 * Revision 1.3  2006/09/25  01:22:14  avi
 * Solaris 10
 *
 * Revision 1.2  2006/08/07  02:47:03  avi
 * safe 1.e-37 test
 *
 * Revision 1.1  2005/10/02  05:51:16  avi
 * Initial revision
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#ifndef _FEATURES_H		/* defined by math.h in linux */
	#include <ieeefp.h>	/* needed for isnan() in Solaris but not linux */
#endif
#include <unistd.h>
#include <Getifh.h>
#include <endianio.h>
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
	fprintf (stderr, "Usage:\t%s <(4dfp) Z-image>\n", program);
	fprintf (stderr, " e.g.,\t%s vce20_z[.4dfp[.img]]\n", program);
	fprintf (stderr, "\toption\n");
	fprintf (stderr, "\t-2\ttwo sided test (default one sided)\n");
	fprintf (stderr, "\t-p\toutput p-values (default output -log10(p))\n");
	fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default input endian)\n");
	fprintf (stderr, "N.B.:\tprobability computed on assumption that voxel values are N(0,1)\n");
	fprintf (stderr, "N.B.:\tundefined (1.e-37, NaN, Inf) voxels in input are output as 1.e-37\n");
	exit (1);
}

static char rcsid[]= "$Id: z2logp_4dfp.c,v 1.6 2010/01/31 07:48:35 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	FILE		*fp_img, *fp_out;
	IFH		ifh;
	char		imgfile[MAXL], outfile[MAXL], imgroot[MAXL], outroot[MAXL];
	char		*str, command[MAXL], program[MAXL];

/**************/
/* processing */
/**************/
	char		control = '\0';
	int		dimension, c, i, j, k;
	int		imgdim[4], orient, isbig;
	float		voxdim[3], voxdimn[3];
	float		*imgz, *imgp;
	double		q;

/*********/
/* flags */
/*********/
	int		status = 0;
	int		p_flag = 0;
	int		two_sided = 0;
	int		defined;
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
				case 'p': p_flag++;		break;
				case '2': two_sided++;		break;
				case '@': control = *str++;	*str = '\0'; break;
			}
		} else switch (k) {
                        case 0: getroot (argv[i], imgroot);	k++; break;
                }
        }
        if (k < 1) usage (program);

/**************************/
/* get input 4dfp z-image */
/**************************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	if (Getifh (imgfile, &ifh)) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';
	dimension = imgdim[0] * imgdim[1] * imgdim[2];
	if (!(imgz = (float *) malloc (dimension * sizeof (float)))) errm (program);
	if (!(imgp = (float *) malloc (dimension * sizeof (float)))) errm (program);

/***************************/
/* prepare 4dfp read/write */
/***************************/
	fprintf (stdout, "Reading: %s\n", imgfile);
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	sprintf (outroot, "%s_%s", imgroot, (p_flag) ? "p" : "log10p");	
	sprintf (outfile, "%s.4dfp.img", outroot);
	fprintf (stdout, "Writing: %s\n", outfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);

/***********/
/* process */
/***********/
	for (j = 0; j < imgdim[3]; j++) {
		if (eread (imgz, dimension, isbig, fp_img)) errr (program, imgfile);
		for (i = 0; i < dimension; i++) {
			q = (double) imgz[i];
			defined = !isnan (q) && imgz[i] != (float) 1.e-37 && finite (q);
			if (defined) {
				if (two_sided) {
					q = fabs (q);
					imgp[i] =         (p_flag) ? 2.*p_z (q) : (-1./M_LN10)*(lnp_z (q) + M_LN2);
				} else {
					if (q > 0) {
						imgp[i] = (p_flag) ?    p_z (q) : (-1./M_LN10)*lnp_z (q);
					} else {
						imgp[i] = (p_flag) ?        1.0 : 0.0;
					}
				}
			} else {
				imgp[i] = 1.e-37;
			}
		}
		if (ewrite (imgp, dimension, control, fp_out)) errw (program, outfile);
	}
	if (fclose (fp_img)) errr (program, imgfile);
	if (fclose (fp_out)) errw (program, outfile);

/*******/
/* ifh */
/*******/
	if (Writeifh (program, outfile, &ifh, control)) errw (program, outroot);

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
	endrec ();

	free (imgz); free (imgp);
	exit (status);
}
