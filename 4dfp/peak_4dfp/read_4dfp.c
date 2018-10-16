/*$Header: /data/petsun4/data1/src_solaris/peak_4dfp/RCS/read_4dfp.c,v 1.4 2010/08/18 04:22:23 avi Exp $*/
/*$Log: read_4dfp.c,v $
 * Revision 1.4  2010/08/18  04:22:23  avi
 * correct endian dependent error
 *
 * Revision 1.3  2007/09/23  03:01:27  avi
 * Solaris 10 and linux compliant
 *
 * Revision 1.2  2006/03/07  04:25:11  avi
 * always output 6 decimal places with cr
 *
 * Revision 1.1  2006/03/03  03:09:53  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <Getifh.h>
#include <endianio.h>

#define MAXL		256

/*************/
/* externals */
/*************/
void vrtflip	(int *iori, int *imgdim, float *centeri, float *mmppixi, float *centert, float *mmppixt);	/* cvrtflip.c */
void imgvalx_	(float *imgt, int *nx, int *ny, int *nz, float *center, float *mmppix,
									float *x, float *v, int *lslice);	/* imgvalx.f */

/********************/
/* global variables */
/********************/
static	char		program[MAXL];
static	char		rcsid[] = "$Id: read_4dfp.c,v 1.4 2010/08/18 04:22:23 avi Exp $";

int main (int argc, char *argv[]) {
	FILE		*imgfp;			/* input file */

/*************/
/* image i/o */
/*************/
	char		*ptr, command[MAXL];
	char		imgroot[MAXL], imgfile[MAXL], ifhfile[MAXL]; 
	int		isbig;

/***************/
/* computation */
/***************/
	float		*imgv;
	IFH		imgifh;
	float		mmppixr[3], centerr[3];
	float		x0[3], value;
	int		dim;
	int		c, i, j, k;

/*********/
/* flags */
/*********/
	int		status = 0;
	int		verbose = 0;
	int		debug = 0;
	
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);

/******************************/
/* get command line arguments */
/******************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-' && k > 2) {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'v': verbose++;	break;
			}
		} else switch (k) {
			case 0: case 1: case 2: x0[k] = atof (argv[i]);	k++; break;
		 	case 3: getroot (argv[i], imgroot);		k++; break;
		}	
	}
	if (k < 4) {
		printf ("%s\n", rcsid);
		printf ("Usage:	%s <flt x0> <flt y0> <flt z0> <4dfp imgroot> [options]\n", program);
		printf ("e.g.,	%s 33.1 -56.2 18. grand_average_222[.4dfp.img]\n", program);
		printf ("\toption\n");
		printf ("\t-v	verbose mode\n");
		exit (1);
	}

/************/
/* read ifh */
/************/
	sprintf (ifhfile, "%s.4dfp.img", imgroot);
	if (verbose) printf ("Reading: %s\n", ifhfile);
	if (Getifh (ifhfile, &imgifh)) errr (program, ifhfile);
	isbig = strcmp (imgifh.imagedata_byte_order, "littleendian");
	if (verbose) printf ("image dimensions \t%10d%10d%10d%10d\n",
			imgifh.matrix_size[0], imgifh.matrix_size[1], imgifh.matrix_size[2], imgifh.matrix_size[3]);
	if (verbose) printf ("image mmppix     \t%10.6f%10.6f%10.6f\n",
			imgifh.mmppix[0], imgifh.mmppix[1], imgifh.mmppix[2]);
	if (verbose) printf ("image center     \t%10.4f%10.4f%10.4f\n",
			imgifh.center[0], imgifh.center[1], imgifh.center[2]);
	if (verbose) printf ("image orientation\t%10d\n", imgifh.orientation);
/*****************************************/
/* virtual flip instead of x4dfp2ecat () */
/*****************************************/
	vrtflip (&imgifh.orientation, imgifh.matrix_size, imgifh.center, imgifh.mmppix, centerr, mmppixr);

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	dim = 1; for (k = 0; k < 3; k++) dim *= imgifh.matrix_size[k];
	if (!(imgv = (float *) calloc (dim, sizeof (float)))) errm (program);
	if (verbose) printf ("Reading: %s\n", imgfile);
	if (!(imgfp = fopen (imgfile, "rb")) || eread (imgv, dim, isbig, imgfp)
	|| fclose (imgfp)) errr (program, imgfile);

	imgvalx_ (imgv, imgifh.matrix_size+0, imgifh.matrix_size+1, imgifh.matrix_size+2, centerr, mmppixr, x0, &value, &k);
	if (verbose) {
		printf ("value=%.6f on slice %d\n", value, k);
	} else {
		printf ("%.6f\n", value);
	}

	free (imgv);
	exit (status);
}
