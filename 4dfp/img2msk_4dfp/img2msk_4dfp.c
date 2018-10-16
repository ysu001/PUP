/*$Id: img2msk_4dfp.c,v 1.7 2007/03/02 07:38:37 avi Exp $*/
/*$Log: img2msk_4dfp.c,v $
 * Revision 1.7  2007/03/02  07:38:37  avi
 * Solaris 10
 *
 * Revision 1.6  2004/09/17  20:21:43  rsachs
 * Cosmetic minutiae.
 *
 * Revision 1.5  2004/09/16  21:07:45  rsachs
 * Replaced 'Get4dfpDimN' with 'get_4dfp_dimo'.
 *
 * Revision 1.4  1998/12/22  00:02:49  avi
 * correct default outroot
 *
 * Revision 1.3  1998/12/21  23:07:28  avi
 * default threshold = 1.0
 *
 * Revision 1.2  1998/12/21  23:03:47  avi
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MAXL 256

/*************/
/* externals */
/*************/
extern void	fimg2msk_ (float *thresh, float *imag, int *nx, int *ny, int *nz, float *mask);	/* fimg2msk.f */

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[] = "$Id: img2msk_4dfp.c,v 1.7 2007/03/02 07:38:37 avi Exp $";
int main (int argc, char **argv) {
/*************/
/* image I/O */
/*************/
        FILE            *fpimg, *fpout;
	IFH		ifh;
        char            imgfile[MAXL], outfile[MAXL];
        char            imgroot[MAXL], outroot[MAXL] = "";
        char            *ptr, command[MAXL], program[MAXL];

/**************/
/* processing */
/**************/
        int             imgdim[4], dimension, orient, isbig;
        float           voxdim[3];
        float           *img3d, *msk3d;
	float		thresh = 1.0;
	char		control = '\0';

/***********/
/* utility */
/***********/
        int             c, i, k;

/*********/
/* flags */
/*********/
        int             status = 0;

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 't': thresh = atof (argv[i]);	*ptr = '\0'; break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
			}
		}
		else switch (k) {
                	case 0: getroot (argv[i], imgroot); k++; break;
                	case 1: getroot (argv[i], outroot); k++; break;
		}
	}
	if (k < 1) {
		printf ("Usage:\timg2msk_4dfp [options] <imgfile> [outfile]\n");
        	printf (" e.g.:\timg2msk_4dfp va1234_mpr -t12.3 va1234_mpr_msk\n");
        	printf ("\toption\n");
        	printf ("\t-t\tspecify threshold (default = 1.0)\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		exit (1);
	}
	sprintf (imgfile, "%s.4dfp.img", imgroot);

	if (!strlen (outroot)) sprintf (outroot, "%s_msk", imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);

/*****************/
/* alloc buffers */
/*****************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';
	dimension = imgdim[0] * imgdim[1] * imgdim[2];
	img3d = (float *) malloc (dimension * sizeof (float));
	msk3d = (float *) malloc (dimension * sizeof (float));
	if (!img3d || !msk3d) errm (program);

/***********/
/* process */
/***********/
	if (!(fpimg = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!(fpout = fopen (outfile, "wb"))) errw (program, outfile);
	fprintf (stdout, "Reading: %s\n", imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);
	for (k = 0; k < imgdim[3]; k++) {
		if (eread (img3d, dimension, isbig, fpimg))	errr (program, imgfile);
		for (i = 0; i < dimension; i++) msk3d[i] = 1.;
		fimg2msk_ (&thresh, img3d, imgdim + 0, imgdim + 1, imgdim + 2, msk3d);
		if (ewrite (img3d, dimension, control, fpout))	errm (program);
	}
	fclose (fpimg);
	fclose (fpout);

/***************/
/* ifh and hdr */
/***************/
	if (Getifh (imgroot, &ifh)) errr (program, imgroot);
	if (Writeifh (program, outfile, &ifh, control)) errw (program, outfile);
	sprintf (command, "ifh2hdr %s", outroot); printf ("%s\n", command);
	status = system (command);

/*******************/
/* create rec file */
/*******************/
	startrece	(outfile, argc, argv, rcsid, control);
	catrec		(imgfile);
	endrec		();

	free (img3d);
	free (msk3d);
	exit (status);
}
