/*$Header: /data/petsun4/data1/src_solaris/zero_lt_4dfp/RCS/zero_gt_4dfp.c,v 1.6 2007/08/18 23:03:11 avi Exp $*/
/*$Log: zero_gt_4dfp.c,v $
 * Revision 1.6  2007/08/18  23:03:11  avi
 * typo
 *
 * Revision 1.5  2007/08/18  23:00:50  avi
 * endian compliant
 *
 * Revision 1.4  2005/10/25  03:27:51  avi
 * allow first argument to be negative (not pseudo-option)
 *
 * Revision 1.3  2004/11/05  21:29:33  rsachs
 * Removed 'Get4dfpDinN'. Installed 'errm','errr','errw','getroot','setprog','get_4dfp_dimo'.
 *
 * Revision 1.2  2001/08/02  01:15:08  avi
 * correct usage
 *
 * Revision 1.1  2001/08/02  00:52:07  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
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

static char rcsid[]= "$Id: zero_gt_4dfp.c,v 1.6 2007/08/18 23:03:11 avi Exp $";
int main (int argc, char **argv) {
/*************/
/* image I/O */
/*************/
	FILE            *fp_img, *fp_out;
	IFH		ifh;
	char            imgfile[MAXL], imgroot[MAXL];
	char            outfile[MAXL], outroot[MAXL] = "";
	char		control = '\0';

/**************/
/* processing */
/**************/
	int             imgdim[4], dimension, orient, isbig;
	float           voxdim[3];
	float           *imgr;
	float		thresh;

/***********/
/* utility */
/***********/
	int             c, i, k;
	char            *ptr, command[MAXL], program[MAXL];

/*********/
/* flags */
/*********/
        int             status = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (i > 1 && *argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case '@': control = *ptr++;	*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0: thresh = atof (argv[i]);    k++;	break;
			case 1: getroot (argv[i], imgroot); k++;	break;
			case 2: getroot (argv[i], outroot); k++;	break;
		}
	}
	if (k < 2) {
		printf ("Usage:\t%s <flt> <(4dfp) image> [outroot] [options]\n", program);
		printf (" e.g.,\t%s 90 pt349_study9to9\n", program);
		printf (" e.g.,\t%s 90 pt349_study9to9 pt349_study9to9z\n", program);
		printf ("\toptions\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:\tdefault output 4dfp root is <(4dfp) image>\"z\"\n");
		printf ("N.B.:\tfirst field can't be used for options because threshold might be negative\n");
		exit (1);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
/***************************************/
/* create output filename if not given */
/***************************************/
	if (!strlen (outroot)) sprintf (outroot, "%sz", imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);

/*****************************/
/* get 4dfp input dimensions */
/*****************************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig)
	||  Getifh (imgfile, &ifh)) errr (program, imgfile);
	dimension = imgdim[0] * imgdim[1] * imgdim[2];
	if (!control) control = (isbig) ? 'b' : 'l';

/*****************/
/* alloc buffers */
/*****************/
	if (!(imgr = (float *) malloc (dimension * sizeof (float)))) errm (program);

/***********/
/* process */
/***********/
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	fprintf (stdout, "Reading: %s\n", imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);
	for (k = 0; k < imgdim[3]; k++) {
		if (eread  (imgr, dimension, isbig, fp_img)) errr (program, imgfile);
		for (i = 0; i < dimension; i++) if (imgr[i] > thresh) imgr[i] = 0.0;
		if (ewrite (imgr, dimension, control, fp_out)) errw (program, outfile);
	}
	fclose (fp_img);
	fclose (fp_out);

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, rcsid, control);
	catrec   (imgfile);
        endrec   ();

/*********************/
/* ifh and hdr files */
/*********************/
	Writeifh (program, outfile, &ifh, control);
	sprintf (command, "ifh2hdr %s", outroot);
	status |= system (command);

        free (imgr);
	exit (status);
}
