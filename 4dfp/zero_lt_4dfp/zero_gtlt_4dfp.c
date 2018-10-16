/*$Header: /data/petsun4/data1/src_solaris/zero_lt_4dfp/RCS/zero_gtlt_4dfp.c,v 1.4 2007/08/19 00:10:26 avi Exp $*/
/*$Log: zero_gtlt_4dfp.c,v $
 * Revision 1.4  2007/08/19  00:10:26  avi
 * endian compliant
 *
 * Revision 1.3  2004/11/08  20:06:55  rsachs
 * Installed 'errm','errr','errw','getroot','setprog','get_4dfp_dimo'. Removed 'Get4dfpDimN'.
 *
 * Revision 1.2  2002/03/25  22:33:47  avi
 * better command line argument recovery
 *
 * Revision 1.1  2002/03/25  21:52:33  avi
 * Initial revision
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>

#define MAXL 256

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

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[]= "$Id: zero_gtlt_4dfp.c,v 1.4 2007/08/19 00:10:26 avi Exp $";
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
	float		tlo = 0.0, thi = 0.0;

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
		strcpy (command, argv[i]); ptr = command;
		if (*argv[i] == '-' && k) {
			while (c = *ptr++) switch (c) {
				case '@': control = *ptr++;	*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0: getrange (command, &tlo, &thi); k++; break;
		 	case 1: getroot  (command, imgroot);	k++; break;
		 	case 2: getroot  (command, outroot);	k++; break;
		}
	}
	if (k < 2) {
		printf ("Usage:\t%s <flt[to<flt>]> <(4dfp) image> [outroot] [options]\n", program);
		printf (" e.g.,\t%s -30to90 pt349_study9to9\n", program);
		printf ("\toptions\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:\t%s zeroes voxel values within the specified range\n", program);
		printf ("N.B.:\tdefault output 4dfp root is <(4dfp) image>\"z\"\n");
		printf ("N.B.:\tfirst field can't be used for options because lower range might be negative\n");
		exit (1);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (!strlen (outroot)) sprintf (outroot, "%sz", imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);
/*****************************/
/* get 4dfp input dimensions */
/*****************************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig)
	||  Getifh (imgfile, &ifh)) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';

	dimension = imgdim[0] * imgdim[1] * imgdim[2];
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
		for (i = 0; i < dimension; i++) if (imgr[i] > tlo && imgr[i] < thi) imgr[i] = 0.0;
		if (ewrite (imgr, dimension, control, fp_out)) errw (program, outfile);
	}
	fclose (fp_img);
	fclose (fp_out);

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, rcsid, control);
	sprintf (command, "Zeroed range %.4f to %.4f\n", tlo, thi); printrec (command);
	catrec  (imgfile);
        endrec  ();

/*********************/
/* ifh and hdr files */
/*********************/
	Writeifh (program, outfile, &ifh, control);
	sprintf (command, "ifh2hdr %s", outroot);
	status |= system (command);

        free (imgr);
	exit (status);
}

