/*$Header: /data/petsun4/data1/src_solaris/sqrt_4dfp/RCS/sqrt_4dfp.c,v 1.7 2008/03/14 02:24:00 avi Exp $*/
/*$Log: sqrt_4dfp.c,v $
 * Revision 1.7  2008/03/14  02:24:00  avi
 * linux compliant
 *
 * Revision 1.6  2006/09/24  05:42:19  avi
 * Solaris 10
 *
 * Revision 1.5  2006/08/07  02:39:37  avi
 * safe 1.e-37 test
 *
 * Revision 1.4  2005/09/10  06:07:56  avi
 * correct handling of undefined voxels
 *
 * Revision 1.3  2005/01/22  23:04:45  avi
 * -E
 *
 * Revision 1.2  2005/01/22  22:46:26  avi
 * updated utilities
 *
 * Revision 1.1  1999/06/29  20:47:33  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <float.h>
#ifndef _FEATURES_H		/* defined by math.h in linux */
	#include <sunmath.h>	/* needed for isnormal() in Solaris but not linux */
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

static char rcsid[]= "$Id: sqrt_4dfp.c,v 1.7 2008/03/14 02:24:00 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	FILE            *fp_img, *fp_out;
	IFH		ifh;
	char            imgfile[MAXL];
	char            outfile[MAXL];
	char            imgroot[MAXL] = "";
	char            outroot[MAXL] = "";
	char            *str, command[MAXL], program[MAXL];

/**************/
/* processing */
/**************/
	int             imgdim[4], orient, isbig;
	float           voxdim[3];
	float           *imgt;
	float		amax = -FLT_MAX, amin = FLT_MAX;
	int             dimension, c, i, k;
	char		control = '\0';

/*********/
/* flags */
/*********/
	int             defined;
	int             status = 0;
	int		E_flag = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);

/************************/
/* process command line */
/************************/
        for (k = 0, i = 1; i < argc; i++) {
                if (*argv[i] == '-') {
		strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case 'E': E_flag++;		break;
				case '@': control = *str++;	*str++ = '\0'; break;
			}
		} else switch (k) {
                        case 0: getroot (argv[i], imgroot);	k++; break;
                        case 1: getroot (argv[i], outroot);	k++; break;
                }
        }
        if (k < 1) {
		fprintf (stderr, "Usage:\t%s <(4dfp) image> [outroot]\n", program);
		fprintf (stderr, "e.g.,\t%s vce20_mpr\n", program);
		fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default input endian)\n");
		fprintf (stderr, "\t-E\toutput undefined voxels as 1.0e-37 (default 0.0)\n");
		fprintf (stderr, "N.B.:\tdefault output filename = <image>_sqrt\t\n");
		exit (1);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
/***************************************/
/* create output filename if not given */
/***************************************/
	if (!strlen (outroot)) sprintf (outroot, "%s_sqrt", imgroot);	
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
			defined = isnormal (imgt[i]) && imgt[i] != (float) 1.e-37 && imgt[i] >= 0.0;
			if (E_flag) {
				if (imgt[i] == (float) 1.e-37) continue;
				imgt[i] = (defined) ? sqrt (imgt[i]) : 1.e-37;
			} else {
				imgt[i] = (defined) ? sqrt (imgt[i]) : 0.0;
			}
			if (defined && imgt[i] > amax) amax = imgt[i];
			if (defined && imgt[i] < amin) amin = imgt[i];
		}
		if (ewrite (imgt, dimension, control, fp_out)) errw (program, outfile);
	}
	fclose (fp_img);
	fclose (fp_out);

/*******/
/* ifh */
/*******/
	Writeifh (program, outfile, &ifh, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s -r%.0fto%.0f", outroot, amin, amax);
	printf ("%s\n", command);
	status |= system (command);

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, rcsid, control);
	catrec    (imgfile);
	endrec    ();

	free (imgt);
	exit (status);
}
