/*$Id: endian_4dfp.c,v 1.5 2007/04/03 03:44:59 avi Exp $*/
/*$Log: endian_4dfp.c,v $
 * Revision 1.5  2007/04/03  03:44:59  avi
 * Linux gcc compatible
 *
 * Revision 1.4  2007/02/28  06:52:46  avi
 * Solaris 10
 * create hdr using ifh2hdr
 *
 * Revision 1.3  2006/08/07  03:06:52  avi
 * correct 1.e-37 test
 *
 * Revision 1.2  2006/04/01  06:50:18  avi
 * option -t
 *
 * Revision 1.1  2006/03/26  01:33:36  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#ifndef fpclassify
	#include <sunmath.h>
	#include <ieeefp.h>
#endif
#include <rec.h>
#include <Getifh.h>
#include <endianio.h>

#define MAXL 	256

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[] = "$Id: endian_4dfp.c,v 1.5 2007/04/03 03:44:59 avi Exp $";
int main (int argc, char *argv[]) {
	FILE		*imgfp;				/* input and output */
	IFH		ifh;
	char		imgroot[MAXL], imgfile[MAXL];	/* input 4dfp file name */
	char		recfile[MAXL], tmpfile[MAXL];
	char		control = '\0';			/* output endian */

/***************/
/* 4dfp arrays */
/***************/
	float		*imgt;				/* frame buffer */
	float		voxdim[3];			/* voxel dimensions in mm */
	int		imgdim[4], vdim, orient, isbig, osbig;

/***********/
/* utility */
/***********/
	int		c, i, j, k, n;
	char		*ptr, command[MAXL], program[MAXL];
	double		q, r, s, t;

/*********/
/* flags */
/*********/
	int		defined, swab_flag;
	int		varlog_flag = 0;
	int		status = 0;

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case '@': control = *ptr++;	*ptr = '\0'; break;
				case 't': varlog_flag++;	break;
			}
		}
		else switch (k) {
			case 0:	getroot (argv[i], imgroot); k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:	%s <(4dfp) image>\n", program);
		printf ("\toption\n");
		printf ("\t-@<b|l|c> make <(4dfp) image> big, little or CPU endian\n");
		printf ("\t-t	perform var(log(fabs(.))) test\n");
		printf ("N.B.:	<(4dfp) image> may be overwritten\n", program);
		printf ("N.B.:	absent option -@ %s only reports state of <(4dfp) image>\n", program);
		exit (1);
	}

/***********************************/
/* get input 4dfp image dimensions */
/***********************************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig)) errr (program, imgfile);
	vdim = imgdim[0]*imgdim[1]*imgdim[2];

	if (varlog_flag) {
/**********************************************/
/* compute variance of image value logarithms */
/**********************************************/
		if (!(imgt = (float *) calloc (vdim, sizeof (float)))) errm (program);

		fprintf (stdout, "computing var(log(fabs(.))) test %s frame", imgfile);
		if (!(imgfp = fopen (imgfile, "r"))) errw (program, imgfile);
		s = t = n = 0;
		for (i = 0; i < imgdim[3]; i++) {printf (" %d", i + 1); fflush (stdout);
			if (gread ((char *) imgt, sizeof (float), vdim, imgfp, isbig)) errr (program, imgfile);
			for (j = 0; j < vdim; j++) {
				q = imgt[j];
				defined = imgt[j] != (float) 1.e-37 && isnormal (q);
				if (!defined) continue;
				r = log (fabs (q));
				s += r;
				t += r*r;
				n++;
			}
		}
		fclose (imgfp); free (imgt);
		s /= n; t /= n; t -= s*s;
		printf ("\nmean(log(fabs(.)))=%f var(log(fabs(.)))=%f nvox=%d\n", s, t, n);
	}
	if (!control) exit (0);

	osbig = (CPU_is_bigendian ()) ? !(control == 'l' || control == 'L') : (control == 'b' || control == 'B');
	swab_flag = (osbig != 0) != (isbig != 0);
	if (!swab_flag) {
		exit (0);
	} else {
		vdim = imgdim[0]*imgdim[1]*imgdim[2];
		if (!(imgt = (float *) calloc (vdim, sizeof (float)))) errm (program);

		fprintf (stdout, "swapping endian %s frame", imgfile);
		if (!(imgfp = fopen (imgfile, "r+b"))) errw (program, imgfile);
/**********************/
/* process all frames */
/**********************/
		for (i = 0; i < imgdim[3]; i++) {printf (" %d", i + 1); fflush (stdout);
			if (fseek (imgfp, (long) i*vdim*sizeof (float), SEEK_SET)
			||  gread ((char *) imgt, sizeof (float), vdim, imgfp, isbig))		errr (program, imgfile);
			if (fseek (imgfp, (long) i*vdim*sizeof (float), SEEK_SET)
			|| gwrite ((char *) imgt, sizeof (float), vdim, imgfp, control))	errw (program, imgfile);
		}
		printf ("\n");
		free (imgt);
	}

/*******/
/* ifh */
/*******/
	Getifh (imgfile, &ifh);
	Writeifh (program, imgfile, &ifh, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s", imgroot);
	printf ("%s\n", command);
	status |= system (command);

/*******/
/* rec */
/*******/
	sprintf (recfile, "%s.rec", imgfile);
	sprintf (tmpfile, "%s0",    recfile);
	if (rename (recfile, tmpfile)) errw (program, tmpfile);
	startrece (imgfile, argc, argv, rcsid, control);
	catrec    (tmpfile);
	endrec    ();
	status |= remove (tmpfile);

	exit (status);
}

