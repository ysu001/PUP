/*$Header: /data/petsun4/data1/src_solaris/maskimg_4dfp/RCS/maskimg_4dfp.c,v 1.17 2012/06/27 22:33:19 avi Exp $*/
/*$Log: maskimg_4dfp.c,v $
 * Revision 1.17  2012/06/27  22:33:19  avi
 * MAXL->512 to accommodate long filenames
 *
 * Revision 1.16  2010/02/11  06:20:42  avi
 * assign option -R to suppress vreation of rec file
 * option -N now codes NaN operation
 *
 * Revision 1.15  2006/09/25  02:33:24  avi
 * Solaris 10
 *
 * Revision 1.14  2006/03/23  05:00:36  avi
 * variable endian competence
 *
 * Revision 1.13  2006/03/16  07:42:51  avi
 * option -A (read mask image values as absolute value)
 *
 * Revision 1.12  2004/09/16  20:16:28  rsachs
 * Replaced 'Get4dfpDimN' with 'get_4dfp_dimo'.
 *
 * Revision 1.11  2003/10/05  21:09:40  avi
 * -1 option (use only first frame opf mask)
 *
 * Revision 1.10  2002/03/22  02:14:18  avi
 * correct -e usage
 *
 * Revision 1.9  2001/03/13  22:38:28  avi
 * -e option also reports within-mask voxel count
 *
 * Revision 1.8  2001/03/13  05:24:03  avi
 * -e option (report mean within_mask value)
 *
 * Revision 1.7  1999/06/28  23:15:02  avi
 * -p option
 *
 * Revision 1.6  1999/02/21  08:17:26  avi
 * correct error in -v and -t evaluation
 *
 * Revision 1.5  1999/02/11  00:27:41  avi
 * -R
 * Revision 1.4  1998/12/19  07:05:10  avi
 * several revisions
 * Revision 1.3  1998/10/12  20:36:29  mcavoy
 * Revision 1.1  1997/09/30  21:00:19  tscull
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MAXL 512

static char rcsid[]= "$Id: maskimg_4dfp.c,v 1.17 2012/06/27 22:33:19 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	FILE            *imgfp[3];
	IFH		ifh;		/* only for first input image */
	char            imgroot[3][MAXL], imgfile[3][MAXL];

/**************/
/* processing */
/**************/
	int             imgdim[3][4], isbig[3], orient[3], dimension, nvox;
	float           voxdim[3][3];
	float           *imgt[2];
	float		threshold = 0.0, value, pct, imgmax;
	char		control = '\0';

/***********/
/* utility */
/***********/
	char            *ptr, command[MAXL], program[MAXL];
	int             c, i, j, k;
	double		q;

/*********/
/* flags */
/*********/
	int             test = 0;
	int             getmean = 0;
	int             uniform = 0;
	int             pctflag = 0;
	int             absflag = 0;
	int             onemask = 0;
	int             status = 0;
	int             opr = 0;
	int		saverec = 1;

	fprintf (stdout, "%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);

/************************/
/* process command line */
/************************/
	for (j = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'e': getmean++;				break;
				case '1': onemask++;				break;
				case 'N': opr = c;				break;
				case 'A': absflag++;				break;
				case 'R': saverec = 0;				break;
				case '@': control = *ptr++;			*ptr = '\0'; break;
				case 'v': uniform++; value = atof (ptr);	*ptr = '\0'; break;
				case 'p': pctflag++; pct = atof (ptr);		*ptr = '\0'; break;
				case 't': threshold = atof (ptr);		*ptr = '\0'; break;
			}
		} else {
                	getroot (argv[i], imgroot[j++]);
		}
	}
	if (j < 3) {
		printf ("Usage:\t%s <(4dfp) imgfile> <(4dfp) mskfile> <(4dfp) outfile>\n", program);
        	printf (" e.g.:\t%s -t23.2 va1234_mpr mask va1234_mpr_msk\n", program);
        	printf ("\toption\n");
        	printf ("\t-N\treplace NaN in <imgfile> with corresponding <mskfile> value\n");
        	printf ("\t-e\treport to stdout mean <imgfile> within-mask value\n");
        	printf ("\t-1\tapply first frame of <mskfile> to all frames of <imgfile>\n");
		printf ("\t-R\tsuppress creation of rec file\n");
        	printf ("\t-v<flt>\tspecify <outfile> uniform within-mask value\n");
        	printf ("\t-p<flt>\tspecify <mskfile> threshold as percent of <mskfile> max\n");
        	printf ("\t-t<flt>\tspecify <mskfile> threshold directly (default = 0.0)\n");
        	printf ("\t-A\tthreshold mask by absolute value of <mskfile>\n");
		printf ("\t-@<b|l>\toutput big or little endian (default <imgfile> endian)\n");
        	printf ("N.B.:\t<imgfile> and <mskfile> may be the same\n");
		exit (1);
	}
	if (uniform) printf ("uniform replacement value %-.4f\n", value);

/*****************************/
/* get input 4dfp dimensions */
/*****************************/
	for (j = 0; j < 3; j++) {
		sprintf (imgfile[j], "%s.4dfp.img", imgroot[j]);
		if (j < 2) {
			if (get_4dfp_dimoe (imgfile[j], imgdim[j], voxdim[j], orient + j, isbig + j) < 0)
				errr (program, imgfile[j]);
		}
		if (j == 0) if (Getifh (imgfile[j], &ifh)) errr (program, imgroot[j]);
		if (j == 1) {
			status = orient[j] != orient[0];
			i = (onemask) ? 3 : 4;
			for (k = 0; k < i; k++) status |= (imgdim[j][k] != imgdim[0][k]);
			for (k = 0; k < 3; k++) status |= (fabs (voxdim[j][k] - voxdim[0][k]) > 0.0001);
			if (status) {
				fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot[0], imgroot[1]);
				exit (-1);
			}
		}
	}

/********************************/
/* alloc buffers and open files */
/********************************/
	dimension = imgdim[0][0] * imgdim[0][1] * imgdim[0][2];
	for (j = 0; j < 2; j++) {
		if (!(imgt[j] = (float *) malloc (dimension * sizeof (float)))) errm (program);
		if (!(imgfp[j] = fopen (imgfile[j], "rb"))) errr (program, imgfile[j]);
		fprintf (stdout, "Reading: %s\n", imgfile[j]);
	}
	if (!(imgfp[2] = fopen (imgfile[j], "wb"))) errw (program, imgfile[j]);
	fprintf (stdout, "Writing: %s\n", imgfile[2]);

/***********/
/* process */
/***********/
	if (!control) control = (isbig[0]) ? 'b' : 'l';
	if (pctflag) {
		imgmax = -FLT_MAX;
		for (k = 0; k < ((onemask) ? 1 : imgdim[1][3]); k++) {
			if (eread (imgt[1], dimension, isbig[1], imgfp[1])) errr (program, imgfile[1]);
			for (i = 0; i < dimension; i++) if (imgt[1][i] > imgmax) imgmax = imgt[1][i];
		}
		threshold = pct * imgmax / 100.;
		rewind (imgfp[1]);
	}
	printf ("threshold %-.4f\n", threshold);

	if (onemask) {
		if (eread (imgt[1], dimension, isbig[1], imgfp[1])) errr (program, imgfile[1]);
	}
	for (k = 0; k < imgdim[0][3]; k++) {
		for (j = 0; j < ((onemask) ? 1 : 2); j++) {
			if (eread (imgt[j], dimension, isbig[j], imgfp[j])) errr (program, imgfile[j]);
		}
		for (q = nvox = i = 0; i < dimension; i++) {
			switch (opr) {
			case 'N' :
				if (isnan (imgt[0][i])) imgt[0][i] = imgt[1][i];
				break;
			default:
				test = (absflag) ? fabs (imgt[1][i]) > threshold : imgt[1][i] > threshold;
				if (test) {
					if (uniform) imgt[0][i] = value;
				} else {
					imgt[0][i] = 0.0;
				}
				break;
			}
			if (test) {
				q += imgt[0][i];
				nvox++;
			}
		}
		if (nvox) q /= nvox;
		if (getmean) {
			printf ("%s frame %3d within-mask voxel count: %10d\n",   imgfile[0], k + 1, nvox);
			printf ("%s frame %3d within-mask mean value:  %10.4f\n", imgfile[0], k + 1, q);
		}
		if (ewrite (imgt[0], dimension, control, imgfp[2])) errw (program, imgfile[2]);
	}
	for (j = 0; j < 3; j++) fclose (imgfp[j]);

/*******************/
/* create rec file */
/*******************/
if (saverec) {
	startrece (imgfile[2], argc, argv, rcsid, control);
	if (absflag) printrec ("absolute value ");
	sprintf (command, "%s threshold=%f\n", imgroot[1], threshold); printrec (command);
	catrec (imgfile[0]);
	if (strcmp (imgroot[0], imgroot[1])) catrec (imgfile[1]);
	endrec ();
}

/*******/
/* ifh */
/*******/
	if (Writeifh (program, imgroot[2], &ifh, control)) errw (program, imgroot[2]);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s", imgroot[2]);
	printf ("%s\n", command);
	status |= system (command);

	free (imgt[0]); free (imgt[1]);
	exit (status);
}
