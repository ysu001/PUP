/*$Header: /data/petsun4/data1/src_solaris/cs2ap_4dfp/RCS/cs2ap_4dfp.c,v 1.4 2010/08/29 03:54:27 avi Exp $*/
/*$Log: cs2ap_4dfp.c,v $
 * Revision 1.4  2010/08/29  03:54:27  avi
 * correct filname argument errors
 *
 * Revision 1.3  2010/08/29  02:57:37  avi
 * extensive rewrite
 * linux/endian compliant
 * separate output of amplitude and phase volumes
 *
 * Revision 1.2  2001/12/21  02:49:08  avi
 * -t (threshold)
 * FWHM defaults to 0
 * output volume lables in rec file
 *
 * Revision 1.1  1999/03/10  02:23:30  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define IMGN		2		/* number of input images */
#define IMGO		2		/* number of output images */
#define MAXL		256
#define FWHM		0.0
#define THRESH		0.0
#define ALPHA		0.1		/* relative height of smoothing kernel cutoff */

/*************/
/* externals */
/*************/
extern void	imgblur3d_ (float *fwhm, float *alpha, float *voxsiz, float *img0, int *imgdim, float *img1);	/* fimgblur.f */

/***********/
/* globals */
/***********/
static char rcsid[] = "$Id: cs2ap_4dfp.c,v 1.4 2010/08/29 03:54:27 avi Exp $";
int main (int argc, char *argv[]) {
	FILE		*imgfp;
	IFH		ifh;
	char		*str, command[MAXL], program[MAXL];
	char		imgroot[IMGN + 1][MAXL], imgfile[IMGN][MAXL], outroot[MAXL], outfile[MAXL];
	int		imgdim[IMGN][4], vdim, orient, isbig[IMGN];
	char		control = '\0';	
	float		voxdim[IMGN][3];		/* voxel dimensions in mm */
	float		*imgt, *imgo;
	float		fwhm = FWHM, thresh = THRESH, alpha = ALPHA;

/***********/
/* utility */
/***********/
	int		c, i, j, k;
	float		degperad;

/*********/
/* flags */
/*********/
	int		debug = 0;
	int		status = 0;

	printf ("%s\n", rcsid);
	if (!(str = strrchr (argv[0], '/'))) str = argv[0]; else str++;
	strcpy (program, str);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case 't': thresh = atof (str);	*str = '\0'; break;
				case 'w': fwhm   = atof (str);	*str = '\0'; break;
			}
		} else {
			if (k < IMGN + 1) getroot (argv[i], imgroot[k++]);
		}
	}
	if (k < 3) {
		printf ("Usage:\t%s <(4dfp) cos_img> <(4dfp) sin_img> <(4dfp) outroot>\n", program);
		printf ("	[option]\n");
		printf ("	-t<flt>	specify amplitude threshold for phase map (default = %6.4f)\n", THRESH);
		printf ("	-w<flt>	specify pre-blur FWHM in mm (default = %6.4f)\n", FWHM);
		printf ("	-@<b|l>\toutput big or little endian (default input endian)\n");
		exit (1);
	}

/*****************************************************/
/* read images and check for dimensional consistency */
/*****************************************************/
	for (j = 0; j < IMGN; j++) {
		if (get_4dfp_dimoe (imgroot[j], imgdim[j], voxdim[j], &orient, isbig + j)) errr (program, imgfile[j]);
		if (!j) {
			if (!control) control = (isbig) ? 'b' : 'l';
			if (Getifh (imgroot[j], &ifh)) errr (program, imgroot[j]);
			orient = ifh.orientation;
			vdim = imgdim[0][0] * imgdim[0][1] * imgdim[0][2];
			if (!(imgt = (float *) malloc (IMGN*vdim*sizeof (float)))
			||  !(imgo = (float *) malloc (IMGO*vdim*sizeof (float)))) errm (program);
		} else {
			status = orient - ifh.orientation;
			for (k = 0; k < 3; k++) status |= (imgdim[j][k] != imgdim[0][k]);
			if (status) {
				fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgfile[0], imgfile[j]);
				exit (-1);
			}
		}
		sprintf (imgfile[j], "%s.4dfp.img", imgroot[j]);
		fprintf (stdout, "Reading: %s\n", imgfile[j]);
		if (!(imgfp = fopen (imgfile[j], "rb")) || eread (imgo, vdim, isbig[j], imgfp)
		|| fclose (imgfp)) errr (program, imgfile[j]);
		if (fwhm == 0.0) {
			for (k = 0; k < vdim; k++) (imgt + j*vdim)[k] = imgo[k];
		} else {
			imgblur3d_ (&fwhm, &alpha, voxdim[j], imgo, imgdim[j], imgt + j*vdim);
		}
	}

/***********/
/* compute */
/***********/
	for (k = 0; k < vdim; k++) {
		imgo[k] = sqrt (imgt[k]*imgt[k] + (imgt + vdim)[k]*(imgt + vdim)[k]);
		(imgo + vdim)[k] = (imgo[k] < thresh) ? 0.0 : atan2 ((imgt + vdim)[k], imgt[k]);
	}

/*********/
/* write */
/*********/
	for (j = 0; j < IMGO; j++) {
		switch (j) {
			case 0: sprintf (outroot, "%s_amp", imgroot[2]); break;
			case 1: sprintf (outroot, "%s_pha", imgroot[2]); break;
		}
		sprintf (outfile, "%s.4dfp.img", outroot);
		fprintf (stdout, "Writing: %s\n", outfile);
		if (!(imgfp = fopen (outfile, "wb")) || ewrite (imgo + j*vdim, vdim, control, imgfp)
		|| fclose (imgfp)) errw (program, outfile);
	
/***************/
/* ifh hdr rec */
/***************/
		if (Writeifh (program, outfile, &ifh, control)) errw (program, outroot);
		sprintf (command, "ifh2hdr %s", outroot);
		status |= system (command);

		startrece (outfile, argc, argv, rcsid, control);
		catrec    (imgfile[0]);
		catrec    (imgfile[1]);
		endrec    ();
	}

	free (imgt); free (imgo);
	exit (0);
}
