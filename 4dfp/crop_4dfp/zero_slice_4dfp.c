/*$Header: /data/petsun4/data1/src_solaris/crop_4dfp/RCS/zero_slice_4dfp.c,v 1.2 2010/08/20 23:21:29 avi Exp $*/
/*$Log: zero_slice_4dfp.c,v $
 * Revision 1.2  2010/08/20  23:21:29  avi
 * restore correct operation also with previous (4 argument) usage
 *
 * Revision 1.1  2010/08/20  04:43:56  avi
 * Initial revision
 *
 * Revision 1.2  2004/09/30  20:47:22  rsachs
 * Removed 'Get4dfpDimN'; replaced it with 'get_4dfp_dimo'.
 *
 * Revision 1.1  1999/06/28  22:44:24  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <unistd.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>
#include <flip_4dfp.h>

#define MAXL 256

/*************/
/* externals */
/*************/
extern void	fzero_slice_ (char *dimstr, float *imgt,
		int *nx, int *ny, int *nz, int *istart, int *iend);	/* fzero_slice.f */

int x4dfp2analyze (float *imag, int *dim, int orientation) {
	switch (orientation) {
		case 4: flipx (imag, dim+0, dim+1, dim+2);	/* sagittal */
		case 3:	flipz (imag, dim+0, dim+1, dim+2);	/* coronal */
		case 2:	flipy (imag, dim+0, dim+1, dim+2);	/* transverse */
			return 0;
		default: return -1;				/* none of the above */
	}
}

void getirange (char *command, int *minval, int *maxval) {
        char	*str;

	str = strstr (command, "to");
	if (str) {
		*str = '\0';
		*minval = atoi (command);
		*maxval = atoi (str + 2);
	} else {
		*minval = 0;
		*maxval = atoi (command);
	}
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void usage (char *program) {
	printf ("Usage:\t%s <(4dfp) image>\n", program);
	printf (" e.g.,\t%s vce20_mpr -z1to3\n", program);
	printf ("   or,\t%s vce20_mpr <x|y|z> istart iend [outroot]\n", program);
	printf ("\toption\n");
	printf ("\t-<x|y|z><int>to<int>\tspecify x y z limits (single required argument mode)\n");
	printf ("\t-f\tinterpret slice numbers using 4dfp<->analyze flips\n");
	printf ("\t-o\tspecify output fileroot (default = <image>z)\n");
	printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
	printf ("N.B.:\tslices count from 1\n");
	printf ("N.B.:\ttwo usages are supported: 1 or 4 required arguments\n");
	exit (1);
}

static char rcsid[]= "$Id: zero_slice_4dfp.c,v 1.2 2010/08/20 23:21:29 avi Exp $";
int main (int argc, char **argv) {
/*************/
/* image I/O */
/*************/
	FILE		*fp_img, *fp_out;
	IFH		ifh;
	char		imgfile[MAXL], outfile[MAXL], imgroot[MAXL], outroot[MAXL] = "";
	char		control = '\0';
	int		isbig;

/**************/
/* processing */
/**************/
	int		imgdim[4], orient, dimension;
	float		voxdim[3];
	float		*imgt;
	int		ix0 = 0, ix1 = 0, iy0 = 0, iy1 = 0, iz0 = 0, iz1 = 0;
	int		istart, iend;
	char		dirstr[4] = "";

/***********/
/* utility */
/***********/
	int		c, i, k;
	char		*str, command[MAXL], program[MAXL];

/*********/
/* flags */
/*********/
	int		doflip = 0;
	int		status = 0;
	int		zerox = 0, zeroy = 0, zeroz = 0;


	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
		strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case 'f': doflip++;			break;
				case 'x': getirange (str, &ix0, &ix1);	*str = '\0'; break;
				case 'y': getirange (str, &iy0, &iy1);	*str = '\0'; break;
				case 'z': getirange (str, &iz0, &iz1);	*str = '\0'; break;
				case 'o': getroot (str, outroot);	*str = '\0'; break;
				case '@': control = *str++;		*str = '\0'; break;
			}
		} else switch (k) {
			case 0: getroot (argv[i], imgroot);	k++; break;
			case 1: strncpy (dirstr, argv[i], 1);	k++; break;
			case 2: istart = atoi (argv[i]);	k++; break;
			case 3: iend = atoi (argv[i]);		k++; break;
			case 4: getroot (argv[i], outroot);	k++; break;
		}
	}
	if (k == 1) {
		zerox = ix0 || ix1;
		zeroy = iy0 || iy1;
		zeroz = iz0 || iz1;
	} else if (k >= 4) {
		switch (dirstr[0]) {
			case 'x': zerox++; ix0 = istart; ix1 = iend;	break;
			case 'y': zeroy++; iy0 = istart; iy1 = iend;	break;
			case 'z': zeroz++; iz0 = istart; iz1 = iend;	break;
			default: usage (program);			break;
		}
	} else {
		usage (program);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
/**************************************************/
/* parse out output filename extension if present */
/* or create output filename if not given         */
/**************************************************/
	if (!strlen (outroot)) sprintf (outroot, "%sz", imgroot); 
	sprintf (outfile, "%s.4dfp.img", outroot);

/*****************************/
/* get 4dfp input dimensions */
/*****************************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0
	||  Getifh (imgfile, &ifh)) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';

	if (zerox && (ix1 > imgdim[0] || ix1 < 1)) ix1 = imgdim[0];
	if (zeroy && (iy1 > imgdim[1] || iy1 < 1)) iy1 = imgdim[1];
	if (zeroz && (iz1 > imgdim[2] || iz1 < 1)) iz1 = imgdim[2];
	if (zerox && (ix0 > imgdim[0] || ix0 < 1)) ix0 = 1;
	if (zeroy && (iy0 > imgdim[1] || iy0 < 1)) iy0 = 1;
	if (zeroz && (iz0 > imgdim[2] || iz0 < 1)) iz0 = 1;
	if (ix1 < ix0) {k = ix1; ix1 = ix0; ix0 = k;}
	if (iy1 < iy0) {k = iy1; iy1 = iy0; iy0 = k;}
	if (iz1 < iz0) {k = iz1; iz1 = iz0; iz0 = k;}
	if (zerox) printf ("x slice limits %3d to %3d\n", ix0, ix1);
	if (zeroy) printf ("y slice limits %3d to %3d\n", iy0, iy1);
	if (zeroz) printf ("z slice limits %3d to %3d\n", iz0, iz1);

/*****************/
/* alloc buffers */
/*****************/
	dimension = imgdim[0] * imgdim[1] * imgdim[2];
	if (!(imgt = (float *) malloc (dimension * sizeof (float)))) errm (program);

/************/
/* process  */
/************/
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	fprintf (stdout, "Reading: %s\n", imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);
	for (k = 0; k < imgdim[3]; k++) {
		if (eread (imgt, dimension, isbig, fp_img)) errr (program, imgfile);
		if (doflip) {
			if (x4dfp2analyze (imgt, imgdim, ifh.orientation)) {
				fprintf (stderr, "%s: invalid %s orientation\n", program, imgroot);
				exit (-1);
			}
		}
		if (zerox) fzero_slice_ ("x", imgt, imgdim + 0, imgdim + 1, imgdim + 2, &ix0, &ix1);
		if (zeroy) fzero_slice_ ("y", imgt, imgdim + 0, imgdim + 1, imgdim + 2, &iy0, &iy1);
		if (zeroz) fzero_slice_ ("z", imgt, imgdim + 0, imgdim + 1, imgdim + 2, &iz0, &iz1);
		if (doflip) x4dfp2analyze (imgt, imgdim, ifh.orientation);
		if (ewrite (imgt, dimension, control, fp_out)) errw (program, outfile);
        }
        fclose (fp_img);
        fclose (fp_out);

/**************/
/* output ifh */
/**************/
	if (writeifhmce (program, outroot, imgdim, voxdim, ifh.orientation, ifh.mmppix, ifh.center, control)) {
		errw (program, outroot);
	}

/***************/
/* run ifh2hdr */
/***************/
	sprintf (command, "ifh2hdr %s", outroot); status |= system (command);

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, rcsid, control);
	if (doflip) printrec ("indices subjected to 4dfp<->analyze orientation dependent flips\n");
	if (ix1) {sprintf (command, "x slices %3d to %3d zeroed\n", ix0, ix1); printrec (command);}
	if (iy1) {sprintf (command, "y slices %3d to %3d zeroed\n", iy0, iy1); printrec (command);}
	if (iz1) {sprintf (command, "z slices %3d to %3d zeroed\n", iz0, iz1); printrec (command);}
	catrec (imgfile);
	endrec ();

        free (imgt);
        exit (status);
}
