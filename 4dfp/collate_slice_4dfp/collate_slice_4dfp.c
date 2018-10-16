/*$Header: /data/petsun4/data1/src_solaris/collate_slice_4dfp/RCS/collate_slice_4dfp.c,v 1.8 2007/09/19 02:05:44 avi Exp $*/
/*$Log: collate_slice_4dfp.c,v $
 * Revision 1.8  2007/09/19  02:05:44  avi
 * linux, Solaris 10 and endian compliant
 *
 * Revision 1.7  2004/11/10  22:10:44  rsachs
 * Fixed problems with the *.ifh file.
 *
 * Revision 1.6  2004/11/09  21:53:31  rsachs
 * Installed 'setprog' & 'writeifh'.
 *
 * Revision 1.5  2003/02/20  18:56:53  avi
 * fix several bugs (including properly stopping at N_IN limit)
 *
 * Revision 1.4  2003/02/18  08:13:29  avi
 * generalize to multiple (up to 6) input files
 * report voxel value min and max
 *
 * Revision 1.3  2001/06/11  23:27:13  avi
 * update usage
 *
 * Revision 1.2  2001/06/11  22:04:07  avi
 * verbose mode
 *
 * Revision 1.1  2001/02/18  03:56:50  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>

#define N_IN	6		/* maximum number of input images */
#define NOUT	1
#define MAXL	256

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

float slice_mean (float *image, int n) {
	float t;
	int i;

	for (t = i = 0; i < n; i++) t += image[i];
	return (n) ? t / n : 0.;
}

static char rcsid[] = "$Id: collate_slice_4dfp.c,v 1.8 2007/09/19 02:05:44 avi Exp $";
int main (int argc, char *argv[]) {
	FILE			*imgfp[N_IN], *outfp;
	char			imgroot[N_IN + NOUT][MAXL], imgfile[MAXL], ifhfile[MAXL], outfile[MAXL];
	IFH			ifh[N_IN];
	int			isbig[N_IN];
	char			control = '\0';

/***********/
/* utility */
/***********/
	char			*str, command[MAXL], program[MAXL];
	int			c, i, j, k;

/**********************/
/* process arithmetic */
/**********************/
	float			*imgt, v, vmin, vmax;
	int			nslin[N_IN], n_in, nslout, nslout0, slicedim;

/*********/
/* flags */
/*********/
	int			debug = 0;
	int			verbose = 0;
	int			status = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) if (*argv[i] == '-') {
		strcpy (command, argv[i]); str = command;
		while (c = *str++) switch (c) {
			case 'v': verbose++; break;
			case '@': control = *str++;	*str = '\0'; break;
		}
	} else {
		if (k < N_IN + NOUT) {
			getroot (argv[i], imgroot[k]); k++;
		} else {
			fprintf (stderr, "%s: maximum allowed number of input images is %d\n", program, N_IN);
			exit (-1);
		}
	}	
	if (k < 2) {
		fprintf (stderr, "Usage:\t%s <4dfp img1> <4dfp img2> ... <4dfp imgn> <4dfp imgout>\n", program);
		fprintf (stderr, "\toption\n");
		fprintf (stderr, "\t-v\tverbose mode\n");
		fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default CPU endian)\n");
		exit (1);
	}

	n_in = k - 1;
	sprintf (outfile, "%s.4dfp.img", imgroot[n_in]);
	startrecle (outfile, argc, argv, rcsid, control);

/*************************************************/
/* read images and check dimensional consistency */
/*************************************************/
	for (nslout0 = k = 0; k < n_in; k++) {
		sprintf (ifhfile, "%s.4dfp.ifh", imgroot[k]);
		if (Getifh (ifhfile, ifh + k)) errr (program, ifhfile);
		nslout0 += ifh[k].matrix_size[2];
		isbig[k] = strcmp (ifh[k].imagedata_byte_order, "littleendian");
		if (k) {
			status = ifh[k].orientation != ifh[0].orientation;
			for (i = 0; i < 3; i++) {
				status |= (fabs (ifh[k].scaling_factor[i] - ifh[0].scaling_factor[i]) > 0.0001);
			}
			for (i = 0; i < 4; i++) if (i != 2) status |= (ifh[k].matrix_size[i] != ifh[0].matrix_size[i]);
			if (status) {
				fprintf (stderr, "%s: %s %s dimension/orientation mismatch\n",
				program, imgroot[0], imgroot[k]);
				exit (-1);
			}
		}
		sprintf (imgfile, "%s.4dfp.img", imgroot[k]);
		printf ("Reading: %s\n", imgfile);
		if (!(imgfp[k] = fopen (imgfile, "rb"))) errr (program, imgfile);
		catrec (imgfile);
	}

/***********/
/* collate */
/***********/
	slicedim = ifh[0].matrix_size[0] * ifh[0].matrix_size[1];
	if (!(imgt = (float *) malloc (slicedim * sizeof (float)))) errm (program);

	printf ("Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "wb"))) errw (program, outfile);

/******************/
/* loop on volume */
/******************/
	vmin = 1.e35; vmax = -vmin;
	for (j = 0; j <ifh[0].matrix_size[3]; j++) {
		if (verbose) printf ("processing volume %d\n", j + 1);
		for (k = 0; k < n_in; k++) nslin[k] = 0;
		nslout = 0; while (nslout < nslout0) {
			for (k = 0; k < n_in; k++) if (nslin[k] < ifh[k].matrix_size[2]) {
				if (eread (imgt, slicedim, isbig[k], imgfp[k])) errr (program, imgroot[k]);
				nslin[k]++;
				v = slice_mean (imgt, slicedim);
				for (i = 0; i < slicedim; i++) {
					if (imgt[i] < vmin) vmin = imgt[i];
					if (imgt[i] > vmax) vmax = imgt[i];
				}
				if (verbose) printf ("%-60s\tslice%3d\tmean=%10.2f\n", imgroot[k], nslin[k], v);
				if (ewrite (imgt, slicedim, control, outfp)) errw (program, outfile);
				nslout++;
			}
		}
	}
	printf ("data range %.2f to %.2f\n", vmin, vmax);
	for (k = 0; k < n_in; k++) fclose (imgfp[k]);
	fclose (outfp);

/***************************/
/* create interfile header */
/***************************/
	ifh[0].matrix_size[2] = nslout0;
	ifh[0].scaling_factor[2] /= n_in;
	sprintf (outfile, "%s.4dfp.img", imgroot[n_in]);
	writeifhe (program, outfile, ifh[0].matrix_size, ifh[0].scaling_factor, ifh[0].orientation, control);
	endrec ();

/**************************/
/* create ANALYZE 7.5 hdr */
/**************************/
	sprintf (command, "ifh2hdr %s -r%d", imgroot[n_in], (int) vmax + 0.5); system (command);

/************/
/* clean up */
/************/
	free (imgt);
        exit (status);
}
