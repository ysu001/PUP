/*$Header: /data/petsun4/data1/src_solaris/paste_4dfp/RCS/paste_4dfp.c,v 1.9 2008/10/17 02:53:48 avi Exp $*/
/*$Log: paste_4dfp.c,v $
 * Revision 1.9  2008/10/17  02:53:48  avi
 * read number of frames to append as third field in input list
 * allow mixed endian input stacks
 *
 * Revision 1.8  2007/08/10  22:26:07  avi
 * correct intolerance of input list blank lines introduced in last revision
 *
 * Revision 1.7  2007/07/31  06:24:20  avi
 * Solaris 10 and endian compliant
 * creates hdr
 *
 * Revision 1.6  2005/07/30  01:09:54  avi
 * eliminate limit to number of input files
 * default period 7->1
 *
 * Revision 1.5  2005/02/04  04:54:28  avi
 * allow list files to omit specification of start frame
 *
 * Revision 1.4  2004/10/08  22:17:02  rsachs
 * Installed 'errm','errr','errw','getroot','setprog','writeifh'. Removed 'Writeifh'.
 *
 * Revision 1.3  2000/02/21  05:02:27  avi
 * clean and update code
 *
 * Revision 1.2  1998/09/25  04:45:23  avi
 * append mode
 *
 * Revision 1.1  1998/09/25  03:45:57  avi
 * Revision 1.1  1997/11/26  19:40:53  tscull
 * Initial revision
 **/
/*_________________________________________________________________
  Program:	paste_4dfp

  Description:	Cut and paste 4dfp stack fragments.

  Authors:	Avi Snyder

  History:	21-Apr-97.	
_________________________________________________________________*/
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <unistd.h>	/* R_OK */
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define NFRAME		1
#define MAXL		256

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;

	if (ptr = strchr (string, '#')) 	*ptr = '\0';
	if (ptr = strchr (string, '\n'))	*ptr = '\0';
	i = m = 0;
	while (m < maxp) {
		while (!isgraph ((int) string[i]) && string[i]) i++;
		if (!string[i]) break;
		srgv[m++] = string + i;
		while (isgraph ((int) string[i])) i++;
		if (!string[i]) break;
		string[i++] = '\0';
	}
	return m;
}

typedef struct {
	char		imgfile[256];
	int		frame0;
	int		nframe;
	int		isbig;
} RUN_INFO;

static char rcsid[] = "$Id: paste_4dfp.c,v 1.9 2008/10/17 02:53:48 avi Exp $";
int main (int argc, char *argv[]) {
	FILE			*imgfp, *outfp;
	char			*srgv[256];		/* input string field pointers for split () */
	char			*ptr, command[MAXL], string[MAXL], program[MAXL];
	char			imgfile[MAXL], imgroot[MAXL], outroot[MAXL], outfile[MAXL];
	int			c, i, j, k, l, m, orient, orient0;
	char			control = '\0';
/**************/
/* stack list */
/**************/
	FILE			*lstfp;			/* input image list */
	char			lstfile[MAXL];
	RUN_INFO		*run_info;
	int			istack, nstack;

/*************************/
/* input (virtual) image */
/*************************/
	IFH			ifh;			/* one ifh for all images */
	float			*imga;
	float			mmppix[3], voxdim[3];
	int			imgdim[4];
	int			xdim, ydim, zdim, vdim, idim, dimension;
	int			nframe = NFRAME;		/* frames per epoch */
	int			nframetot = 0;			/* accumulated frames in append mode */

/****************/
/* output image */
/****************/
	float			amax = -FLT_MAX, amin = FLT_MAX;
	float			*imgs;

/*********/
/* flags */
/*********/
	int			debug = 0;
	int			status = 0;
	int			append = 0;

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
/******************************/
/* get command line arguments */
/******************************/
        for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'a': append++;		break;
				case 'p': nframe = atoi (ptr);	break;
				case '@': control = *ptr++;	*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0: strcpy (lstfile, argv[i]);	k++; break;
			case 1: getroot (argv[i], outroot);	k++; break;
		}
	}
	if (k < 2) {
		printf ("Usage:\t%s <inlist> <outfile>\n", program);
		printf ("\toption\n");
		printf ("\t-a\tappend successive epochs (default average)\n");
		printf ("\t-p<int>\tspecify period in frames (default=%d)\n", NFRAME);
		printf ("\t-@<b|l>\toutput big or little endian (default initial input endian)\n");
		exit (1);
	}

/*******************/
/* scan input list */
/*******************/
	if (!(lstfp = fopen (lstfile, "r"))) errr (program, lstfile);
	nstack = 0; while (fgets (string, (int) MAXL, lstfp)) {
		m = split (string, srgv, 256);	
		if (!m) continue;			/* blank line */
		nstack++;
	}
	run_info = (RUN_INFO *) calloc (nstack, sizeof (RUN_INFO));
	rewind (lstfp);
	for (istack = 0; istack < nstack; istack++) {
		fgets (string, (int) MAXL, lstfp);
		m = split (string, srgv, 256);	
		if (!m) {istack--; continue;}		/* blank line */
		for (k = i = 0; i < m; i++) {
			switch (k) {
				case 0:	getroot (srgv[i], imgroot);			k++; break;
				case 1:	run_info[istack].frame0 = atoi (srgv[i]) - 1;	k++; break;
				case 2:	run_info[istack].nframe = atoi (srgv[i]);	k++; break;
			}
		}
		sprintf (run_info[istack].imgfile, "%s.4dfp.img", imgroot);
	}
	fclose (lstfp);

/*************************************************/
/* check image and voxel dimensional consistency */
/*************************************************/
	for (istack = 0; istack < nstack; istack++) {
		if (get_4dfp_dimoe (run_info[istack].imgfile, imgdim, voxdim, &orient, &run_info[istack].isbig)) {
			errr (program, run_info[istack].imgfile);
		}
		if (append) {
			if (!run_info[istack].nframe) run_info[istack].nframe = nframe;
		} else {
			run_info[istack].nframe = nframe;
		}
		if (imgdim[3] < run_info[istack].frame0 + run_info[istack].nframe) {
			fprintf (stderr, "%s: insufficient frames in %s\n", program, run_info[istack].imgfile);
			exit (-1);
		}
		if (!istack) {
			xdim = imgdim[0];
			ydim = imgdim[1];
			zdim = imgdim[2];
			idim = xdim * ydim * zdim;
			orient0 = orient;
			for (k = 0; k < 3; k++) mmppix[k] = voxdim[k];
			if (Getifh (run_info[istack].imgfile, &ifh)) errr (program, run_info[istack].imgfile);
		} else {
			status = orient != orient0;
			status |= xdim != imgdim[0] || ydim != imgdim[1] || zdim != imgdim[2];
			for (k = 0; k < 3; k++) status |= (fabs (mmppix[k] - voxdim[k]) > 0.0001);
			if (status) {
				fprintf (stderr, "%s: %s %s dimensional inconsistency\n", program,
					run_info[0].imgfile, run_info[istack].imgfile);
				exit (-1);
			}
		}
	}

/************************************/
/* construct output image file name */
/************************************/
	if (!control) control = (run_info[0].isbig) ? 'b' : 'l';
	sprintf (outfile, "%s.4dfp.img", outroot);

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, rcsid, control);
	sprintf (string, "imglst\t%s\n", lstfile); printrec (string);
	catrec (lstfile);
	printrec ("endimglst\n");

/*************************/
/* allocate image memory */
/*************************/
	dimension = idim * nframe;
	imga = (float *) malloc (idim      * sizeof (float));	/* input frame buffer */
	imgs = (float *) calloc (dimension,  sizeof (float));	/* accumulator (average mode) */
	if (!imga || !imgs) errm (program);

	fprintf (stdout, "Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);

	for (i = 0; i < nstack; i++) {			/* stack loop */
		fprintf (stdout, "Reading: %s frames %3d to %3d\n",
			run_info[i].imgfile, run_info[i].frame0 + 1, run_info[i].frame0 + run_info[i].nframe);
		if (!(imgfp = fopen (run_info[i].imgfile, "r"))
		|| fseek (imgfp, (long) run_info[i].frame0 * idim * sizeof (float), SEEK_SET)) errr (program, run_info[i].imgfile);

		for (k = 0; k < run_info[i].nframe; k++) {
			if (eread (imga, idim, run_info[i].isbig, imgfp)) errr (program, run_info[i].imgfile);
			if (append) {
				for (j = 0; j < idim; j++) {
					if (imga[j] > amax) amax = imga[j];
					if (imga[j] < amin) amin = imga[j];
				}
				if (ewrite (imga, idim, control, outfp)) errw (program, outfile);
				nframetot++;
			} else {
				for (j = 0; j < idim; j++) (imgs + k*idim)[j] += imga[j];
			}
		}					/* end frame loop */
		if (fclose (imgfp)) errr (program, run_info[i].imgfile);
	}						/* end stack loop */

	if (append) {
		ifh.matrix_size[3] = nframetot;
	} else {
/*****************************************/
/* normalize and write stack accumulator */
/*****************************************/
		for (k = 0; k < dimension; k++) {
			imgs[k] /= nstack;
			if (imgs[k] > amax) amax = imgs[k];
			if (imgs[k] < amin) amin = imgs[k];
		}
		ifh.matrix_size[3] = nframe;
		if (ewrite (imgs, dimension, control, outfp)) errw (program, outfile);
	}
	if (fclose (outfp)) errw (program, outfile);

/*************/
/* write ifh */
/*************/
	printf ("output stack dimensions%10d%10d%10d%10d\n", xdim, ydim, zdim, ifh.matrix_size[3]);
	printf ("output stack mmppix    %10.6f%10.6f%10.6f\n", ifh.mmppix[0], ifh.mmppix[1], ifh.mmppix[2]);
	Writeifh (program, outfile, &ifh, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s -r%.0fto%.0f", outroot, amin, amax);
	printf ("%s\n", command);
	status |= system (command);

/*****************************************/
/* include unique contributing rec files */
/*****************************************/
	for (i = 0; i < nstack; i++) {
		for (k = j = 0; j < i; j++) k |= !strcmp (run_info[j].imgfile, run_info[i].imgfile);
		if (!k) catrec (run_info[i].imgfile);
	}
	endrec ();

	free (imgs); free (imga);
	free (run_info);
	exit (status);
}
