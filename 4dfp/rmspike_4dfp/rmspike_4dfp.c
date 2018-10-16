/*$Header: /data/petsun4/data1/src_solaris/rmspike_4dfp/RCS/rmspike_4dfp.c,v 2.14 2007/04/11 01:54:36 avi Exp $*/
/*$Log: rmspike_4dfp.c,v $
 * Revision 2.14  2007/04/11  01:54:36  avi
 * eliminate f_init() and f_exit() in preparatiion for linux port
 *
 * Revision 2.13  2006/09/27  20:36:44  avi
 * Solaris 10
 *
 * Revision 2.12  2005/07/07  00:42:53  avi
 * remove all references to direct IFH i/o
 * correct dangerous typo
 * use standard writeifh() and also call ifh2hdr
 *
 * Revision 2.11  2004/10/11  22:05:31  rsachs
 * Installed 'errm','errr','errw','getroot','get_4dfp_dimo','setprog'. Removed 'Get4dfpDimN'.
 *
 * Revision 2.10  2001/06/29  00:10:29  avi
 * MAXF -> 2048
 *
 * Revision 2.9  2001/05/17  07:01:59  avi
 * -F (arbitrary non-functional frame run format) option
 * eliminate -V (read pre-ifh style 4dfp stacks) option
 *
 * Revision 2.8  2001/02/02  17:34:37  avi
 * correct misplaced exit (0); => $status == 255 bug
 *
 * Revision 2.7  2000/11/14  07:32:43  avi
 * extensive code update
 * eliminate -s option (save variance volume)
 * always output ifh even if reading old (-V option) input
 * use rec subroutines
 *
 * Revision 2.6  1997/05/23  00:49:24  yang
 * new rec macros
 *
 * Revision 2.5  1997/04/28  20:49:23  yang
 * working Solaris version.
 *
 * Revision 2.4  1996/08/23  00:57:11  avi
 * -r switch to call freorder for -u16 packed data
 *
 * Revision 2.3  1996/07/29  01:31:42  avi
 * x_spike and y_spike
 *
 * Revision 2.2  1996/07/10  06:05:26  avi
 * rec file saving
 *
 * Revision 2.1  1996/07/06  01:16:10  avi
 * Rev 2 initial checkin
 **/
/*________________________________________________________________________________
Module:		rmspike_4dfp.c
Date:		April 12, 1996 original program
		Revised 04-Jul-96 AZS
Authors:	Avi Snyder, Bob McKinstry, Tom Yang
Description:	Remove intensity spike from a 4dfp stack 
		and write 4dfp stack as "_rmsp" file.
__________________________________________________________________________________*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>			/* R_OK */
#include <unistd.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MAXL 	256				/* file name characters */		
#define MAXF	2048				/* frames per run */

/*************/
/* externals */
/*************/
extern int	expandf (char *string, int len);				/* expandf.c */
extern void	frmspiker_ (float *stack, int *nx, int *ny, int *nz, int *nv, int *ix, int *iy);	/* frmspike.f */

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[] = "$Id: rmspike_4dfp.c,v 2.14 2007/04/11 01:54:36 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	FILE		*fp_img, *fp_out;
	char		imgroot[MAXL], imgfile[MAXL];
	char		outroot[MAXL], outfile[MAXL];

	float		*imgs;			/* mean slice variance */
	float		*imgt;			/* frame */
	float		*imgu;			/* frame average */
	float		*imgv;			/* frame variance */
	char		format[MAXF] = "";	/* (expanded) run format */

/***********/
/* utility */
/***********/
	int		c, i, j, k, ix, iy, one = 1, orient;
	char		*ptr, command[MAXL], program[MAXL], textstr[MAXL];
	short		*sort;				/* order buffer for sorting */
	float		ftmp;

/********************/
/* image dimensions */
/********************/
	float		voxdim[3];			/* voxel dimensions in mm */
	int		imgdim[4];
	int		nf_anat = 0;
	int		nf_func;			/* number of functional frames */
	int		xdim, ydim, zdim, vdim, sdim;	/* image dimensions */
	int		isbig;
	int		x_spike = 0;			/* restrict search to this coordinate */
	int		y_spike = 0;			/* restrict search to this coordinate */
	char		control = '\0';

/*********/
/* flags */
/*********/
	int		rmspike = 0;			/* spike found flag */
	int		status = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'n':
					nf_anat = atoi (ptr);
					*ptr = '\0';
					printf ("%s: nf_anat=%d\n", program, nf_anat);
					break;
				case 'x':
					x_spike = atoi (ptr);
					*ptr = '\0';
					printf ("%s: x_spike=%d\n", program, x_spike);
					break;
				case 'y':
					y_spike = atoi (ptr);
					*ptr = '\0';
					printf ("%s: y_spike=%d\n", program, y_spike);
					break;
				case 'F': strcpy (format, ptr);	*ptr = '\0'; break;
				case '@': control = *ptr++;	*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0:	getroot (argv[i], imgroot); k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:\t%s <file_4dfp>\n", program);
		printf ("\toption\n");
		printf ("\t-n<int>\tspecify number of anatomy frames\n");
		printf ("\t-x<int>\trestrict search to specified column\n");
		printf ("\t-y<int>\trestrict search to specified row\n");
		printf ("\t-F<str>\tspecify whole run functional/non-functional format\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf (" e.g.,\t%s -n3 -x33 test_b1.4dfp.img\n", program);
		printf (" e.g.,\t%s -x33 -F\"45(1x6+)\" test_b1\n", program);
		exit (1);
	}
	sprintf (imgfile, "%s.4dfp.img", imgroot);

/***********************************/
/* get input 4dfp image dimensions */
/***********************************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';

	xdim	= imgdim[0];
	ydim	= imgdim[1];
	zdim	= imgdim[2];
	sdim	= xdim * ydim;
	vdim	= xdim * ydim * zdim;
	printf ("ANALYZE dimensions  %10d%10d%10d%10d\n",   imgdim[0], imgdim[1], imgdim[2], imgdim[3]);
	printf ("ANALYZE mmppix      %10.6f%10.6f%10.6f\n", voxdim[0], voxdim[1], voxdim[2]);

	if (!(imgt = (float *) calloc (vdim, sizeof (float)))) errm (program);
	if (!(imgu = (float *) calloc (vdim, sizeof (float)))) errm (program);
	if (!(imgv = (float *) calloc (vdim, sizeof (float)))) errm (program);
	if (!(imgs = (float *) calloc (sdim, sizeof (float)))) errm (program);
	if (!(sort = (short *) calloc (sdim, sizeof (short)))) errm (program);

/****************************/
/* expand run format format */
/****************************/
	if (!strlen (format)) {
		nf_func = imgdim[3] - nf_anat;
		sprintf (format, "%dx%d+", nf_anat, nf_func);
	}
	if (expandf (format, MAXF)) exit (-1);
	for (nf_anat = i =  0;  i < imgdim[3]; i++) if (format[i] == 'x') nf_anat++; else format[i] = '+';
	if ((nf_func = imgdim[3] - nf_anat) <= 0) {
		fprintf (stderr, "%s: no specified functional frames in %s\n", program, imgfile);
		exit (-1);
	}
	printf ("run format = %s\n", format);

/*********************************************************/
/* read 4dfp stack and compute mean and variance volumes */
/*********************************************************/
	fprintf (stdout, "Reading: %s\n", imgfile);
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);

	for (i = 0; i < imgdim[3]; i++) {
		if (eread (imgt, vdim, isbig, fp_img)) errr (program, imgfile);
		if (format[i] != 'x') for (k = 0; k < vdim; k++) {
			imgu[k] += imgt[k];
			imgv[k] += imgt[k] * imgt[k];
		}
	}
	for (k = 0; k < vdim; k++) {
		imgu[k] /= nf_func;
		imgv[k] /= nf_func;
		imgv[k] -= imgu[k] * imgu[k];
	}

/********************************/
/* average variance over slices */
/********************************/
	for (k = 0; k < sdim; k++) {
		for (j = 0; j < imgdim[2]; j++) imgs[k] += (imgv + j * sdim)[k];
		imgs[k] /= imgdim[2];
	}

/****************************************/
/* sort values in summed variance slice */
/****************************************/
	for (k = 0; k < sdim; k++) {
		if (x_spike && k % xdim != x_spike - 1) imgs[k] = 0.;	/* consider only specified column */
		if (y_spike && k / xdim != y_spike - 1) imgs[k] = 0.;	/* consider only specified row */
		sort[k] = k;
	}
	for (k = 0; k < sdim - 1; k++) {
		for (j = k + 1; j < sdim; j++) {
			if (imgs[j] > imgs[k]) {
				ftmp = imgs[k]; imgs[k] = imgs[j]; imgs[j] = ftmp;
				i    = sort[k]; sort[k] = sort[j]; sort[j] = i;
			}
		}
	}
	for (k = 0; k < 10; k++) printf ("%10.2f%10d%10d%10d\n", imgs[k], sort[k], sort[k] % xdim, sort[k] / xdim);
	ix = 1 + sort[0] % xdim;
	iy = 1 + sort[0] / xdim;
	rmspike = (imgs[0] > 2. * imgs[1] && ix != 1 && ix != xdim && iy != 1 && iy != ydim);
	if (rmspike) {
		sprintf (textstr, "Spike removed at (x,y) = (%d,%d)\n", ix, iy);
	} else {
		sprintf (textstr, "No spike found in %s\n", imgfile);
	}
	printf ("%s", textstr);

/**********************************************/
/* remove spike and write rmspiked 4dfp stack */
/**********************************************/
	sprintf (outroot, "%s_rmsp", imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);
	fprintf (stdout, "Writing: %s\n", outfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	rewind (fp_img);
	for (i = 0; i < imgdim[3]; i++) {
		if (eread  (imgt, vdim, isbig, fp_img)) errr (program, imgfile);
		if (rmspike) frmspiker_ (imgt, &xdim, &ydim, &zdim, &one, &ix, &iy);
		if (ewrite (imgt, vdim, control, fp_out)) errw (program, outfile);

	}
	fclose (fp_out);
	fclose (fp_img);

/******************/
/* create recfile */
/******************/
	startrece (outfile, argc, argv, rcsid, control);
	printrec (textstr);
	catrec (imgfile);
	endrec ();

/*********************/
/* ifh and hdr files */
/*********************/
	if (writeifhe (program, outfile, imgdim, voxdim, orient, control)) errw (program, outroot);
	sprintf (command, "ifh2hdr %s", outroot);
	status |= system (command);
	
	free (imgt);
	free (imgu);
	free (imgv);
	free (imgs);
	free (sort);
	exit (0);
}
