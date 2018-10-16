/* $Header: /data/petsun4/data1/src_solaris/interp_4dfp/RCS/interp_4dfp.c,v 1.12 2008/07/19 04:27:40 avi Exp $ */
/* $Log: interp_4dfp.c,v $
 * Revision 1.12  2008/07/19  04:27:40  avi
 * reorganize i/o so only one slice is in memory at any time
 * allow TR_slice_in to be input as 0 (interpreted as evenly spaced slices)
 * additional rec info
 *
 * Revision 1.11  2007/05/03  17:53:11  avi
 * Solaris 10; endian compliant
 *
 * Revision 1.10  2005/01/29  05:48:13  avi
 * cleaning and code standardization
 *
 * Revision 1.9  2004/11/22  21:07:59  rsachs
 * Replaced 'Get4dfpDimN' with 'get_4dfp_dimo'.
 *
 * Revision 1.8  2003/09/02  18:56:22  avi
 * -d option
 *
 * Revision 1.7  2002/11/26  01:15:09  avi
 * correct slice dependent shift computation for case of odd slices
 *
 * Revision 1.6  2002/05/06  17:56:51  jzacks
 * Deleted ghost "fclose(fp)", which was causing
 * stochastic coredumps.
 *
 * Revision 1.5  2002/04/03  19:19:49  jzacks
 * new spline calls, using 1d rather than 3d
 *
 * Revision 1.4  2002/04/03  19:02:23  jzacks
 * revised user message
 *
 * Revision 1.3  2002/04/03  18:59:22  jzacks
 * after checking slice offset, before changing user messages
 *
 * Revision 1.2  2002/02/23  20:16:01  jzacks
 *
 * Revision 1.1  2002/02/23  19:49:09  jzacks
 * Initial revision
 * */
/******************************************************************************
Program: interp_4dfp
  
Purpose:	Interpolate frames to resample to an arbitrary acquisition time.
		Simultaneously frame align.
Date:		3/22/02
By:		Jeff Zacks & Avi Snyder
******************************************************************************/
#include	<stdio.h>
#include	<math.h>
#include	<stdlib.h>
#include	<string.h>
#include	<endianio.h>
#include	<Getifh.h>
#include        <rec.h>

#define MAXL	256
#define MARGIN	16
#define SLCDIR	1

/*************/
/* externals */
/*************/
int	npad_ (int *n, int *m);			/* npad.f or imgpad.f in librms */
void	spline_  (float *data, int *n, float *d2da);				/* spline.f */
void	splintv_ (float *data, int *n, float *d2da, float *z, float *v);	/* spline.f */

/***********/
/* globals */
/***********/
static char	program[MAXL];
static char	rcsid[] = "$Id: interp_4dfp.c,v 1.12 2008/07/19 04:27:40 avi Exp $";

/*************/
/* utilities */
/*************/
void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

int main (int argc,char *argv[]) {
	char		*ptr, command[MAXL], string[MAXL];
	char		imgroot[MAXL], imgfile[MAXL], outroot[MAXL], outfile[MAXL];
	char		control = '\0';			/* output endian */
	int		status = 0;

	int		c, i, k, z, ix, iy, iz, xdim, ydim, zdim, vdim, slcdim, ivox,
			tdim_in, tdim_out, tdim_in_pad, dim_in, dim_out;
	int		slcdir = SLCDIR, imgdim[4], orient, isbig,
			margin = MARGIN, padwidth;

	float		*delta, TR_slice_in, TR_vol_in, TR_vol_out;
	float		*imgt, *imgo, *ts_in_pad, *ts_out;
	float		*d2tf, t, v, voxdim[3];

	double		q;

	FILE		*fp_img, *fp_out;

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'd': slcdir = atoi (ptr);	*ptr = '\0'; break;
				case '@': control = *ptr++;	*ptr = '\0'; break;
			}
		} else switch (k) {
		 	case 0: getroot (argv[i], imgroot);	k++; break;
		 	case 1: TR_vol_in = atof (argv[i]);	k++; break;
		 	case 2: TR_slice_in = atof (argv[i]);	k++; break;
		 	case 3: TR_vol_out = atof (argv[i]);	k++; break;
		}	
	}
	if (k < 4) {
		fprintf (stderr, "Usage:\t%s <(4dfp) image> <TR_vol_in> <TR_slice_in> <TR_vol_out>\n", program);
		fprintf (stderr, "e.g.:\t%s bold_run[.4dfp[.img]] 2.25 .136 2.5\n", program);
		fprintf (stderr, "\toptions\n");
		fprintf (stderr, "\t-d<0|1> specify slice acquisition direction (0:Inf->Sup; 1:Sup->Inf) (default=%d)\n",
					SLCDIR);
		fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default input endian)\n");
		fprintf (stderr, "N.B.: if <TR_slice_in> is input as 0 slices are spaced evenly on TR_vol\n");
		exit (1);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	sprintf (outroot, "%s_ntrp",imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);

	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgroot);
	if (!control) control = (isbig) ? 'b' : 'l';
	xdim =    imgdim[0];
	ydim =    imgdim[1];
	zdim =    imgdim[2];
	tdim_in = imgdim[3];
	tdim_out = ceil (tdim_in*TR_vol_in/TR_vol_out);
	tdim_in_pad = npad_ (&tdim_in, &margin);
	padwidth = (tdim_in_pad - tdim_in)/2;
	printf ("%s: input frame TR %.4f sec output frame TR %.4f sec\n", program, TR_vol_in, TR_vol_out);
	printf ("%s: input frame count %d  output frame count %d\n", program, tdim_in, tdim_out);
	printf ("%s: padded time series length = %d, pad width = %d\n", program, tdim_in_pad, padwidth);

	if (!(delta = (float *) malloc (zdim*sizeof (float)))) errm (program);
/******************************/
/* compute slice time offsets */
/******************************/
	if (TR_slice_in == 0) TR_slice_in = TR_vol_in/zdim;
	q = 0.5*TR_slice_in;
	for (i = 0; i < 2; i++) {
		switch (slcdir) {
			case 0:	for (z = i; z < zdim; z += 2) {
				delta[z] = q;
				q += TR_slice_in;
				} break;
			case 1:	for (z = zdim - 1 - i; z >= 0; z -= 2) {
				delta[z] = q;
				q += TR_slice_in;
				} break;
		}
	}
	for (z = 0; z < zdim; z++) {
		printf ("%s: time_shift for plane %d  = %f\n", program, z, delta[z]);
	}

/******************************/
/* allocate memory for images */
/******************************/
	slcdim	= xdim*ydim;
	vdim	= xdim*ydim*zdim;
	dim_in 	= slcdim*tdim_in;
	dim_out = slcdim*tdim_out;
	if (!(imgt = (float *) calloc (dim_in,  sizeof (float)))) errm (program);
	if (!(imgo = (float *) calloc (dim_out, sizeof (float)))) errm (program);

/************************/
/* allocate time series */
/************************/
	if (!(ts_out =    (float *) calloc (tdim_out,    sizeof (float)))) errm (program);
	if (!(ts_in_pad = (float *) calloc (tdim_in_pad, sizeof (float)))) errm (program);
	if (!(d2tf =      (float *) calloc (tdim_in_pad, sizeof (float)))) errm (program);

/****************/
/* initizee i/o */
/****************/
	printf ("Reading: %s\n", imgfile);
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	printf ("Writing: %s\n", outfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	for (k = 0; k < zdim; k++) if (ewrite (imgo, dim_out, control, fp_out)) errw (program, outfile);
	rewind (fp_out);

/***************/
/* interpolate */
/***************/
	printf ("Interpolating slice:");
	for (iz = 0; iz < zdim; iz++) {
		printf (" %d", iz + 1); fflush (stdout);
		for (k = 0; k < imgdim[3]; k++) {
			if (fseek (fp_img, (long) sizeof(float)*(k*vdim + iz*slcdim), SEEK_SET)
			||  eread (imgt + k*slcdim, slcdim, isbig, fp_img)) errr (program, imgfile);
		}
		for (iy = 0; iy < ydim; iy++) {
		for (ix = 0; ix < xdim; ix++) {
			ivox = ix + xdim*iy;
			for (i = 0; i < padwidth; i++) ts_in_pad[i] = (imgt + ivox)[0];
			for (i = 0; i < tdim_in; i++) {
				ts_in_pad[i + padwidth] = (imgt + ivox)[i*slcdim];
			}
			for (i = 0; i < tdim_in_pad - tdim_in - padwidth; i++) {
				ts_in_pad[tdim_in_pad - i - 1] = (imgt + ivox)[slcdim*(tdim_in - 1)];
			}
			spline_ (ts_in_pad, &tdim_in_pad, d2tf);
			for (i = 0; i < tdim_out; i++) {
				t = i*((float) tdim_in/(float) tdim_out) + (float) padwidth - delta[iz]/TR_vol_in;
				splintv_ (ts_in_pad, &tdim_in_pad, d2tf, &t, &v);
				ts_out[i] = v;
			}
			for (i = 0; i < tdim_out; i++) (imgo + ivox)[i*slcdim] = ts_out[i];
		}}
		for (k = 0; k < tdim_out; k++) {
			if (fseek (fp_out, (long) sizeof(float)*(k*vdim + iz*slcdim), SEEK_SET)
			||  ewrite (imgo + k*slcdim, slcdim, control, fp_out)) errw (program, outfile);
		}
	}
	printf ("\n");
	if (fclose (fp_img)) errr (program, imgfile);
	if (fclose (fp_out)) errw (program, outfile);

/***************/
/* ifh and hdr */
/***************/
	imgdim[3] = tdim_out;
	if (writeifhe (program, outfile, imgdim, voxdim, orient, control)) errw (program, outroot);
	sprintf (command, "ifh2hdr -r2000 %s", outroot); status = system (command);

/************/
/* rec file */
/************/
	startrece (outfile, argc, argv, rcsid, control);
	sprintf (string, "TR_vol %.4f -> %.4f sec\n", TR_vol_in, TR_vol_out); printrec (string);
	sprintf (string, "frame count %d -> %d\n", tdim_in, tdim_out); printrec (string);
	for (z = 0; z < zdim; z++) {
		sprintf (string, "time_shift for plane %d  = %f\n", z, delta[z]);
		printrec (string);
	}
	catrec (imgfile);
	endrec ();

	free (delta); free (imgt); free (imgo);
	free (ts_in_pad); free (ts_out); free (d2tf);
	exit (0);
}
