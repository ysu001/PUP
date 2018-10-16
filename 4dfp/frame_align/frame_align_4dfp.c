/*$Header: /data/petsun4/data1/src_solaris/frame_align/RCS/frame_align_4dfp.c,v 1.23 2011/07/26 04:47:57 avi Exp $ 
  $Log: frame_align_4dfp.c,v $
 * Revision 1.23  2011/07/26  04:47:57  avi
 * use new subroutine slice_sequence ()
 * option -N enables Siemens slice sequence ordering
 * default is non-reverse slice direction
 *
 * Revision 1.22  2010/04/01  20:08:17  avi
 * silly typo
 *
 * Revision 1.21  2010/04/01  20:05:58  avi
 * address data on disk using long arithmetic
 *
 * Revision 1.20  2009/08/28  03:43:18  avi
 * option -S (correct sequentially acquired slices)
 *
 * Revision 1.19  2007/09/19  01:06:51  avi
 * correctly use endian invariant 4dfp i/o
 *
 * Revision 1.18  2007/09/19  00:57:02  avi
 * Solaris 10 and linux compliant
 *
 * Revision 1.17  2005/01/28  23:35:10  avi
 * correct failure to create output ifh with proper outroot code
 *
 * Revision 1.16  2005/01/28  06:03:52  avi
 * extensive cleaning
 * eliminate libmri and non-standard i/o
 *
 * Revision 1.15  2003/06/28  03:04:12  avi
 * correct bug to enable -d option
 *
 * Revision 1.14  2003/05/08  22:13:42  avi
 * -d (compensate for ascending or descending slice sequence)
 *
 * Revision 1.13  2003/05/08  20:44:44  avi
 * code cleaning and modernization
 *
 * Revision 1.12  2002/11/20  07:08:24  avi
 * typos
 *
 * Revision 1.11  2002/11/20  07:03:56  avi
 * use rec.o subroutines
 *
 * Revision 1.10  2002/11/20  05:48:09  avi
 * correct computation of slice timing offsets (zdelta[])
 * to fix error when slice count is odd
 *
 * Revision 1.9  2001/05/24  04:31:29  avi
 * rec file first line second field now output filename
 *
 * Revision 1.8  2000/09/18  19:37:26  jmo
 * Fix bug in setting N.
 *
 * Revision 1.7  2000/09/14  20:04:39  jmo
 * Fix bug that occurred when number of frames were greater than 256.
 *
 * Revision 1.5  1999/07/22  19:59:24  jmo
 * Revision 1.4  1999/07/22  04:55:38  jmo
 * Revision 1.3  1999/07/19  21:05:13  jmo
 * Add command line arguments to specify TRs
 *
 * Revision 1.2  1998/01/22  23:54:15  jmo
 * Revision 1.1  1998/01/22  23:01:45  jmo
 * Initial revision
 * */
/*************************************************************
Program: frame_align_4dfp
Purpose: Interpolate frames to a common reference point in time.
Date: December 29, 1997
By: John Ollinger
***************************************************************/
#include	<stdlib.h>
#include	<stdio.h>
#include	<math.h>
#include	<stdlib.h>
#include	<string.h>
#include	<Getifh.h>
#include	<endianio.h>
#include	<rec.h>

#define MAXL		256
#define FORWARD		1
#define REVERSE		-1
#define TR_VOL		2.36

/*************/
/* externals */
/*************/
void realft (double *seq,int tdim,int direction);	/* rfft1d.c */
void slice_sequence (int interleaved, int force0, int reverse, int *seq, int nslice);	/* slice_sequnce.c */

/***********/
/* globals */
/***********/
static char rcsid[] = "$Id: frame_align_4dfp.c,v 1.23 2011/07/26 04:47:57 avi Exp $";

int main (int argc, char *argv[]) {
	FILE		*fpr, *fpw;
	char		imgroot[MAXL], imgfile[MAXL], outroot[MAXL], outfile[MAXL];
	char		*ptr, command[MAXL], program[MAXL];

	int		i, k, z, lenslc, lenvol, dir,
			skip, lenseq, t, zdim, tdim, N, tn,
			seqoff, pix, *sequence;

	long		zoff, toff, dskptr;

	float		*delta, *seqimg, theta, real, imag, a, b,
			TR_slice = 0.0,
			TR_vol = TR_VOL;

	double		*seq, *seqf, q;

/********************/
/* image properties */
/********************/
	int		imgdim[4], orient, isbig;
	float		voxdim[3];
	char		control;

/*********/
/* flags */
/*********/
	int		interleaved = 1;
	int		force0 = 1;
	int		reverse = 0;
	int		status = 0;
	int		debug = 0;

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);
	if (debug) printf ("M_PI=%f\n", M_PI);

	if (argc < 3) {
		fprintf (stderr,"Usage: %s <(4dfp) input> <frames_to_skip> [options]\n", program);
		fprintf (stderr,"\t[option]\n");
		fprintf (stderr,"\t-N		enable ");
		fprintf (stderr, "interleaved order 2,4,6,...,1,3,5,... for even total slice counts\n");
		fprintf (stderr,"\t-S\t\tspecify sequential slice acquisition (default interleaved)\n");
		fprintf (stderr,"\t-d <0|1>\tspecify slice acquisition direction (0:Inf->Sup; 1:Sup->Inf) (default=%d)\n", 0);
		fprintf (stderr,"\t-TR_vol <flt>	specify frame TR in sec (default=%.2f)\n", TR_VOL);
		fprintf (stderr,"\t-TR_src <flt>	specify slice TR in sec (default=TR_vol/nslice)\n");
		fprintf (stderr,"e.g.:	%s bold_run.4dfp.img 4\n", program);
		fprintf (stderr,"	%s bold_run.4dfp.img 4 -TR_vol 2.5\n", program);
		fprintf (stderr,"	%s bold_run.4dfp.img 4 -TR_vol 2.5 -TR_slc .136\n", program);
		fprintf (stderr,"N.B.:	space between option and value\n");
		exit (1);
	}
	getroot (argv[1], imgroot);
	sprintf (outroot, "%s_faln", imgroot);
	skip = atoi (argv[2]);
	for (i = 3; i < argc; i++) {
		if (!strcmp (argv[i], "-S"))		{interleaved = 0;}
		if (!strcmp (argv[i], "-d"))		{i++; reverse	= atoi(argv[i]);}
		if (!strcmp (argv[i], "-N"))		{force0 = 0;}
		if (!strcmp (argv[i], "-TR_vol"))	{i++; TR_vol	= atof(argv[i]);}
		if (!strcmp (argv[i], "-TR_slc"))	{i++; TR_slice	= atof(argv[i]);}
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);
/**************************************************************/
/* read input file dimensions and optionally compute TR_slice */
/**************************************************************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	control = (isbig) ? 'b' : 'l';
	lenslc = imgdim[0]*imgdim[1];
	lenvol = lenslc*imgdim[2];
	lenseq = lenslc*imgdim[3];
	zdim = imgdim[2];
	tdim = imgdim[3];
	if (TR_slice == 0.0) TR_slice = TR_vol/zdim;

	if (!(fpr = fopen (imgfile,"r"))) errr (program, imgfile);
	printf ("Reading: %s\n", imgfile);
	if (!(fpw = fopen (outfile,"w"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);

	if (!(seqimg = (float *) malloc (lenseq * sizeof (float)))) errm (program);
	N = 32; while (N < tdim) N *= 2;
	if (debug) printf ("N=%d\n", N);
	if (!(seqf = (double *) malloc ((2*N+3) * sizeof (double)))) errm (program);
	seq = seqf;
	++seq;
	if (!(delta =  (float *) calloc (zdim, sizeof (float)))) errm (program);
	if (!(sequence = (int *) calloc (zdim, sizeof (int))))   errm (program);

	startrece (outfile, argc, argv, rcsid, control);
	sprintf  (command, "Slice TR: %f sec; Frame TR: %f sec\n", TR_slice, TR_vol); printrec (command);
	fprintf  (stdout,  "Slice TR: %f sec; Frame TR: %f sec\n", TR_slice, TR_vol);

/******************************/
/* compute slice time offsets */
/******************************/
	slice_sequence (interleaved, force0, reverse, sequence, zdim);
	for (z = 0; z < zdim; z++) delta[sequence[z]] = (z + 0.5)*TR_slice;
	for (z = 0; z < zdim; z++) {
		sprintf (command, "time_shift for plane %2d  = %f\n", z, delta[z]);
		printrec (command);
	}

/******************/
/* loop on planes */
/******************/
	for (z = 0, zoff = 0; z < zdim; z++, zoff += lenslc) {
		if (debug) printf ("processing plane %d\n", z + 1);

/****************************************/
/* assemble time-sequence for one plane */
/****************************************/
		for (t = 0, toff = 0, seqoff = 0; t < tdim; t++, toff += lenvol, seqoff += lenslc) {
			dskptr = sizeof (float) * (zoff + toff);
			if (fseek (fpr, dskptr, SEEK_SET)
			||  eread (seqimg + seqoff, lenslc, isbig, fpr)) errr (program, imgfile);
		}

/************************************/
/* shift time-course for each voxel */
/************************************/
		for (pix = 0; pix < lenslc; pix++) {
			for (t = 0, toff = pix; t < tdim; t++, toff += lenslc)	seq[t] = seqimg[toff];
/******************************************************/
/* skip first few frames (magnetization steady-state) */
/******************************************************/
			for (t = 0;      t < skip;   t++)			seq[t] = seq[skip];
			for (t = tdim;   t < 2*tdim; t++)			seq[t] = seq[2*tdim-t-1];
			for (t = 2*tdim; t < 2*N;    t++)			seq[t] = seq[0];
			dir = FORWARD;
			realft (seqf, N, dir);
			for (t = 1, tn = 2; t < N; t++, tn += 2) {
				theta = M_PI*delta[z]*t/((tdim-1)*TR_vol);
				real = (float) cos (theta);
				imag = (float) sin (theta);
				a = seq[tn];
				b = seq[tn+1];
				seq[tn]   = a*real - b*imag;
				seq[tn+1] = a*imag + b*real;
			}
			dir = REVERSE;
			realft (seqf, N, dir);
			for (t = skip, toff = pix + skip*lenslc; t < tdim; t++, toff += lenslc) {
				seqimg[toff] = (float) seq[t]/(float) N;
			}
		}

/********************************************/
/* write shifted sequences back out to disk */
/********************************************/
		for (t = toff = seqoff = 0; t < tdim; t++, toff += lenvol, seqoff += lenslc) {
			dskptr = sizeof (float) * (zoff + toff);
			if (fseek (fpw, (long) dskptr, SEEK_SET)
			|| ewrite (seqimg + seqoff, lenslc, control, fpw)) errw (program, outfile);
		}
	}	/* z */
	fclose (fpw);
	fclose (fpr);

/************/
/* copy ifh */
/************/
	sprintf (command, "/bin/cp %s.4dfp.ifh %s.4dfp.ifh", imgroot, outroot);
	status |= system (command);

/*******/
/* rec */
/*******/
	catrec (imgfile);
	endrec ();

	free (seqimg); free (seqf); free (sequence); free (delta);
	exit (0);
}
