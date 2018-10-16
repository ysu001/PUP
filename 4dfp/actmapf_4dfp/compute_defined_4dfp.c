/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/compute_defined_4dfp.c,v 1.5 2012/08/10 22:12:04 avi Exp $*/
/*$Log: compute_defined_4dfp.c,v $
 * Revision 1.5  2012/08/10  22:12:04  avi
 * typos
 *
 * Revision 1.4  2012/08/10  22:09:42  avi
 * option -z
 * replace isnormal() with isnan()
 *
 * Revision 1.3  2008/06/01  04:07:04  avi
 * linux compliant
 *
 * Revision 1.2  2006/10/01  01:24:30  avi
 * Solaris 10
 *
 * Revision 1.1  2006/08/10  03:34:34  avi
 * Initial revision
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <unistd.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>
#include <conc.h>

#define MAXL      256	

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: compute_defined_4dfp.c,v 1.5 2012/08/10 22:12:04 avi Exp $";
static char	program[MAXL];
static int	debug = 0;
static int	conc_flag = 0;

int write_4dfp (char *imgroot, char *trailer, float *imgt, int *imgdim, IFH *ifh, int argc, char *argv[], char control);	/* below */

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void usage (char* program) {
	printf ("Usage:\t%s <4dfp|conc input>\n", program);
	printf ("	option\n");
	printf ("	-z	count zero voxels as undefined (default defined)\n");
	printf ("N.B.:\t%s computes a mask of voxels that are defined over all frames\n", program);
	exit (1);
}

int main (int argc, char **argv) {
	CONC_BLOCK	conc_block;			/* conc i/o control block */
	FILE 		*imgfp;
	IFH		ifh;
	char 		imgroot[MAXL], imgfile[MAXL];
	char 		outroot[MAXL] = "", outfile[MAXL], trailer[MAXL];

        int  		imgdim[4], vdim, isbig;
	float		*imgt, *imgo;
	char		control = '\0';

/***********/
/* utility */
/***********/
	char 		command[MAXL], program[MAXL], *ptr;
	int		c, i, k;
	double		q;

/*********/
/* flags */
/*********/
	int		zero_is_undefined = 0;

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]);
			ptr = command;
			while (c = *ptr++) switch (c) {
				case 'z': zero_is_undefined++;		break;
			}
		}
		else switch (k) {
			case 0:	getroot (argv[i], imgroot);
				conc_flag = (strstr (argv[i], ".conc") == argv[i] + strlen (imgroot));
								k++; break;
		}	
	}
	if (k < 1) usage (program);

/*******************************************/
/* get 4dfp dimensions and open read/write */
/*******************************************/
	if (conc_flag) {
		conc_init (&conc_block, program);
		conc_open (&conc_block, imgroot);
		strcpy (imgfile, conc_block.lstfile);
		for (k = 0; k < 4; k++) imgdim[k] = conc_block.imgdim[k];
		if (Getifh (conc_block.imgfile0[0], &ifh)) errr (program, conc_block.imgfile0[0]);
		isbig = conc_block.isbig;
	} else {
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		if (Getifh (imgfile, &ifh)) errr (program, imgfile);
		for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];
		isbig = strcmp (ifh.imagedata_byte_order, "littleendian");
		if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
		printf ("Reading: %s\n", imgfile);
	}
	if (!control) control = (isbig) ? 'b' : 'l';
	vdim = imgdim[0]*imgdim[1]*imgdim[2];
	printf ("Reading: %s\n", imgfile);

	imgt = (float *) calloc (vdim, sizeof (float));
	imgo = (float *) calloc (vdim, sizeof (float));
	if (!imgo || !imgt) errm (program);
	for (i = 0; i < vdim; i++) imgo[i] = 1.0;	/* initialze voxels as defined */

	printf ("processing volume");
	for (k = 0; k < imgdim[3]; k++) {printf(" %d", k + 1); fflush (stdout);
		if (conc_flag) {
			conc_read_vol (&conc_block, imgt);
		} else {
			if (eread (imgt, vdim, isbig, imgfp)) errr (program, imgfile);
		}
		for (i = 0; i < vdim; i++) {
			q = (double) imgt[i];
			if (isnan (q) || imgt[i] == (float) 1.e-37)	imgo[i] = 0.0;
			if (zero_is_undefined && imgt[i] == 0.0)	imgo[i] = 0.0;
		}
	}
	printf("\n");

	imgdim[3] = 1;
	write_4dfp (imgroot, "dfnd", imgo, imgdim, &ifh, argc, argv, control);

	free (imgo); free (imgt);
	exit (0);
}

int write_4dfp (char *imgroot, char *trailer, float *imgt, int *imgdim, IFH *ifh, int argc, char *argv[], char control) {
	FILE		*imgfp;
	char		*ptr, imgfile[MAXL], outfile[MAXL], command[MAXL];
	int		i, tdim, status;

	tdim = imgdim[0]*imgdim[1]*imgdim[2]*imgdim[3];

/********************/
/* assemble outroot */
/********************/
	if (ptr = strrchr (imgroot, '/')) ptr++; else ptr = imgroot;
	sprintf (outfile, "%s_%s.4dfp.img", ptr, trailer);

/*********/
/* write */
/*********/
	printf ("Writing: %s\n", outfile);
	if (!(imgfp = fopen (outfile, "wb")) || ewrite (imgt, tdim, control, imgfp)
	|| fclose (imgfp)) errw (program, outfile);

/*******/
/* ifh */
/*******/
	writeifhmce (program, outfile, imgdim,
		ifh->scaling_factor, ifh->orientation, ifh->mmppix, ifh->center, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr -r1 %s", outfile);
	printf ("%s\n", command); status |= system (command);

/*******/
/* rec */
/*******/
	sprintf (imgfile, (conc_flag ? "%s.conc" : "%s.4dfp.img"), imgroot);
	startrecle (outfile, argc, argv, rcsid, control);
	catrec (imgfile);
	endrec ();

	sprintf (imgfile, "%s.conc", imgroot);
	return status;
}
