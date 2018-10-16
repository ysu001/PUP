/*$Id: asciito4dfp.c,v 1.1 2007/09/08 21:00:58 avi Exp $*/
/*$Log: asciito4dfp.c,v $
 * Revision 1.1  2007/09/08  21:00:58  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <rec.h>
#include <Getifh.h>
#include <endianio.h>

#define MAXL 	256
#define MAXS	4096			/* maximum length of input line */

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;

	if (ptr = strchr (string, '#')) *ptr = '\0';
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

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[] = "$Id: asciito4dfp.c,v 1.1 2007/09/08 21:00:58 avi Exp $";
int main (int argc, char *argv[]) {
	FILE		*txtfp, *imgfp;
	IFH		ifh;
	char		txtfile[MAXL];			/* input ascii text */
	char		outroot[MAXL], outfile[MAXL];	/* output 4dfp */
	char		control = '\0';			/* output endian */

/***************/
/* 4dfp arrays */
/***************/
	float		*imgt;				/* frame buffer */
	float		voxdim[3]= {1., 1., 1.};	/* placeholder voxel dimensions in mm */
	int		imgdim[4] = {0, 1, 1, 0}, vdim, orient = 2;
	int		npts, ncol;

/***********/
/* utility */
/***********/
	int		c, i, j, k, m;
	char		*ptr, command[MAXL], string[MAXL], program[MAXL], *srgv[MAXL];

/*********/
/* flags */
/*********/
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
			}
		}
		else switch (k) {
			case 0:	strcpy (txtfile, argv[i]);	k++; break;
			case 1:	getroot (argv[i], outroot);	k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:	%s <text file> <(4dfp) out>\n", program);
		printf ("\toption\n");
		printf ("\t-@<b|l> output big or little endian (default CPU endian)\n");
		printf ("N.B.:	columns in <text file> map to voxels in <(4dfp) out>\n");
		printf ("N.B.:	'#' in <text file> introduce comments\n");
		printf ("N.B.:	<text file> lines beginning with '#' are included in <(4dfp) out>.img.rec\n");
		exit (1);
	}

/****************/
/* read txtfile */
/****************/
	printf ("Reading: %s\n", txtfile);
	if (!(txtfp = fopen (txtfile, "r"))) errr (program, txtfile);
	npts = 0; while (fgets (string, MAXS, txtfp)) {
		if (!(m = split (string, srgv, MAXL))) continue;
		if (!npts) {
			ncol = m;
		} else {
			if (m != ncol) {
				fprintf (stderr, "%s: %s format error\n", program, txtfile);
				exit (-1);
			}
		}
		npts++;
	}
	printf ("npts=%d ncol=%d\n", npts, ncol);
	rewind (txtfp);

	if (!(imgt = (float *) malloc (ncol * sizeof (float)))) errm (program);
	sprintf (outfile, "%s.4dfp.img", outroot);
	if (!(imgfp = fopen (outfile, "wb"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);

/**********************/
/* process all frames */
/**********************/
	while (fgets (string, MAXS, txtfp)) {
		if (!(m = split (string, srgv, MAXL))) continue;
		for (i = 0; i < ncol; i++) imgt[i] = atof (srgv[i]);
		if (gwrite ((char *) imgt, sizeof (float), ncol, imgfp, control)) errw (program, outfile);
	}
	if (fclose (imgfp)) errw (program, outfile);
	free (imgt);

/*******/
/* ifh */
/*******/
	imgdim[0] = ncol;
	imgdim[3] = npts;
	writeifhe (program, outroot, imgdim, voxdim, orient, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s", outroot);
	printf ("%s\n", command);
	status |= system (command);

/*******/
/* rec */
/*******/
	startrece (outfile, argc, argv, rcsid, control);
	rewind (txtfp);
	while (fgets (string, MAXS, txtfp)) {
		if (string[0] == '#') printrec (string);
	}
	endrec ();

	if (fclose (txtfp)) errr (program, txtfile);
	exit (status);
}

