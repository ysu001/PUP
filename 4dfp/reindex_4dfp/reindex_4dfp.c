/*$Header: /data/petsun4/data1/src_solaris/reindex_4dfp/RCS/reindex_4dfp.c,v 1.7 2010/03/31 01:33:21 avi Exp $*/
/*$Log: reindex_4dfp.c,v $
 * Revision 1.7  2010/03/31  01:33:21  avi
 * remove mmppix and center from output ifh
 *
 * Revision 1.6  2007/08/16  04:15:15  avi
 * index rotration
 * ifh includes mmppix and center
 *
 * Revision 1.5  2007/08/15  23:54:53  avi
 * Solaris10, endian and linux compliant
 *
 * Revision 1.4  2004/09/23  20:28:54  rsachs
 * Installed 'errm','errr','errw','setprog' & 'get_4dfp_dimo'.
 *
 * Revision 1.3  2004/05/02  02:03:32  avi
 * -o option
 *
 * Revision 1.2  2000/09/08  07:01:24  avi
 * correct ifh filename
 *
 * Revision 1.1  2000/09/08  05:21:54  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MAXL	256

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static void usage (char *program) {
	printf ("Usage:\t%s <(4dfp> input> <index1> <index2> [options]\n", program);
	printf ("\toption\n");
	printf ("\t-v\tverbose mode\n");
	printf ("\t-o<str>\tspecify 4dfp output root (default = <input>_r<index1><index2>)\n");
	printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
	printf ("e.g.,\t%s my4Dstack 3 4\n", program);
	printf ("N.B.:\t%s swaps specified indices\n", program);
	printf ("N.B.:\t<index1> and <index2> must be unequal integers in the range 1-4 except as follows\n");
	printf ("\t<index1> == 4 and <index2> == 0: right rotate indices (first index <-  last index)\n");
	printf ("\t<index1> == 0 and <index2> == 4:  left rotate indices ( last index <- first index)\n");
	exit (1);
}

static char rcsid[] = "$Id: reindex_4dfp.c,v 1.7 2010/03/31 01:33:21 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	FILE 		*fp_img, *fp_out;
	IFH		ifh;
	char            imgroot[MAXL], imgfile[MAXL];
	char            outroot[MAXL] = "", outfile[MAXL];

/**************/
/* processing */
/**************/
	int		imgdim[4], outdim[4], dimension, orient, isbig;
	float		imgvox[3], outvox[3];
	float		*fptr, *imgt, *imgo;
	char		control = '\0';

/***********/
/* utility */
/***********/
	int		c, i, j, k, l, m;
	int		jndex, kindex;
	char		*ptr, command[MAXL], program[MAXL];

/*********/
/* flags */
/*********/
	int		status = 0;
	int		verbose = 0;
	int		i1, i2, ii;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'v': verbose++;			break;
				case 'o': getroot (ptr, outroot);	*ptr = '\0'; break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0:	getroot (argv[i], imgroot);	k++; break;
			case 1:	i1 = atoi (argv[i]);		k++; break;
			case 2:	i2 = atoi (argv[i]);		k++; break;
		}	
	}
	if (k < 3) usage (program);
	if (i1 == i2 || i1 < 0 || i1 > 4 || i2 < 0 || i2 > 4) usage (program);

	if (!strlen (outroot)) sprintf (outroot, "%s_r%d%d", imgroot, i1, i2);
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);

/*********************/
/* initiate file I/O */
/*********************/
	if (get_4dfp_dimoe (imgfile, imgdim, imgvox, &orient, &isbig) < 0) errr (program, imgfile);
	Getifh (imgfile, &ifh);
	if (!control) control = (isbig) ? 'b' : 'l';
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	fprintf (stdout, "Reading: %s\n", imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);

/***********/
/* process */
/***********/
	ii = 10*i1 + i2;
	if (ii == 34 || ii == 43) {
		dimension = imgdim[0] * imgdim[1];
		if (!(imgt = (float *) malloc (dimension * sizeof (float)))) errm (program);
		for (k = 0; k < imgdim[2]; k++) for (l = 0; l < imgdim[3]; l++) {
			m = k + imgdim[2] * l;
			fseek (fp_img, (long) m*dimension * sizeof (float), SEEK_SET);
			if (eread  (imgt, dimension, isbig,   fp_img))	errr (program, imgfile);
			if (ewrite (imgt, dimension, control, fp_out))	errw (program, outfile);
		}
		free (imgt);
	} else if (ii == 12 || ii == 21) {
		dimension = imgdim[0] * imgdim[1];
		if (!(imgt = (float *) malloc (dimension * sizeof (float)))) errm (program);
		if (!(imgo = (float *) malloc (dimension * sizeof (float)))) errm (program);
		for (l = 0; l < imgdim[3]; l++) for (k = 0; k < imgdim[2]; k++) {
			if (eread  (imgt, dimension, isbig, fp_img)) errr (program, imgfile);
			fptr = imgo;
			for (i = 0; i < imgdim[0]; i++) for (j = 0; j < imgdim[1]; j++) {
				m = i + imgdim[0] * j;
				*fptr++ = imgt[m];
			}
			if (ewrite (imgo, dimension, control, fp_out)) errw (program, outfile);
		}
		free (imgt);
		free (imgo);
	} else if (ii == 40) {
		if (!(imgt = (float *) malloc (imgdim[3] * sizeof (float)))) errm (program);
		printf ("Reindexing slice");
		dimension = imgdim[0] * imgdim[1] * imgdim[2];
		for (k = 0; k < dimension; k++) {
			if (!(k % (imgdim[0]*imgdim[1]))) printf (" %d", k/(imgdim[0]*imgdim[1]) + 1); fflush (stdout);
			for (l = 0; l < imgdim[3]; l++) {
				m = l*dimension + k;
				fseek (fp_img, (long) m*sizeof (float), SEEK_SET);
				if (eread  (imgt + l, 1, isbig, fp_img)) errr (program, imgfile);
			}
			if (ewrite (imgt, imgdim[3], control, fp_out)) errw (program, outfile);
		}
		printf ("\n");
		free (imgt);
	} else if (ii == 04) {
		if (!(imgt = (float *) calloc (imgdim[0], sizeof (float)))) errm (program);
		dimension = imgdim[1] * imgdim[2] * imgdim[3];
		for (k = 0; k < dimension; k++) if (ewrite  (imgt, imgdim[0], control, fp_out)) errw (program, outfile);
		printf ("Reindexing slice");
		for (k = 0; k < dimension; k++) {
			if (!(k % (imgdim[1]*imgdim[2]))) printf (" %d", k/(imgdim[1]*imgdim[2]) + 1); fflush (stdout);
			if (eread (imgt, imgdim[0], isbig, fp_img)) errr (program, imgfile);
			for (l = 0; l < imgdim[0]; l++) {
				m = l*dimension + k;
				fseek (fp_out, (long) m*sizeof (float), SEEK_SET);
				if (ewrite  (imgt + l, 1, control, fp_out)) errw (program, outfile);
			}
		}
		printf ("\n");
		free (imgt);
	} else {
		fprintf (stderr, "%s: index %d %d swapping not yet supported\n", program, i1, i2);
		unlink (outfile);
		exit (-1);
	}
	fclose (fp_img);
	fclose (fp_out);

/*******************/
/* create rec file */
/*******************/
	for (k = 0; k < 3; k++) outvox[k] = imgvox[k];
	startrece (outfile, argc, argv, rcsid, control);
	if (ii == 40) {
		sprintf (command, "indices right rotated\n"); printrec (command);
		for (k = 0; k < 4; k++) outdim[k] = imgdim[(k + 3) % 4];
	} else if (ii == 04) {
		sprintf (command, "indices left rotated\n"); printrec (command);
		for (k = 0; k < 4; k++) outdim[k] = imgdim[(k + 1) % 4];
	} else {
		sprintf (command, "indices %d and %d swapped\n", i1, i2); printrec (command);
		for (k = 0; k < 4; k++) outdim[k] = imgdim[k];
		outdim[i1 - 1] = imgdim[i2 - 1];
		outdim[i2 - 1] = imgdim[i1 - 1];
		if (i1 < 4 && i2 < 4) {
			outvox[i1 - 1] = imgvox[i2 - 1];
			outvox[i2 - 1] = imgvox[i1 - 1];
		}
	}
	catrec (imgfile);
	endrec ();

/*******/
/* ifh */
/*******/
	writeifhe (program, outfile, outdim, outvox, orient, control);
/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s", outroot);
	printf ("%s\n", command); status |= system (command);
	exit (status);
}
