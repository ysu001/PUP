/*$Header: /data/petsun4/data1/src_solaris/S2T_4dfp/RCS/C2T_4dfp.c,v 1.8 2007/08/31 04:49:31 avi Exp $*/
/*$Log: C2T_4dfp.c,v $
 * Revision 1.8  2007/08/31  04:49:31  avi
 * linux and X86 Solaris compliant
 *
 * Revision 1.7  2007/01/18  06:52:56  avi
 * explicit ifh not needed
 *
 * Revision 1.6  2007/01/18  03:32:16  avi
 * correct voxdimo[] code
 *
 * Revision 1.5  2007/01/18  03:24:35  avi
 * output correctly ordered voxdim
 *
 * Revision 1.4  2004/12/03  22:01:18  rsachs
 * Corrected the output ifh orientation.
 *
 * Revision 1.3  2004/11/23  20:51:02  rsachs
 * Replaced code with call to 'writeifh'. Installed 'setprog'.
 *
 * Revision 1.2  2004/03/11  05:34:56  avi
 * eliminate -lrms dependencies
 *
 * Revision 1.1  2003/01/08  06:35:06  avi
 * Initial revision
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <unistd.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MAXL 256

/*************/
/* externals */
/*************/
extern void	c2tf_ (float *imgm, int *nx, int *nz, int *ny, float *imgn);	/* reorder.f */
extern void	flipz (float *imgf, int *pnx, int* pny, int *pnz);		/* cflip.c */

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[]= "$Id: C2T_4dfp.c,v 1.8 2007/08/31 04:49:31 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	FILE            *fp_img, *fp_out;
	char            imgfile[MAXL], outfile[MAXL];
	char            imgroot[MAXL], outroot[MAXL] = "";
	char            *ptr, command[MAXL], program[MAXL];

/**************/
/* processing */
/**************/
	int             imgdim[4], imgdimo[4], orient, isbig;
	float           voxdim[3], voxdimo[3];
	float           *img1, *img2;
	int             dimension;
	char		control = '\0';

/***********/
/* utility */
/***********/
	int             c, i, k;
	int		status;

	printf ("%s\n", rcsid);
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
		} else switch (k) {
			case 0: getroot (argv[i], imgroot);	k++; break;
			case 1: getroot (argv[i], outroot);	k++; break;
		}
	}
	if (k < 1) {
		printf ("Usage:\t%s <(4dfp) image> [(4dfp) outroot]\n", program);
		printf (" e.g.,\t%s vm6c_b1\n", program);
		printf (" e.g.,\t%s vm6c_b1 vm6c_b1T\n", program);
		printf ("\toption\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:	default output root = <imgroot>\"T\"\n");
		exit (1);
	}
	if (!strlen (outroot)) sprintf (outroot, "%sT", imgroot); 

/*****************************/
/* get 4dfp input dimensions */
/*****************************/
	if (get_4dfp_dimoe (imgroot, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgroot);
	if (!control) control = (isbig) ? 'b' : 'l';
	if (orient != 3) fprintf (stderr, "%s warning: %s orientation (%d) not coronal\n", program, imgroot, orient);
	dimension = imgdim[0] * imgdim[1] * imgdim[2];

/*****************/
/* alloc buffers */
/*****************/
	if (!(img1 = (float *) malloc (dimension * sizeof (float)))
	||  !(img2 = (float *) malloc (dimension * sizeof (float)))) errm (program);

/************/
/* process  */
/************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	fprintf (stdout, "Reading: %s\n", imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);
	for (k = 0; k < imgdim[3]; k++) {
		if (eread (img1, dimension, isbig, fp_img)) errr (program, imgfile);
		c2tf_ (img1, imgdim + 0, imgdim + 1, imgdim + 2, img2);
		flipz (img2, imgdim + 0, imgdim + 2, imgdim + 1);
		if (ewrite (img2, dimension, control, fp_out)) errw (program, outfile);
	}
	fclose (fp_img);
	fclose (fp_out);

/*******************/
/* create ifh file */
/*******************/
	imgdimo[0] = imgdim[0]; voxdimo[0] = voxdim[0];
	imgdimo[1] = imgdim[2]; voxdimo[1] = voxdim[2];
	imgdimo[2] = imgdim[1]; voxdimo[2] = voxdim[1];
	imgdimo[3] = imgdim[3];
	writeifhe (program, outfile, imgdimo, voxdimo, 2, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s", outroot); printf ("%s\n", command);
	status = system (command);

/*******************/
/* create rec file */
/*******************/
	sprintf   (outfile, "%s.4dfp.img", outroot);
        startrece (outfile, argc, argv, rcsid, control);
        catrec    (imgfile);
        endrec    ();

        free (img1); free (img2);
        exit (status);
}
