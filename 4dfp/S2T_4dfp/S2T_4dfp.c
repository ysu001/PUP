/*$Header: /data/petsun4/data1/src_solaris/S2T_4dfp/RCS/S2T_4dfp.c,v 1.6 2007/08/31 05:03:24 avi Exp $*/
/*$Log: S2T_4dfp.c,v $
 * Revision 1.6  2007/08/31  05:03:24  avi
 * linux and X86 Solaris compliant
 *
 * Revision 1.5  2007/01/18  06:53:37  avi
 * use endianio
 *
 * Revision 1.4  2004/12/03  21:58:36  rsachs
 * Corrected the orientation in the ifh file.
 *
 * Revision 1.3  2004/11/23  21:18:23  rsachs
 * Replaced code with call to 'writeifh'. Installed 'setprog'.
 *
 * Revision 1.2  2004/03/11  08:06:49  avi
 * correct orientatation warning
 *
 * Revision 1.1  2004/03/11  06:08:49  avi
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
extern void	s2tf_ (float *imgm, int *nx, int *nz, int *ny, float *imgn);	/* reorder.f */
extern void	flipx (float *imgf, int *pnx, int* pny, int *pnz);		/* cflip.c */
extern void	flipz (float *imgf, int *pnx, int* pny, int *pnz);		/* cflip.c */

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[]= "$Id: S2T_4dfp.c,v 1.6 2007/08/31 05:03:24 avi Exp $";
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
		printf ("Usage:\t%s <(4dfp) imgroot> [(4dfp) outroot]\n", program);
		printf (" e.g.,\t%s vm6c_mpr\n", program);
		printf (" e.g.,\t%s vm6c_mpr vm6c_mprT\n", program);
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
	if (orient != 4) fprintf (stderr, "%s warning: %s orientation (%d) not sagittal\n", program, imgroot, orient);
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
		s2tf_ (img1, imgdim + 0, imgdim + 1, imgdim + 2, img2);
		flipx (img2, imgdim + 2, imgdim + 0, imgdim + 1);
		flipz (img2, imgdim + 2, imgdim + 0, imgdim + 1);
		if (ewrite (img2, dimension, control, fp_out)) errw (program, outfile);
	}
	fclose (fp_img);
	fclose (fp_out);

/*******************/
/* create ifh file */
/*******************/
	imgdimo[0] = imgdim[2]; voxdimo[0] = voxdim[2];
	imgdimo[1] = imgdim[0]; voxdimo[1] = voxdim[0];
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
