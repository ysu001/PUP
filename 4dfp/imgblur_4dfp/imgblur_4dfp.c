/*$Header: /data/petsun4/data1/src_solaris/imgblur_4dfp/RCS/imgblur_4dfp.c,v 1.9 2009/07/18 01:30:43 avi Exp $*/
/*$Log: imgblur_4dfp.c,v $
 * Revision 1.9  2009/07/18  01:30:43  avi
 * Solaris10/Linux compliant
 *
 * Revision 1.8  2004/09/20  21:39:45  rsachs
 * Installed 'errm','errr','errw','setprog'. Replace 'Get4dfpDimN' with 'get_4dfp_dimo'.
 *
 * Revision 1.7  2000/11/17  02:20:47  avi
 * correct usage
 *
 * Revision 1.6  1999/05/17  05:32:26  avi
 * correct selective xyz blur command line parsing
 *
 * Revision 1.5  1999/01/29  06:42:50  avi
 * correct algorithmic muddle created by tscull
 * Revision 1.4  1999/01/29  06:31:59  avi
 * Revision 1.3  1998/10/12  19:02:52  mcavoy
 * Revision 1.1  1997/01/16  04:46:53  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <unistd.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>

#define MAXL 256

/*************/
/* externals */
/*************/
extern void imgblurx_ (float *fwhm, float *alpha, float *voxsizx, float *img0, int *dimx, int *dimy, int *dimz, float *img1);
extern void imgblury_ (float *fwhm, float *alpha, float *voxsizy, float *img0, int *dimx, int *dimy, int *dimz, float *img1);
extern void imgblurz_ (float *fwhm, float *alpha, float *voxsizz, float *img0, int *dimx, int *dimy, int *dimz, float *img1);

static char rcsid[] = "$Id: imgblur_4dfp.c,v 1.9 2009/07/18 01:30:43 avi Exp $";

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

int main (int argc, char **argv) {
/*************/
/* image I/O */
/*************/
	FILE            *imgfp, *outfp;
	IFH		ifh;
	char            imgfile[MAXL];
	char            outfile[MAXL];
	char            imgroot[MAXL];
	char            outroot[MAXL];
	char		control = '\0';

/**************/
/* processing */
/**************/
	int             imgdim[4], dimension, orient, isbig;
	float           voxdim[3];
	float           *imgt, *imgo;
	float		fwhm, alpha = 0.02;

/***********/
/* utility */
/***********/
        int             c, i, k;
	char            *ptr, command[MAXL], program[MAXL];

/*********/
/* flags */
/*********/
	int             status = 0;
	int		xflag = 0, yflag = 0, zflag = 0, select_flag = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) if (*argv[i] == '-') {
		strcpy (command, argv[i]); ptr = command;
		while (c = *ptr++) switch (c) {
			case 'x': xflag++; select_flag++; 	break;
			case 'y': yflag++; select_flag++; 	break;
			case 'z': zflag++; select_flag++; 	break;
			case '@': control = *ptr++;		*ptr = '\0'; break;
		}
	} else switch (k) {
		case 0:	getroot (argv[i], imgroot);	k++; break;
		case 1: fwhm = atof (argv[i]);		k++; break;
	}
	if (k < 2) {
		printf ("Usage:\t%s [options] <image_file> <FWHM_in_mm>\n", program);
		printf (" e.g.,\t%s -yz vc345 5.5\n", program);
		printf ("\toption\n");
		printf ("\t-x	slective x blur\n");
		printf ("\t-y	slective y blur\n");
		printf ("\t-z	slective z blur\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:	default blur is 3D isotropic\n");
		exit (1);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
/**************************/
/* create output fileroot */
/**************************/
	strcpy (command, "_b");
	if (select_flag) {
		if (xflag) strcat (command, "x");
		if (yflag) strcat (command, "y");
		if (zflag) strcat (command, "z");
	} else {
		xflag  = yflag  = zflag  = 1;
	}
    	sprintf (outroot, "%s%s%d", imgroot, command, (int)(10*fwhm + 0.499999));
	sprintf (outfile, "%s.4dfp.img", outroot);

/*****************************/
/* get input 4dfp dimensions */
/*****************************/
	if (Getifh (imgfile, &ifh)) errr (program, imgfile);
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig)) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';
	dimension = imgdim[0] * imgdim[1] * imgdim[2];
	imgt = (float *) malloc (dimension * sizeof (float));
	imgo = (float *) malloc (dimension * sizeof (float));
	if (!imgt || !imgo) errm (program);

	if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!(outfp = fopen (outfile, "wb"))) errw (program, outfile);
	fprintf (stdout, "Reading: %s\n", imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);
/************/
/* process  */
/************/
	for (i = 0; i < imgdim[3]; i++) {
		if (eread  (imgt, dimension, isbig, imgfp)) errr (program, imgfile);
		if (xflag) {
			imgblurx_ (&fwhm, &alpha, voxdim+0, imgt, imgdim+0, imgdim+1, imgdim+2, imgo);
			for (k = 0; k < dimension; k++) imgt[k] = imgo[k];
		}
		if (yflag) {
			imgblury_ (&fwhm, &alpha, voxdim+1, imgt, imgdim+0, imgdim+1, imgdim+2, imgo);
			for (k = 0; k < dimension; k++) imgt[k] = imgo[k];
		}
		if (zflag) {
			imgblurz_ (&fwhm, &alpha, voxdim+2, imgt, imgdim+0, imgdim+1, imgdim+2, imgo);
			for (k = 0; k < dimension; k++) imgt[k] = imgo[k];
		}
		if (ewrite (imgo, dimension, control, outfp)) errw (program, outfile);
	}
	fclose (imgfp);
	fclose (outfp);

/***************/
/* ifh hdr rec */
/***************/
	if (Writeifh (program, outfile, &ifh, control)) errw (program, outroot);
	sprintf (command, "ifh2hdr %s", outroot);
	status |= system (command);

	startrece (outfile, argc, argv, rcsid, control);
	catrec    (imgfile);
	endrec    ();

	free (imgt);
	free (imgo);
	exit (status);
}
