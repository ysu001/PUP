/*$Header: /data/petsun4/data1/src_solaris/gauss_4dfp/RCS/hsphere_4dfp.c,v 1.5 2007/04/17 05:26:50 avi Exp $*/
/*$Log: hsphere_4dfp.c,v $
 * Revision 1.5  2007/04/17  05:26:50  avi
 * gcc compliant
 *
 * Revision 1.4  2006/09/25  17:43:17  avi
 * Solaris 10
 *
 * Revision 1.3  2004/10/11  22:00:29  rsachs
 * Installed 'errm','errr','errw','getroot','get_4dfp_dimo'. Removed 'Get4dfpDimN'.
 *
 * Revision 1.2  2001/07/11  06:31:36  avi
 * correct margin computation
 *
 * Revision 1.1  2001/07/11  04:31:18  avi
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
#include <librms.h>

#define MAXL      256	

/*************/
/* externals */
/*************/
extern void	hsphere3d_ (float *imag, int *nx, int *ny, int *nz, float *mmppix, float *radius);		/* hsphere3d.f */

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[] = "$Id: hsphere_4dfp.c,v 1.5 2007/04/17 05:26:50 avi Exp $";
int main (int argc, char **argv) {
	FILE 		*fp_img, *fp_out;
	IFH		ifh;
	char            imgroot[MAXL], imgfile[MAXL];
	char            outroot[MAXL], outfile[MAXL];

        int  		imgdim[4];
        float	 	voxdim[3];	
	float		radius;
	float		*image3d;
	float		*padimage3d;
	int		nx, ny, nz;
	int		nxp, nyp, nzp;
	int		dimension, margin, orient, isbig;
	char		control = '\0';

/***********/
/* utility */
/***********/
	char 		*ptr, program[MAXL], command[MAXL];
	float		val;
	int		c, i, k;

/*********/
/* flags */
/*********/
	int		debug = 0;
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
			}
		}
		else switch (k) {
			case 0:	getroot (argv[i], imgroot);	k++; break;
			case 1:	radius = atof (argv[i]);	k++; break;
		}	
	}
	if (k < 2) {
		printf ("Usage:\t%s <file_4dfp> <flt>\n", program);
		printf (" e.g.,\t%s np1234_zmap_222 5\n", program);
		printf ("N.B.:	second argument specifies hard sphere radius in mm\n");
		printf ("\toptions\n");
		exit (1);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	sprintf (outroot, "%s_%.0fmm", imgroot, radius);
	sprintf (outfile, "%s.4dfp.img", outroot);

/***********************/
/* get 4dfp dimensions */
/***********************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';
	if (Getifh (imgfile, &ifh)) errr (program, imgfile);
	nx = imgdim[0];
	ny = imgdim[1];
	nz = imgdim[2];
	dimension = nx * ny * nz;

/********************/
/* allocate buffers */
/********************/
	margin = 2.0 * radius / voxdim[0]; nxp = npad_ (&nx, &margin);
	margin = 2.0 * radius / voxdim[1]; nyp = npad_ (&ny, &margin);
	margin = 2.0 * radius / voxdim[2]; nzp = npad_ (&nz, &margin);
	printf ("image dimensions %d %d %d padded to %d %d %d \n", nx, ny, nz, nxp, nyp, nzp);
        image3d = (float *) malloc (dimension * sizeof (float));
        padimage3d = (float *) calloc (nxp * nyp * nzp, sizeof (float));
        if (!padimage3d || !image3d) errm (program);

	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	printf ("Reading: %s\n", imgfile);
	printf ("Writing: %s\n", outfile);
	for (k = 0; k < imgdim[3]; k++) {
		if (eread (image3d, dimension, isbig, fp_img)) errr (program, imgfile);
		imgpad_		(image3d, &nx, &ny, &nz, padimage3d, &nxp, &nyp, &nzp);
		hsphere3d_	(padimage3d, &nxp, &nyp, &nzp, voxdim, &radius);
		imgdap_		(image3d, &nx, &ny, &nz, padimage3d, &nxp, &nyp, &nzp);
		if (ewrite (image3d, dimension, control, fp_out)) errw (program, outfile);
	}
	fclose (fp_img);
	fclose (fp_out);

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, rcsid, control);
	sprintf   (command, "Convolved with hard sphere of radius %.4f mm\n", radius); printrec (command);
	catrec    (imgfile);
	endrec    ();

/***************/
/* ifh and hdr */
/***************/
	if (Writeifh (program, outfile, &ifh, control)) errw (program, outroot);
	sprintf (command, "ifh2hdr %s", outroot);
	status |= system (command);

	free (padimage3d);
	free (image3d);
	exit (status);
}
