/*$Header: /data/petsun4/data1/src_solaris/peak_4dfp/RCS/burn_sphere_4dfp.c,v 1.13 2010/11/08 23:51:53 avi Exp $*/
/*$Log: burn_sphere_4dfp.c,v $
 * Revision 1.13  2010/11/08  23:51:53  avi
 * extend -s functionality to single voxels (orad <= 0)
 *
 * Revision 1.12  2010/11/02  23:41:39  avi
 * option -s
 *
 * Revision 1.11  2010/08/12  03:23:45  avi
 * eliminate limit on input list length
 *
 * Revision 1.10  2009/07/06  05:52:51  avi
 * greatly accelerate by search within computed bounding box
 *
 * Revision 1.9  2009/07/06  03:56:44  avi
 * option -l
 *
 * Revision 1.8  2007/09/23  03:15:22  avi
 * linux compliant
 *
 * Revision 1.7  2007/06/06  22:05:15  avi
 * usage update
 *
 * Revision 1.6  2007/06/06  21:53:12  mohanar
 * Fixed endian issue
 *
 * Revision 1.5  2006/09/26  02:36:15  avi
 * Solaris 10
 *
 * Revision 1.4  2006/08/05  01:49:17  avi
 * enable standard altas space specification of output
 *
 * Revision 1.3  2006/03/07  04:14:59  avi
 * single voxel burn
 *
 * Revision 1.2  2005/05/23  05:21:09  avi
 * first 3 arguments cannot be options
 *
 * Revision 1.1  2005/05/23  04:32:34  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>	/* R_OK */
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>

#define MAXL		256
#define BVAL		1.
#define ORAD		6.

/*************/
/* externals */
/*************/
void vrtflip_ 	(int *orientation, int *imgdim, float *center, float *mmppix, float *centerr, float *mmppixr);

/********************/
/* global variables */
/********************/
static	char		program[MAXL];
static	char		rcsid[] = "$Id: burn_sphere_4dfp.c,v 1.13 2010/11/08 23:51:53 avi Exp $";

float **calloc_float2 (int n1, int n2) {
	int	i;
	float	**a;

	if (!(a = (float **) malloc (n1 * sizeof (float *)))) errm (program);
	if (!(a[0] = (float *) calloc (n1 * n2, sizeof (float)))) errm (program);
	for (i = 1; i < n1; i++) a[i] = a[0] + i*n2;
	return a;
}

void free_float2 (float **a) {
	free (a[0]);
	free (a);
}

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

int main (int argc, char *argv[]) {
	FILE		*fp;			/* utility */
	FILE		*imgfp;			/* input image */
	FILE		*outfp;			/* output image */

/*************/
/* image i/o */
/*************/
	char		*ptr, command[MAXL], *srgv[MAXL];
	char		imgroot[MAXL], imgfile[MAXL], ifhfile[MAXL], lstfile[MAXL] = "";
	char		outroot[MAXL], outfile[MAXL];

/***************/
/* computation */
/***************/
	float		bval = BVAL, orad = 6.;
	float		*imgv;
	IFH		imgifh;
	float		mmppixr[3], centerr[3];
	float		fndex[3], x[3];
	float		**x0, x1[3];
	int		nlst = 0;
	int		ix, iy, iz, dim, nx, ny, nz, yshift = 0, isbig;
	int		ixl, ixh, iyl, iyh, izl, izh;
	int		c, i, j, k, m;
	double		q, dist;
	char		control = '\0';

/*********/
/* flags */
/*********/
	int		read = 0;
	int		status = 0;
	int		debug = 0;
	int		super = 0;
	int		sum_flag = 0;

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);

/******************************/
/* get command line arguments */
/******************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-' && k > 2) {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'a': super++;		break;
				case 's': sum_flag++;		break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
				case 'o': orad = atof (ptr);		*ptr = '\0'; break;
				case 'l': strcpy (lstfile, ptr);	*ptr = '\0'; break;
				case 'v': bval = atof (ptr);		*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0: case 1: case 2: x1[k] = atof (argv[i]);		k++; break;
		 	case 3: getroot (argv[i], imgroot);			k++; break;
		 	case 4: getroot (argv[i], outroot);			k++; break;
		}
	}
	if (k < 5) {
		printf ("Usage:	%s <flt x0> <flt y0> <flt z0> <4dfp imgroot> <4dfp outroot> [options]\n", program);
		printf ("e.g.:	%s 33.1 -56.2 18. grand_average_222[.4dfp.img] -v2 -o7.5\n", program);
		printf ("or:	%s 33.1 -56.2 18. 222                          -v2 -o7.5\n", program);
		printf ("\toption\n");
		printf ("\t-a\tsuperimpose sphere on image (default duplicate input format with zero background)\n");
		printf ("\t-s\tsum overlapping spheres (default overwrite)\n");
		printf ("\t-v<flt>\tspecify burn in value (default=%.4f)\n", BVAL);
		printf ("\t-o<flt>\tspecify sphere radius in mm (default=%.4f)\n", ORAD);
		printf ("\t-l<lst>\tread sphere coordinates from specified list (command line coords ignored)\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:\tsetting output radius to zero causes a single pixel burn\n");
		printf ("N.B.:\twithout -a only the input ifh (or standard atlas string) is required\n", ORAD);
		printf ("N.B.:\tspecifying <4dfp imgroot> as \"333[.n]\" \"222\" or \"111\" generates standard atlas space output\n");
		printf ("N.B.:\tif the 4dfp image does not exist the default output endianness is CPU endian\n");
		exit (1);
	}

/*******************/
/* parse list file */
/*******************/
	if (strlen (lstfile)) {
		if (!(fp = fopen (lstfile, "r"))) errr (program, lstfile);
		nlst = 0; while (fgets (command, MAXL, fp)) {
			m = split (command, srgv, MAXL); 
			if (!m) continue;				/* skip blank lines */
			if (m == 3) {
				nlst++;
			} else {
				fprintf (stderr, "%s: %s line ignored\n", program, lstfile);
			}
		}
		x0 = calloc_float2 (nlst, 3);
		rewind (fp);
		j = 0; while (fgets (command, MAXL, fp)) {
			m = split (command, srgv, MAXL); 
			if (!m || m != 3) continue;
			for (k = 0; k < 3; k++) x0[j][k] = atof (srgv[k]);
			j++;
		}
		fclose (fp);
	} else {
		nlst = 1;
		x0 = calloc_float2 (nlst, 3);
		for (k = 0; k < 3; k++) x0[0][k] = x1[k];
	}

	imgifh.orientation = 2; isbig = CPU_is_bigendian() ;
	if (0) {
	} else if (!strncmp (imgroot, "333", 3)) {	/* see t4imgs_4dfp.c */
		if (imgroot[3] == '.') yshift = atoi (imgroot + 4);
		imgifh.matrix_size[0] = 48;	imgifh.matrix_size[1] = 64;		imgifh.matrix_size[2] = 48;
		imgifh.scaling_factor[0] = 3.0;	imgifh.scaling_factor[1] = 3.0;		imgifh.scaling_factor[2] = 3.0;
		imgifh.mmppix[0] = 3.0;		imgifh.mmppix[1] = -3.0;		imgifh.mmppix[2] = -3.0;
		imgifh.center[0] = 73.5;	imgifh.center[1] = -3.0*(29 - yshift); 	imgifh.center[2] = -84.0;
	} else if (!strcmp (imgroot, "222")) {		/* see t4imgs_4dfp.c */
		imgifh.matrix_size[0] = 128;	imgifh.matrix_size[1] = 128;		imgifh.matrix_size[2] = 75;
		imgifh.scaling_factor[0] = 2.0;	imgifh.scaling_factor[1] = 2.0;		imgifh.scaling_factor[2] = 2.0;
		imgifh.mmppix[0] = 2.0;		imgifh.mmppix[1] = -2.0;		imgifh.mmppix[2] = -2.0;
		imgifh.center[0] = 129.;	imgifh.center[1] = -129.; 		imgifh.center[2] = -2.0*41;
	} else if (!strcmp (imgroot, "111")) {		/* see t4imgs_4dfp.c */
		imgifh.matrix_size[0] = 176;	imgifh.matrix_size[1] = 208;		imgifh.matrix_size[2] = 176;
		imgifh.scaling_factor[0] = 1.0;	imgifh.scaling_factor[1] = 1.0;		imgifh.scaling_factor[2] = 1.0;
		imgifh.mmppix[0] = 1.0;		imgifh.mmppix[1] = -1.0;		imgifh.mmppix[2] = -1.0;
		imgifh.center[0] = 89.;		imgifh.center[1] = -85.; 		imgifh.center[2] = -101.;
	} else {
/************/
/* read ifh */
/************/
		sprintf (ifhfile, "%s.4dfp.img", imgroot);
		printf ("Reading: %s\n", ifhfile);
		if (Getifh (ifhfile, &imgifh)) errr (program, ifhfile);
		isbig = strcmp (imgifh.imagedata_byte_order, "littleendian");
		read++;
	}
	printf ("image dimensions \t%10d%10d%10d\n",
		imgifh.matrix_size[0], imgifh.matrix_size[1], imgifh.matrix_size[2]);
	printf ("image mmppix     \t%10.6f%10.6f%10.6f\n",
		imgifh.mmppix[0], imgifh.mmppix[1], imgifh.mmppix[2]);
	printf ("image center     \t%10.4f%10.4f%10.4f\n",
		imgifh.center[0], imgifh.center[1], imgifh.center[2]);
	printf ("image orientation\t%10d\n", imgifh.orientation);
	if (!control) control = (isbig) ? 'b' : 'l';

/*****************************************/
/* virtual flip instead of x4dfp2ecat () */
/*****************************************/
	vrtflip_ (&imgifh.orientation, imgifh.matrix_size, imgifh.center, imgifh.mmppix, centerr, mmppixr);

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	dim = 1; for (k = 0; k < 3; k++) dim *= imgifh.matrix_size[k];
	if (!(imgv = (float *) calloc (dim, sizeof (float)))) errm (program);
	if (super) {
		printf ("Reading: %s\n", imgfile);
		if (!(imgfp = fopen (imgfile, "rb")) || eread (imgv, dim, isbig, imgfp)
		|| fclose (imgfp)) errr (program, imgfile);
	}

/******************/
/* burn in sphere */
/******************/
	nx = imgifh.matrix_size[0];
	ny = imgifh.matrix_size[1];
	nz = imgifh.matrix_size[2];
	for (j = 0; j < nlst; j++) {
		for (k = 0; k < 3; k++) fndex[k] = -0.5 + (x0[j][k] + centerr[k])/mmppixr[k];
		if (orad <= 0) {
			ix = fndex[0]; iy = fndex[1]; iz = fndex[2];
			if (fndex[0] >= 0 && ix < nx && fndex[1] >= 0 && iy < ny && fndex[2] >= 0 && iz < nz) {
				i = ix + nx*(iy + ny*iz);
				if (sum_flag) {
					imgv[i] += bval;
				} else {
					imgv[i]  = bval;
				}
			}
		} else {
			ixl = fndex[0] - (0.5 + orad/fabs(mmppixr[0])); if (ixl < 0)  ixl = 0;
			ixh = fndex[0] + (0.5 + orad/fabs(mmppixr[0])); if (ixh > nx) ixh = nx;
			iyl = fndex[1] - (0.5 + orad/fabs(mmppixr[1])); if (iyl < 0)  iyl = 0;
			iyh = fndex[1] + (0.5 + orad/fabs(mmppixr[1])); if (iyh > ny) iyh = ny;
			izl = fndex[2] - (0.5 + orad/fabs(mmppixr[2])); if (izl < 0)  izl = 0;
			izh = fndex[2] + (0.5 + orad/fabs(mmppixr[2])); if (izh > nz) izh = nz;
			if (0) printf (" ixl, ixh, iyl, iyh, izl, izh, %d %d %d %d %d %d\n", ixl, ixh, iyl, iyh, izl, izh);
			for (iz = izl; iz < izh; iz++) {
			for (iy = iyl; iy < iyh; iy++) {
			for (ix = ixl; ix < ixh; ix++) {	
				i = ix + nx*(iy + ny*iz);
				fndex[0] = (float) (ix + 1);
				fndex[1] = (float) (iy + 1);
				fndex[2] = (float) (iz + 1);
				for (q = k = 0; k < 3; k++) {
					x[k] = fndex[k]*mmppixr[k] - centerr[k];
					q += (x[k] - x0[j][k])*(x[k] - x0[j][k]);
				}
				dist = sqrt (q);
				if (dist < orad) {
					if (sum_flag) {
						imgv[i] += bval;
					} else {
						imgv[i]  = bval;
					}
					if (0) printf ("ix, iy, iz, %d %d %d\n", ix, iy, iz);
				}
			}}}
		}
	}
	sprintf (outfile, "%s.4dfp.img", outroot);
	printf ("Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "wb")) || ewrite (imgv, dim, control, outfp)
	|| fclose (outfp)) errw (program, outfile);

/*******/
/* rec */
/*******/
	startrece (outfile, argc, argv, rcsid, control);
	if (strlen (lstfile)) {
		printrec ("start coord list\n");
		catrec (lstfile);
		printrec ("end coord list\n");
	} else {
		sprintf (command, "%.4f mm radius sphere burned in at coords %10.4f %10.4f %10.4f\n",
			orad, x0[0][0], x0[0][1], x0[0][2]); printrec (command);
	}
	if (read) catrec (imgfile);
	endrec ();

/*******/
/* ifh */
/*******/
	imgifh.matrix_size[3] = 1;
	writeifhmce (program, outfile, imgifh.matrix_size, imgifh.scaling_factor, imgifh.orientation,
		imgifh.mmppix, imgifh.center, control);

/***************/
/* analyze hdr */
/***************/
	if (super) {
		sprintf (command, "ifh2hdr %s", outroot);
	} else {
		sprintf (command, "ifh2hdr %s -r%d", outroot, (int) bval);
	}
	printf ("%s\n", command); status |= system (command);

	free (imgv);
	free_float2 (x0);
	exit (status);
}
