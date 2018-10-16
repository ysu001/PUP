/*$Header: /data/petsun4/data1/src_solaris/crop_4dfp/RCS/crop_4dfp.c,v 1.10 2010/11/12 05:54:27 avi Exp $*/
/*$Log: crop_4dfp.c,v $
 * Revision 1.10  2010/11/12  05:54:27  avi
 * correct call to errr -> errw in case of write failure
 *
 * Revision 1.9  2010/10/19  03:51:51  avi
 * preseerve mmppix and center in output ifh with option -Z
 *
 * Revision 1.8  2010/08/20  04:59:59  avi
 * correct usage to include option "-@"
 * #include flip_4dfp.h
 *
 * Revision 1.7  2009/04/15  04:28:46  avi
 * option -Z
 *
 * Revision 1.6  2008/05/08  22:07:38  avi
 * correct output ifh
 *
 * Revision 1.5  2008/05/05  04:15:48  avi
 * Solaris 10 and linux compliant
 *
 * Revision 1.4  2004/11/18  21:38:23  rsachs
 * Removed 'Get4dfpDimN','writeifh'. Installed 'get_4dfp_dimo' & call to std. 'writeifh'.
 *
 * Revision 1.3  2004/04/28  04:10:01  avi
 * scroll option
 *
 * Revision 1.2  2004/02/19  01:15:02  avi
 * substitute cflip.o for fflip.o call
 *
 * Revision 1.1  2002/10/19  23:06:09  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>
#include <flip_4dfp.h>

#define MAXL 256

void getirange (char *string, int *minval, int *maxval) {
        char	*str;

	str = strstr (string, "to");
	if (str) {
		*str = '\0';
		*minval = atoi (string);
		*maxval = atoi (str + 2);
	} else {
		*minval = 0;
		*maxval = atoi (string);
	}
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

/********************/
/* global variables */
/********************/
char		program[MAXL];
int		status = 0;
static char	rcsid[] = "$Id: crop_4dfp.c,v 1.10 2010/11/12 05:54:27 avi Exp $";

int	x4dfp2analyze (float *imag, int *dim, int orientation);		/* below */
void	xscroll (float *imgt, float *imgc, int *imgdim, int scrollx);	/* below */
void	yscroll (float *imgt, float *imgc, int *imgdim, int scrolly);	/* below */
void	zscroll (float *imgt, float *imgc, int *imgdim, int scrollz);	/* below */

int main (int argc, char *argv[]) {
/************/
/* 4dfp I/O */
/************/
	FILE		*fp_img, *fp_out;
	IFH		ifh;
	char		outfile[MAXL], imgfile[MAXL];
	char		imgroot[MAXL], outroot[MAXL] = "";
	char		control = '\0';

/****************/
/* image arrays */
/****************/
	float		*img1, *imgc;
	float		voxdim[3];
	int		imgdim[4], outdim[4], ix, iy, iz, index, vdim, isbig;
	int		ix0 = 1, ix1 = 0, iy0 = 1, iy1 = 0, iz0 = 1, iz1 = 0;
	int		vdim_crop, orient;

/***********/
/* utility */
/***********/
	int		c, i, j, k;
	char		*ptr, string[MAXL];

/*********/
/* flags */
/*********/
	int		doflip = 0;
	int		docrop = 1;
	int		scrollx = 0, scrolly = 0, scrollz = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (string, argv[i]); ptr = string;
			while (c = *ptr++) switch (c) {
				case 'f': doflip++;		break;
				case 'Z': docrop = 0;		break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
				case 'x': getirange (ptr, &ix0, &ix1);	*ptr = '\0'; break;
				case 'y': getirange (ptr, &iy0, &iy1);	*ptr = '\0'; break;
				case 'z': getirange (ptr, &iz0, &iz1);	*ptr = '\0'; break;
				case 's': switch (*ptr++) {
					case 'x': scrollx = atoi (ptr); break;
					case 'y': scrolly = atoi (ptr); break;
					case 'z': scrollz = atoi (ptr); break;
				}					*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	getroot (argv[i], imgroot); k++; break;
			case 1:	getroot (argv[i], outroot); k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:\t%s <(4dfp) inroot> [(4dfp) outroot]\n", program);
		printf ("\toption\n");
		printf ("\t-<x|y|z><int>[to[<int>]\tspecify x y z crop limits\n");
		printf ("\t-s<x|y|z><int>\tscroll specified axis by specified number of pixels\n");
		printf ("\t-f\tinterpret specifications under 4dfp<->analyze flips\n");
		printf ("\t-Z\tzero voxels instead of physically cropping\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:\tcrop limit indices count from 1\n");
		printf ("N.B.:\tscrolling is done after cropping\n");
		printf ("N.B.:\tdefault (4dfp) output root is <(4dfp) inroot>\"_crop\"\n");
		exit (1);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0
	||  Getifh (imgfile, &ifh)) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';
	printf ("input  dimensions  %10d%10d%10d%10d\n",   imgdim[0], imgdim[1], imgdim[2], imgdim[3]);
	printf ("voxel  dimensions  %10.6f%10.6f%10.6f\n", voxdim[0], voxdim[1], voxdim[2]);

	if (ix1 > imgdim[0] || ix1 < 1) ix1 = imgdim[0];
	if (iy1 > imgdim[1] || iy1 < 1) iy1 = imgdim[1];
	if (iz1 > imgdim[2] || iz1 < 1) iz1 = imgdim[2];
	if (ix0 > imgdim[0] || ix0 < 1) ix0 = 1;
	if (iy0 > imgdim[1] || iy0 < 1) iy0 = 1;
	if (iz0 > imgdim[2] || iz0 < 1) iz0 = 1;
	if (ix1 < ix0) {k = ix1; ix1 = ix0; ix0 = k;}
	if (iy1 < iy0) {k = iy1; iy1 = iy0; iy0 = k;}
	if (iz1 < iz0) {k = iz1; iz1 = iz0; iz0 = k;}
	printf ("crop limits x %3d to %3d\n", ix0, ix1);
	printf ("crop limits y %3d to %3d\n", iy0, iy1);
	printf ("crop limits z %3d to %3d\n", iz0, iz1);
	printf ("scroll x %3d\n", scrollx);
	printf ("scroll y %3d\n", scrolly);
	printf ("scroll z %3d\n", scrollz);
	outdim[0] = (docrop) ? ix1 - ix0 + 1 : imgdim[0];
	outdim[1] = (docrop) ? iy1 - iy0 + 1 : imgdim[1];
	outdim[2] = (docrop) ? iz1 - iz0 + 1 : imgdim[2];
	outdim[3] = imgdim[3];
	printf ("output dimensions  %10d%10d%10d%10d\n",   outdim[0], outdim[1], outdim[2], outdim[3]);

	vdim =  imgdim[0] * imgdim[1] * imgdim[2];
	if (!(img1 = (float *) malloc (vdim * sizeof (float)))) errm (program);
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	printf ("Reading: %s\n", imgfile);

	if (!strlen (outroot)) sprintf (outroot, "%s_crop", imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);
	vdim_crop = outdim[0] * outdim[1] * outdim[2];
	if (!(imgc = (float *) calloc (vdim_crop, sizeof (float)))) errm (program);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);

/***********/
/* process */
/***********/
	printf ("processing volume");
	for (k = 0; k < imgdim[3]; k++) {printf (" %d", k + 1); fflush (stdout);
		if (eread (img1, vdim, isbig, fp_img)) errr (program, imgfile);
		if (doflip) {
			if (x4dfp2analyze (img1, imgdim, ifh.orientation)) {
				fprintf (stderr, "%s: invalid %s orientation\n", program, imgroot);
				exit (-1);
			}
		}
		i = 0;
		for (iz = iz0 - 1; iz < iz1; iz++) {
		for (iy = iy0 - 1; iy < iy1; iy++) {
		for (ix = ix0 - 1; ix < ix1; ix++) {
			index = ix + imgdim[0]*(iy + imgdim[1]*iz);
			if (docrop) {
				imgc[i++]	= img1[index];
			} else {
				imgc[index]	= img1[index];
			}
		}}}
		if (scrollx) xscroll (img1, imgc, outdim, scrollx);
		if (scrolly) yscroll (img1, imgc, outdim, scrolly);
		if (scrollz) zscroll (img1, imgc, outdim, scrollz);
		if (doflip) x4dfp2analyze (imgc, outdim, ifh.orientation);
		if (ewrite (imgc, vdim_crop, control, fp_out)) errw (program, outfile);
	} printf ("\n"); fflush (stdout);
	fclose (fp_img);
	fclose (fp_out);

/**************/
/* output ifh */
/**************/
	if (docrop) {
		if (writeifhe (program, outroot, outdim, voxdim, ifh.orientation, control)) errw (program, outroot);
	} else {
		Writeifh (program, outroot, &ifh, control);
	}

/***************/
/* run ifh2hdr */
/***************/
	sprintf (string, "ifh2hdr %s", outroot); status |= system (string);

/************/
/* rec file */
/************/
	startrece (outfile, argc, argv, rcsid, control);
	if (doflip)  printrec ("indices subjected to 4dfp<->analyze orientation dependent flips\n");
	if (!docrop) printrec ("voxels zeroed instead of physically cropped\n");
	sprintf (string, "crop limits x %3d to %3d\n", ix0, ix1); printrec (string);
	sprintf (string, "crop limits y %3d to %3d\n", iy0, iy1); printrec (string);
	sprintf (string, "crop limits z %3d to %3d\n", iz0, iz1); printrec (string);
	if (scrollx) {sprintf (string, "scroll x %3d\n", scrollx); printrec (string);}
	if (scrolly) {sprintf (string, "scroll y %3d\n", scrolly); printrec (string);}
	if (scrollz) {sprintf (string, "scroll z %3d\n", scrollz); printrec (string);}
	catrec (imgfile);
	endrec ();

/************/
/* clean up */
/************/
	free (img1); free (imgc);
	exit (status);
}

int x4dfp2analyze (float *imag, int *dim, int orientation) {
	switch (orientation) {
		case 4: flipx (imag, dim+0, dim+1, dim+2);	/* sagittal */
		case 3:	flipz (imag, dim+0, dim+1, dim+2);	/* coronal */
		case 2:	flipy (imag, dim+0, dim+1, dim+2);	/* transverse */
			return 0;
		default: return -1;				/* none of the above */
	}
}

void xscroll (float *imgt, float *imgc, int *imgdim, int scrollx) {
	int	i, ix, ixp, iy, iz, index;

	while (scrollx < 0) scrollx += imgdim[0];
	i = 0;
	for (iz = 0; iz < imgdim[2]; iz++) {
	for (iy = 0; iy < imgdim[1]; iy++) {
	for (ix = 0; ix < imgdim[0]; ix++) {
		ixp = (ix - scrollx + imgdim[0]) % imgdim[0];
		index = ixp + imgdim[0]*(iy + imgdim[1]*iz);
		imgt[i++] = imgc[index];
	}}}
	for (i = 0; i < imgdim[0]*imgdim[1]*imgdim[2]; i++) imgc[i] = imgt[i];
}

void yscroll (float *imgt, float *imgc, int *imgdim, int scrolly) {
	int	i, ix, iy, iyp, iz, index;

	while (scrolly < 0) scrolly += imgdim[1];
	i = 0;
	for (iz = 0; iz < imgdim[2]; iz++) {
	for (iy = 0; iy < imgdim[1]; iy++) {
	for (ix = 0; ix < imgdim[0]; ix++) {
		iyp = (iy - scrolly + imgdim[1]) % imgdim[1];
		index = ix + imgdim[0]*(iyp + imgdim[1]*iz);
		imgt[i++] = imgc[index];
	}}}
	for (i = 0; i < imgdim[0]*imgdim[1]*imgdim[2]; i++) imgc[i] = imgt[i];
}

void zscroll (float *imgt, float *imgc, int *imgdim, int scrollz) {
	int	i, ix, iy, iz, izp, index;

	while (scrollz < 0) scrollz += imgdim[2];
	i = 0;
	for (iz = 0; iz < imgdim[2]; iz++) {
	for (iy = 0; iy < imgdim[1]; iy++) {
	for (ix = 0; ix < imgdim[0]; ix++) {
		izp = (iz - scrollz + imgdim[2]) % imgdim[2];
		index = ix + imgdim[0]*(iy + imgdim[1]*izp);
		imgt[i++] = imgc[index];
	}}}
	for (i = 0; i < imgdim[0]*imgdim[1]*imgdim[2]; i++) imgc[i] = imgt[i];
}
