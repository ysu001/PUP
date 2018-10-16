/*$Header: /data/petsun4/data1/src_solaris/imato4dfp2/RCS/unpack_4dfp.c,v 1.10 2010/04/03 23:38:59 avi Exp $*/
/*$Log: unpack_4dfp.c,v $
 * Revision 1.10  2010/04/03  23:38:59  avi
 * option -x (squeeze)
 * moderately extensive code rearrangement
 *
 * Revision 1.9  2009/08/22  03:20:50  avi
 * accommodate non-mosaic input
 *
 * Revision 1.8  2008/07/24  01:04:11  avi
 * option -z (reorder slices top to bottom)
 *
 * Revision 1.7  2007/09/19  02:38:50  avi
 * default output endian is input endian
 *
 * Revision 1.6  2007/05/03  19:57:57  avi
 * Solaris 10; endian compliant
 *
 * Revision 1.5  2004/11/15  21:04:48  rsachs
 * Installed 'setprog', 'writeifh'.
 *
 * Revision 1.4  2003/06/26  22:35:21  avi
 * -R option (correct x and y voxel dimensions for pack factor)
 *
 * Revision 1.3  2003/05/06  02:40:59  avi
 * optionally read frame count from ifh slice count
 *
 * Revision 1.2  2003/02/17  06:28:37  avi
 * minor corrections
 *
 * Revision 1.1  2003/02/16  05:14:47  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>

#define MAXL			256
#define NX			64
#define NY			64

/**************/
/* prototypes */
/**************/
extern void	unpackf_ (float *imgt, int *nxp, int *nyp, float *imgu, int *nx, int *ny, int *np);	/* unpack.f */
extern void	flipz (float *imgf, int *pnx, int* pny, int *pnz);		/* cflip.c */

/***********/
/* globals */
/***********/
static char	rcsid[] = "$Id: unpack_4dfp.c,v 1.10 2010/04/03 23:38:59 avi Exp $";
static char	program[MAXL];

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void fsqueeze (float *imgu, int *imgdim, float *imgo, int *outdim) {
	float		*optr;
	int		i, j, sx, sy, ix, iy, iz;
	unsigned long	jndex, m;
	double		q;

	sx = imgdim[0]/outdim[0];
	sy = imgdim[1]/outdim[1];
	optr = imgo;
	for (iz = 0; iz < outdim[2]; iz++) {
	for (iy = 0; iy < outdim[1]; iy++) {
	for (ix = 0; ix < outdim[0]; ix++) {
		jndex = ix*sx + imgdim[0]*(iy*sy + imgdim[1]*iz);
		q = 0.0;
		for (j = 0; j < sy; j++) {
		for (i = 0; i < sx; i++) {
			m = jndex + i + j*imgdim[0];
			q += imgu[m];
		}}
		q /= sx*sy;
		*optr++ = q;
	}}}
}

void usage () {
	printf ("Usage:\t%s <(4dfp) input> <(4dfp) output>\n", program);
	printf (" e.g.,\t%s 030211_EL_b_1 030211_EL_b1\n", program);
	printf ("\toption\n");
	printf ("\t-V\tread frame count from input ifh slice count\n");
	printf ("\t-R\tmultiply output x and y voxsiz by pack factor\n");
	printf ("\t-z\tflipz (unpack slices in reverse order)\n");
	printf ("\t-nx<int>\tspecify unpacked nx (default=%d)\n", NX);
	printf ("\t-ny<int>\tspecify unpacked ny (default=%d)\n", NY);
	printf ("\t-sx<int>\tsqueeze unpacked x dimension by specified factor\n");
	printf ("\t-sy<int>\tsqueeze unpacked y dimension by specified factor\n");
	printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
	exit (1);
}

int main (int argc, char *argv[]) {
/*******/
/* I/O */
/*******/
	FILE			*imgfp, *outfp;
	char			imgroot[MAXL], imgfile[MAXL], ifhfile[MAXL];
	char			outroot[MAXL], outfile[MAXL];
	IFH			imgifh;

/***********/
/* process */
/***********/
	int			imgdim[4], umgdim[4], outdim[4], isbig;
	int			volsiz, slcsiz, iframe, iz, np, nz;
	int			nx = NX, ny = NY, sx = 1, sy = 1, nxo, nyo;
	float			fmin, fmax;
	float			*imgo, *imgu, *imgt, imgvox[3], outvox[3];
	char			control = '\0';

/***********/
/* utility */
/***********/
	char			*ptr, command[MAXL];
	int			c, i, j, k;

/*********/
/* flags */
/*********/
	int			status = 0;
	int			debug = 0;
	int			nz2nv = 0;
	int			mosaic, squeeze;
	int			zflip = 0;
	int			Rflag = 0;	/* correct output voxel dimensions by pack factor */

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'd': debug++;			break;
				case 'R': Rflag++;			break;
				case 'z': zflip++;			break;
				case 'V': nz2nv++;			break;
				case 'n': switch (*ptr++) {
						case 'x': nx = atoi (ptr); break;
						case 'y': ny = atoi (ptr); break;
					  }			*ptr = '\0'; break;
				case 's': switch (*ptr++) {
						case 'x': sx = atoi (ptr); break;
						case 'y': sy = atoi (ptr); break;
					  }			*ptr = '\0'; break;
				case '@': control = *ptr++;	*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	getroot (argv[i], imgroot);	k++; break;
			case 1:	getroot (argv[i], outroot);	k++; break;
		}	
	}
	if (k < 2) usage ();
	if (nx % sx || ny % sy) {
		fprintf (stderr, "%: squeeze factors must evenly divide image dimensions\n");
		usage ();
	}
	nxo = nx/sx; nyo = ny/sy;
	squeeze = nxo != nx || nyo != ny;

/************/
/* read ifh */
/************/
	sprintf (ifhfile, "%s.4dfp.ifh", imgroot);
	fprintf (stdout, "Reading: %s\n", ifhfile);
	if (Getifh (ifhfile, &imgifh)) errr (program, ifhfile);
	isbig = !strcmp (imgifh.imagedata_byte_order, "bigendian");
	if (!control) control = (isbig) ? 'b' : 'l';
	for (k = 0; k < 4; k++) imgdim[k] = imgifh.matrix_size[k];
	for (k = 0; k < 3; k++) imgvox[k] = imgifh.scaling_factor[k];
	printf ("input image dimensions%8d%10d%10d%10d\n", 
			imgdim[0], imgdim[1], imgdim[2], imgdim[3]);
	printf ("voxdim    %10.6f%10.6f%10.6f\n",
			imgvox[0], imgvox[1], imgvox[2]);
	printf ("orientation%9d\n", imgifh.orientation);
	if (imgdim[0] % nx || imgdim[1] % ny) {
		fprintf (stderr, "%s: incommensurate %s mosaic dimensions\n", program, imgroot);
		exit (-1);
	}
	np = (imgdim[0]/nx)*(imgdim[1]/ny);
	mosaic = (np > 1);

	if (nz2nv && mosaic) {
		if (imgdim[3] != 1) {
			fprintf (stderr, "Warning: input ifh volume count not 1\n");
		}
		imgdim[3] = imgdim[2];
	}
	slcsiz = nx*ny;

	volsiz = imgdim[0]*imgdim[1]; if (!mosaic) volsiz *= imgdim[2];
	if (!(imgu = (float *) malloc (volsiz * sizeof (float)))
	||  !(imgt = (float *) malloc (volsiz * sizeof (float)))) errm (program);
	if (squeeze) {
		if (!(imgo = (float *) malloc ((volsiz/(sx*sy)) * sizeof (float)))) errm (program);
	} else {
		imgo = imgu;
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
	printf ("Reading: %s\n", imgfile);
	sprintf (outfile, "%s.4dfp.img", outroot);
	if (!(outfp = fopen (outfile, "wb"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);

	umgdim[0] = nx;
	umgdim[1] = ny;
	umgdim[3] = imgdim[3];

	outdim[0] = nxo;
	outdim[1] = nyo;
	outdim[3] = imgdim[3];
	if (Rflag) {
		for (k = 0; k < 2; k++) outvox[k] = imgvox[k]*imgdim[k]/outdim[k];
	} else {
		outvox[0] = imgvox[0]*sx;
		outvox[1] = imgvox[1]*sy;
		outvox[2] = imgvox[2];
	}

/******************/
/* loop on frames */
/******************/
	fmin = 1.0e37; fmax = -fmin;
	printf ("processing frame");
	nz = volsiz/slcsiz;
	for (iframe = 0; iframe < imgdim[3]; iframe++) {
		if (gread ((char *) imgt, sizeof (float), volsiz, imgfp, isbig)) errr (program, imgfile);
		printf (" %d", iframe + 1); fflush (stdout);
		if (mosaic) {
			unpackf_ (imgt, imgdim+0, imgdim+1, imgu, &nx, &ny, &np);
		} else {
			for (i = 0; i < volsiz; i++) imgu[i] = imgt[i];
		}
/*************************/
/* count nonblank slices */
/*************************/
		for (j = iz = 0; iz < nz; iz++) {
			for (k = i = 0; i < slcsiz; i++) {
				if (imgu[j] > fmax) fmax = imgu[j];
				if (imgu[j] < fmin) fmin = imgu[j];
				k |= (imgu[j++] != 0.0);
			}
			if (!k) break;
		}
		if (iframe) {
			if (iz != outdim[2]) {
				fprintf (stderr, "%s: inconsistent %s mosaic\n", program, imgroot);
				exit (-1);
			}
		} else {
			umgdim[2] = imgdim[2] = outdim[2] = iz;
		}
		if (squeeze) fsqueeze (imgu, umgdim, imgo, outdim);
		if (zflip) flipz (imgo, outdim + 0, outdim + 1, outdim + 2);
		if (ewrite (imgo, outdim[0]*outdim[1]*outdim[2], control, outfp)) errw (program, outfile);
	}
	printf ("\n"); fflush (stdout);
	printf ("data range%10.4f%10.4f\n", fmin, fmax);
	if (fclose (imgfp)) errr (program, imgfile);
	if (fclose (outfp)) errw (program, outfile);

/*******/
/* ifh */
/*******/
	writeifhe (program, outfile, outdim, outvox, imgifh.orientation, control); 
/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s.4dfp.ifh -r%d", outroot, (int) fmax);
	status = system (command);

/************/
/* rec file */
/************/
	sprintf   (outfile, "%s.4dfp.img", outroot);
	startrece (outfile, argc, argv, rcsid, control);
	sprintf   (command, "data range%10.4f%10.4f\n", fmin, fmax); printrec (command);
	if (np == 1) printrec ("non-mosaic condition detected\n");
	if (zflip) printrec ("slices reordered top to bottom\n");
	catrec    (imgfile);
	endrec ();
	sprintf (command, "cat %s.4dfp.ifh", outroot); system (command);

	free (imgt); free (imgu);
	if (squeeze) free (imgo);
        exit (status);
}
