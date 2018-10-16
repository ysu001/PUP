/*$Log: t4imgs_4dfp.c,v $
 * Revision 1.36  2011/04/09  00:19:16  avi
 * input list line buffer length increased from 2*MAXL to 4*MAXL
 * filename buffer lenght increased form MAXL to 2*MAXL
 *
 * Revision 1.35  2010/12/23  03:35:43  avi
 * long -> int in rnan()
 *
 * Revision 1.34  2010/04/05  02:46:03  avi
 * OK mask and img do not match endianness
 *
 * Revision 1.33  2010/04/05  02:12:46  avi
 * option -R (suppress saving rec file)
 *
 * Revision 1.32  2010/03/17  23:03:10  avi
 * correct bug applying dimension check to volume count of mask
 *
 * Revision 1.31  2009/07/21  04:34:03  avi
 * eliminate (time) shift option
 * maintain imgs and imgw accumulators on disk to reduce memory usage during
 * multiframe transforms
 *
 * Revision 1.30  2008/01/15  06:05:35  avi
 * correct bug arising when dimroot == "333*"
 *
 * Revision 1.29  2007/08/10  02:58:40  avi
 * eliminate FORTRAN f_exit() and f_init()
 *
 * Revision 1.28  2007/05/04  01:49:10  avi
 * quiet gcc -Wall
 *
 * Revision 1.27  2007/05/01  05:10:26  avi
 * restore lost endrec
 *
 * Revision 1.26  2007/04/29  05:17:51  avi
 * gcc v3 compliant: eliminate limits.h sunmath.h; define FLT_MAX and float rnan();
 * screen for isnan before performing numeric comparisons
 *
 * Revision 1.25  2006/09/26  23:26:58  avi
 * Solaris 10
 *
 * Revision 1.24  2006/02/13  07:34:00  avi
 * eliminate dependency on libmri
 * create analyze hdr
 *
 * Revision 1.23  2004/09/02  20:19:46  rsachs
 * Installed fns 'errm','errr','errw' & 'getroot'.
 *
 * Revision 1.22  2004/08/23  23:57:16  rsachs
 * Added nearest-neighbor interpolation.
 *
 * Revision 1.21  2001/01/10  00:49:18  avi
 * correct file_error ()
 *
 * Revision 1.20  2001/01/07  05:03:01  avi
 * -s option (3D spline interpolation)
 *
 * Revision 1.19  2001/01/07  00:23:49  avi
 * define -O to yield 111h convention
 *
 * Revision 1.18  1999/10/31  04:49:43  avi
 * mask= and slices= list file options
 *
 * Revision 1.17  1999/10/20  21:01:11  avi
 * NSTACKMAX -> no limit
 *
 * Revision 1.16  1999/10/20  04:17:37  avi
 * correct misidentification (imgreg_4dfp) in file_error ()
 *
 * Revision 1.15  1999/02/10  07:57:32  avi
 * NaN_flag
 *
 * Revision 1.14  1999/02/06  01:05:13  avi
 * more secure dimroot code
 * one-field input list lines allowed
 *
 * Revision 1.13  1999/01/28  22:23:39  avi
 * make sqrt_flag actually have an effect
 *
 * Revision 1.12  1999/01/10  04:07:50  avi
 * t4file2param_ () -> t4scale ()
 *
 * Revision 1.11  1999/01/10  02:36:05  avi
 * call t4file2param to get t4file scale
 *
 * Revision 1.10  1999/01/04  01:48:24  avi
 * call vrtflip_ on input images instead of x4dfp2ecat
 *
 * Revision 1.9  1998/12/31  06:42:37  avi
 * log to screen mmppix and center as in input and out ifh files
 *
 * Revision 1.8  1998/12/31  02:11:11  avi
 * -O111 and -O222 switches must match exactly
 *
 * Revision 1.7  1998/12/30  00:16:51  avi
 * x4dfp2ecat (,, orio) also after transform
 *
 * Revision 1.6  1998/12/29  23:23:21  avi
 * allow output orientation to be not 2
 *
 * Revision 1.5  1998/12/29  07:52:15  avi
 * -O111 -O222 -O333.n -Omy_file
 *
 * Revision 1.4  1998/12/25  08:40:00  avi
 * move code closer to imgreg_4dfp.c without yet adding new functionality.
 *
 * Revision 1.3  1998/06/28  22:06:48  avi
 * -2 option
 *
 * Revision 1.2  1997/10/11  07:04:03  avi
 * -B option (comwrp=to_711-2B_t43d)
 *
 * Revision 1.1  1997/10/11  00:17:27  avi
 * Initial revision
 **/
/*____________________________________________________________________
  Program:	t4imgs_4dfp

  Description:	12 parameter linear transform and average 4dfp stacks.

  Authors:	Avi Snyder

  History:	Jan-01-97.
  Revision:     25-Aug-04. Rimmon Sachs
                Added nearest_neighbor interpolation option.	
____________________________________________________________________*/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <unistd.h>		/* R_OK */
#include <rec.h>
#include <endianio.h>
#include <Getifh.h>

#define FLT_MAX		1e37
#define MAXL		256

/*************/
/* externals */
/*************/
extern void	f_init (void), f_exit (void);	/* FORTRAN i/o */
extern void	t4_init_ (float *t4);					/* t4_sub.f */
extern void	t4_read_ (char *t4file, float *t4);			/* t4_sub.f */
extern void	to_711_2b_ (float *t4);					/* to_711-2B.f */
extern void	pt4ixyz_ (int *imgdim, float *voxdim, float *centert, float *mmppixt, float *t4mat,
						      float *centero, float *mmppixo, float *t4atl, float *t4);		/* ft4ixyz.f */
extern void	ft4ixyz_ (int *mode, float *t4, float *imgt, int *nxt, int *nyt, int *nzt,
						float *imgo, int *nxo, int *nyo, int *nzo);				/* ft4ixyz.f */
extern void	ft4imgn_ (float *t4, float *imgt, int *nxt, int *nyt, int *nzt, float *centert, float *mmppixt,
				     float *imgo, int *nxo, int *nyo, int *nzo, float *centero, float *mmppixo);	/* ft4imgn.f */
extern void	ft4imgo_ (float *t4, float *imgt, int *nxt, int *nyt, int *nzt, float *centert, float *mmppixt,
				     float *imgo, int *nxo, int *nyo, int *nzo, float *centero, float *mmppixo);	/* ft4imgo.f */
extern void	vrtflip_ (int *iori, int *imgdim, float *centeri, float *mmppixi, float *centert, float *mmppixt);	/* ft4imgo.f */
extern void	splineza_ (float *imgt, int *nx, int *ny, int *nz, float *d2zi);	/* spline3dvgh.f */
extern void	splinezf_ (float *imgt, int *nx, int *ny, int *nz, float *d2zi);	/* spline3dvgh.f */
extern void	splinex_  (float *imgt, int *nx, int *ny, int *nz, float *d2xi);	/* spline3dvgh.f */
extern void	spliney_  (float *imgt, int *nx, int *ny, int *nz, float *d2yi);	/* spline3dvgh.f */
extern float	t4scale (char *t4file);					/* t4scale.c */

extern void	flipx (float *imag, int *nx, int *ny, int *nz);		/* cflip.c */
extern void	flipz (float *imag, int *nx, int *ny, int *nz);		/* cflip.c */
extern int 	x4dfp2ecat (float *imag, int *dim, int orientation);	/* below */
extern void	spline3d (float *imgt, int *dim);			/* below */

float rnan (void) {
	union {
		float		r;
		unsigned int	j;
        } word;
	word.j = 0x7fffffff;
	return word.r;
}

typedef struct {
	char		imgfile[2*MAXL];
	char		t4file [2*MAXL];
	char		mskfile[2*MAXL];
	float		scale;
	float		weight;
	int		lo_slice, hi_slice;
	int		orient;		/* ifh convention: 2 = TRA; 3 = COR; 4 = SAG; */
	int		isbig;
} RUN_INFO;

static char rcsid[] = "$Id: t4imgs_4dfp.c,v 1.36 2011/04/09 00:19:16 avi Exp $";
int main (int argc, char *argv[]) {
	FILE			*lstfp;			/* input image list */
	FILE			*imgfp, *mskfp;		/* input image and mask file pointers */
	FILE			*imgsfp, *imgwfp;	/* output file pointers */
	char			lstfile[MAXL], dimroot[MAXL] = "";
	char			imgroot[2*MAXL];
	char			outroot[MAXL], outfile[MAXL], tmpfile[MAXL];

/***********/
/* utility */
/***********/
	char			*ptr, string[4*MAXL], command[MAXL], program[MAXL];
	int			c, i, j, k, m;

/**************/
/* stack list */
/**************/
	char			*srgv[256];		/* input string field pointers */
	RUN_INFO		*stackspc;
	int			nstack;

/***************/
/* input image */
/***************/
	IFH			ifh;
	float			*imgt, *imgm = NULL;	/* initialization only to keep gcc -Wall quiet */
	float			voxdim[3], voxdimm[3];
	int			imgdim[4], imgdimm[4], outdim[4];
	int			xdim, ydim, zdim, vdim, idim, isbig, isbigm, orientm;

/***************/
/* computation */
/***************/
	float			t4atl[16], t4mat[16], t4[16];	/* affine warps */
	float			q, scale;
	float			mmppixt[3], centert[3];

/****************/
/* output image */
/****************/
	float			*imgo, *imgw, *imgs;
	float			mmppixo[3], centero[3];
	float			maxo, mino;
	int			nxo, nyo, nzo, odim, nframe;
	int			orio = 2, yshift = 0;
	char			control = '\0';

/*********/
/* flags */
/*********/
	int			mode = 2048;		/* ENDSLICE enabled */
	int			status;
	int			B_flag = 0;		/* comwrp=to_711-2B_t4 */
	int			odim_flag = 3;		/* output image dimension control switch */
 	int			spline_flag = 0;	/* use spline3dvgh to interpolate instead of imgvalx */
 	int			sqrt_flag = 0;		/* divide summmed image by sqrt (n) as for z-img averaging */
 	int			NaN_flag = 0;		/* enable NaN output if value undefined */
	int 			nearest_neighbor = 0;   /* interpolate */
	int			saverec = 1;

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);
	t4_init_ (t4mat);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while ((c = *ptr++)) switch (c) {
				case 'R': saverec = 0;		break;
				case 'B': B_flag++;		break;
				case 'N': NaN_flag++;		break;
				case 'n': nearest_neighbor++;   break;
				case 's': spline_flag++;	break;
				case 'z': sqrt_flag++;		break;
				case '@': control = *ptr++;	*ptr = '\0'; break;
				case 'O': if (!strcmp (ptr, "333")) {
						odim_flag = 3;
						if (ptr[3] == '.') yshift = atoi (ptr + 4);
					} else if (!strcmp (ptr, "222")) {
						odim_flag = 2;
					} else if (!strcmp (ptr, "111")) {
						odim_flag = 1;
					} else {
						getroot (ptr, dimroot);
						odim_flag = 0;
					}
								*ptr = '\0'; break;
				default:			break;
			}
		} else switch (k) {
			case 0:	strcpy  (lstfile, argv[i]);	k++; break;
			case 1:	getroot (argv[i], outroot);	k++; break;
			default:				break;
		}	
	}
	if (k < 2) {
		printf ("Usage:\t%s [options] <inlist> <outfile>\n", program);
		printf ("\toption\n");
		printf ("\t-z\tnormalize by sqrt(n) rather than n (for z images)\n");
		printf ("\t-s\tinterpolate by 3D cubic spline (default is 3D linear)\n");
		printf ("\t-N\toutput NaN (default 0.0) for undefined values\n");
		printf ("\t-B\tinternally convert to_711-2A_t4->to_711-2B_t4\n");
		printf ("\t-n\tuse nearest neighbor interpolation\n");
		printf ("\t-R\tsuppress creation of rec file\n");
		printf ("\t-O111\toutput in 111 space instead of default 333.0 space\n");
		printf ("\t-O222\toutput in 222 space instead of default 333.0 space\n");
		printf ("\t-O333.n\toutput in 333.n space (y shifted up by n pixels)\n");
		printf ("\t-Omy_image\tduplicate dimensions of my_image.4dfp.ifh\n");
		printf ("\t-@<b|l>\toutput big or little endian (default CPU endian)\n");
		exit (-1);
	}

/*******************/
/* scan input list */
/*******************/
	if (!(lstfp = fopen (lstfile, "r"))) errr (program, lstfile);
	nstack = 0; while (fgets (string, 4*MAXL, lstfp)) nstack++; rewind (lstfp);
	if (!(stackspc = (RUN_INFO *) malloc (nstack*sizeof (RUN_INFO)))) errm (program);
	nstack = nframe = 0;
	while (fgets (string, 4*MAXL, lstfp)) {
		if ((ptr = strchr (string, '#'))) *ptr = '\0';
		i = m = 0;
		while (m < 256) {
			while (!isgraph ((int) string[i]) && string[i]) i++;
			if (!string[i]) break;
			srgv[m] = string + i; m++;
			while (isgraph ((int) string[i])) i++;
			if (!string[i]) break;
			string[i] = '\0'; i++;
		}
		if (m < 1) continue;			/* blank line */

		strcpy (stackspc[nstack].t4file, "");
		strcpy (stackspc[nstack].mskfile, "");
		stackspc[nstack].scale = 1.0;
		stackspc[nstack].weight = 1.0;
		stackspc[nstack].lo_slice = stackspc[nstack].hi_slice = 0;

		for (i = 0; i < m; i++) {
			if (i == 0) getroot (srgv[i], imgroot);
			if (strstr (srgv[i], "img="))
				getroot (srgv[i] + 4, imgroot);
			if (strstr (srgv[i], "t4="))
				strncpy (stackspc[nstack].t4file,  srgv[i] + 3, 2*MAXL - 1);
			if (strstr (srgv[i], "mask="))
				strncpy (stackspc[nstack].mskfile, srgv[i] + 5, 2*MAXL - 1);
			if (strstr (srgv[i], "scale="))
				stackspc[nstack].scale = atof (srgv[i] + 6);
			if (strstr (srgv[i], "weight="))
				stackspc[nstack].weight = atof (srgv[i] + 7);
			if ((ptr = strstr (srgv[i], "slices="))) {
				if ((ptr = strstr (srgv[i], "to"))) {
					stackspc[nstack].hi_slice = atoi (ptr + 2);
					*ptr = '\0';
				}
				stackspc[nstack].lo_slice = atof (srgv[i] + 7);
				if (!stackspc[nstack].hi_slice) stackspc[nstack].hi_slice = stackspc[nstack].lo_slice;
			}
		}

		sprintf (stackspc[nstack].imgfile, "%s.4dfp.img", imgroot);
		if (get_4dfp_dimoe_quiet (imgroot, imgdim, voxdim, &k, &isbig)) errr (program, stackspc[nstack].imgfile);
		if (k < 2 || k > 4) {
			fprintf (stderr, "%s: invalid %s ifh orientation (=%d)\n", program, imgroot, k);
			exit (-1);
		}

		if (strlen (stackspc[nstack].t4file)) {
			if (access (stackspc[nstack].t4file, R_OK)) errr (program, stackspc[nstack].t4file);
		}

		if (strlen (stackspc[nstack].mskfile)) {		
			getroot (stackspc[nstack].mskfile, stackspc[nstack].mskfile);
			strcat (stackspc[nstack].mskfile, ".4dfp.img");
			if (access (stackspc[nstack].mskfile, R_OK)) errr (program, stackspc[nstack].mskfile);
		}

		nframe = (imgdim[3] > nframe) ? imgdim[3]: nframe;
		nstack++;
	}
	fclose (lstfp);

	switch (odim_flag) {
	case 0:	if (Getifh (dimroot, &ifh)) exit (-1);
		nxo = ifh.matrix_size[0];
		nyo = ifh.matrix_size[1];
		nzo = ifh.matrix_size[2];
		for (k = 0; k < 3; k++) {
			mmppixo[k] = ifh.mmppix[k];
			centero[k] = ifh.center[k];
		}
		orio = ifh.orientation;
		break;
	case 1:	nxo = 176;
		nyo = 208;
		nzo = 176;
		mmppixo[0] =  1.0;
		mmppixo[1] = -1.0;
		mmppixo[2] = -1.0;
		centero[0] =   89.;
		centero[1] =  -85.;
		centero[2] = -101.;
		break;
	case 2:	nxo = 128;
		nyo = 128;
		nzo = 75;
		mmppixo[0] =  2.0;
		mmppixo[1] = -2.0;
		mmppixo[2] = -2.0;
                centero[0] =  129.;
		centero[1] = -129.;
		centero[2] = mmppixo[2] * 41;
		break;
	case 3:	nxo = 48;
		nyo = 64;
		nzo = 48;
		mmppixo[0] =  3.0;
		mmppixo[1] = -3.0;
		mmppixo[2] = -3.0;
		centero[0] =  73.5;					/* 24.5 * 3.0 */
		centero[1] =  (29 - yshift) * mmppixo[2];		/*  -29 * 3.0 */
		centero[2] = -84.0;					/*  -28 * 3.0 */
		break;
	}
	odim = nxo * nyo * nzo;			/* output frame voxel count */
	outdim[0] = nxo;
	outdim[1] = nyo;
	outdim[2] = nzo;
	outdim[3] = nframe;
	if (!control) control = (CPU_is_bigendian()) ? 'b' : 'l';

/******************************/
/* construct output file name */
/******************************/
	sprintf  (outfile, "%s.4dfp.img", outroot);
	sprintf  (tmpfile, "%s.4dfp.imgw", outroot);
	if (saverec) {
		startrece (outfile, argc, argv, rcsid, control);
		if (spline_flag) {
			printrec ("resampling by 3D cubic spline interpolation\n");
		} else if (nearest_neighbor) {
			printrec ("resampling by nearest-neighbor interpolation\n");
		} else {
			printrec ("resampling by 3D linear interpolation\n");
		}
 		printrec ("sub\n"); catrec (lstfile); printrec ("endsub\n");
	}

/*******************************/
/* allocate accumulator memory */
/*******************************/
	imgo = (float *) calloc (odim, sizeof (float));
	imgw = (float *) calloc (odim, sizeof (float));
	imgs = (float *) calloc (odim, sizeof (float));
	if (!imgo || !imgw || !imgs) errm (program);
	printf ("Writing: %s\n", outfile);
	if (!(imgsfp = fopen (outfile, "w+b")) || !(imgwfp = fopen (tmpfile, "w+b"))) errw (program, outfile);
	if (imgdim[3] > 1) for (j = 0; j < imgdim[3]; j++) {
		if (fwrite (imgs, sizeof (float), odim, imgsfp) < odim
		||  fwrite (imgw, sizeof (float), odim, imgwfp) < odim) errw (program, outfile);
	}

	maxo = -FLT_MAX; mino = FLT_MAX;
	for (i = 0; i < nstack; i++) {
		get_4dfp_dimoe_quiet (stackspc[i].imgfile, imgdim, voxdim, &stackspc[i].orient, &stackspc[i].isbig);
		Getifh (stackspc[i].imgfile, &ifh);
		xdim = imgdim[0];
		ydim = imgdim[1];
		zdim = imgdim[2];
		vdim = imgdim[3];
		idim = xdim * ydim * zdim;

		if (!(imgfp = fopen (stackspc[i].imgfile, "r"))) errr (program, stackspc[i].imgfile);
		fprintf (stdout, "Reading: %s\n", stackspc[i].imgfile);

		scale = stackspc[i].scale;
		if (strlen (stackspc[i].t4file)) {
			printf ("Reading: %s\n", stackspc[i].t4file);
			t4_read_ (stackspc[i].t4file, t4atl);
			scale *= t4scale (stackspc[i].t4file);
			if (saverec) {
				printrec ("t4\n");
				catrec (stackspc[i].t4file);
			}
		} else {
			t4_init_ (t4atl);
		}
		if (B_flag) {
			to_711_2b_ (t4atl);
			printf ("modify warp to_711-2A->to_711-2B\n");
		}

		if (!(imgt = (float *) malloc (idim * ((spline_flag) ? 4 : 1) * sizeof (float)))
		||  !(imgm = (float *) malloc (idim * sizeof (float)))) errm (program);
		if (strlen (stackspc[i].mskfile)) {
			if (get_4dfp_dimoe_quiet (stackspc[i].mskfile, imgdimm, voxdimm, &orientm, &isbigm))
					errr (program, stackspc[i].mskfile);
			status = (orientm != stackspc[i].orient);
			for (k = 0; k < 3; k++) {
				status |= (imgdimm[k] != imgdim[k]);
				status |= (fabs (voxdimm[k] - voxdim[k]) > 1.e-5);
			}
			if (status) {
				fprintf (stderr, "%s: %s %s dimension mismatch\n",
					program, stackspc[i].mskfile, stackspc[i].imgfile);
				exit (-1);
			}
			fprintf (stdout, "Reading: %s\n", stackspc[i].mskfile);
			if (!(mskfp = fopen (stackspc[i].mskfile, "r"))
			|| eread (imgm, idim, isbigm, mskfp)
			|| fclose (mskfp)) errr (program, stackspc[i].mskfile);
			if (saverec) catrec (stackspc[i].mskfile);
		} else {
			for (k = 0; k < idim; k++) imgm[k] = 1.0;
		}
		if (stackspc[i].lo_slice) {
			for (k = 0; k < xdim * ydim * (stackspc[i].lo_slice - 1); k++) imgm[k] = 0.;
		}
		if (stackspc[i].hi_slice) {
			for (k = xdim * ydim * stackspc[i].hi_slice; k < idim; k++) imgm[k] = 0.;
		}

		printf ("input orient=%d", stackspc[i].orient);
		printf (" %s\n", ((stackspc[i].isbig) ? "bigendian" : "littleendian"));
		printf ("input stack dimensions %10d%10d%10d%10d\n", imgdim[0], imgdim[1], imgdim[2], imgdim[3]);
		printf ("input stack mmppix     %10.6f%10.6f%10.6f\n", ifh.mmppix[0], ifh.mmppix[1], ifh.mmppix[2]);
		printf ("input stack center     %10.4f%10.4f%10.4f\n", ifh.center[0], ifh.center[1], ifh.center[2]);
		printf ("imgfile=%s\n", stackspc[i].imgfile);
		if (strlen (stackspc[i].mskfile))	printf ("maskfile=%s\n", stackspc[i].mskfile);
		if (strlen (stackspc[i].t4file))	printf ("t4file=%s\n",   stackspc[i].t4file);
		if (stackspc[i].lo_slice)		printf ("slices=%dto%d  ",
								stackspc[i].lo_slice, stackspc[i].hi_slice);
		printf ("scale=%f  weight=%f\n", scale, stackspc[i].weight);

		printf ("frame:");
		for (j = 0; j < imgdim[3]; j++) {
			printf (" %d", j + 1); fflush (stdout);
			if (eread (imgt, idim, stackspc[i].isbig, imgfp)) errr (program, stackspc[i].imgfile);
			if (spline_flag) spline3d (imgt, imgdim);
			for (k = 0; k < idim; k++) if (!(imgm[k] > 0.)) imgt[k] = rnan ();
/*****************************************/
/* virtual flip instead of x4dfp2ecat () */
/*****************************************/
			vrtflip_ (&stackspc[i].orient, imgdim, ifh.center, ifh.mmppix, centert, mmppixt);
/*********************/
/* execute transform */
/*********************/
			if (spline_flag) {
	pt4ixyz_ (imgdim, voxdim, centert, mmppixt, t4mat, centero, mmppixo, t4atl, t4);
	ft4ixyz_ (&mode, t4, imgt, imgdim+0, imgdim+1, imgdim+2, imgo, outdim+0, outdim+1, outdim+2);
			} else if (nearest_neighbor) {
	ft4imgn_ (t4atl, imgt, &xdim, &ydim, &zdim, centert, mmppixt, imgo, &nxo, &nyo, &nzo, centero, mmppixo);
			} else {
	ft4imgo_ (t4atl, imgt, &xdim, &ydim, &zdim, centert, mmppixt, imgo, &nxo, &nyo, &nzo, centero, mmppixo);
			}
/*******************************/
/* restore 4dfp voxel indexing */
/*******************************/
			x4dfp2ecat (imgo, outdim, orio);

/********************************************/
/* bring accumulator frame back into memory */
/********************************************/
			if (imgdim[3] > 1) {
				if (fseek (imgsfp, (long) j*odim*sizeof (float), SEEK_SET)
				||  fseek (imgwfp, (long) j*odim*sizeof (float), SEEK_SET)
				||  fread (imgs, sizeof (float), odim, imgsfp) < odim
				||  fread (imgw, sizeof (float), odim, imgwfp) < odim) errr (program, outfile);
			}
			for (k = 0; k < odim; k++) {
				if (!isnan (imgo[k])) {
					imgs[k] += imgo[k] * scale * stackspc[i].weight;
					imgw[k] += stackspc[i].weight;
				}
			}
/******************************************/
/* store accumulator frames as CPU endian */
/******************************************/
			if (imgdim[3] > 1) {
				if (fseek (imgsfp, (long) j*odim*sizeof (float), SEEK_SET)
				||  fseek (imgwfp, (long) j*odim*sizeof (float), SEEK_SET)
				||  fwrite (imgs, sizeof (float), odim, imgsfp) < odim
				||  fwrite (imgw, sizeof (float), odim, imgwfp) < odim) errw (program, outfile);
			}
		}			/* end frame loop */
		printf ("\n");
		fclose (imgfp);
		free (imgt); free (imgm);
		if (saverec) catrec (stackspc[i].imgfile);
	}				/* end stack loop */

/*******************************/
/* normalize stack accumulator */
/*******************************/
	rewind (imgsfp); rewind (imgwfp);
	maxo = -FLT_MAX; mino = FLT_MAX;
	for (j = 0; j < imgdim[3]; j++) {
		if (imgdim[3] > 1) {
			if (fseek (imgsfp, (long) j*odim*sizeof (float), SEEK_SET)
			||  fseek (imgwfp, (long) j*odim*sizeof (float), SEEK_SET)
			||  fread (imgs, sizeof (float), odim, imgsfp) < odim
			||  fread (imgw, sizeof (float), odim, imgwfp) < odim) errr (program, outfile);
		}
		for (k = 0; k < odim; k++) {
			if (imgw[k] > 0.) {
				q = (sqrt_flag) ? sqrt (imgw[k]) : imgw[k];
				imgs[k] /= q;
			} else {
				imgs[k] = (NaN_flag) ? rnan () : 0.0;
			}
			if (!isnan (imgs[k])) {
				if (imgs[k] < mino) mino = imgs[k];
				if (imgs[k] > maxo) maxo = imgs[k];
			}
		}
		if (fseek (imgsfp, (long) j*odim*sizeof (float), SEEK_SET)
		||  ewrite (imgs, odim, control, imgsfp)) errw (program, outfile);
	}
	if (fclose (imgsfp) || fclose (imgwfp)) errw (program, outfile);
	if (remove (tmpfile)) errw (program, tmpfile);

	voxdim[0] = (mmppixo[0] > 0.) ? mmppixo[0]: -mmppixo[0];
	voxdim[1] = (mmppixo[1] > 0.) ? mmppixo[1]: -mmppixo[1];
	voxdim[2] = (mmppixo[2] > 0.) ? mmppixo[2]: -mmppixo[2];
	printf ("output orient=%d\n", orio);
	printf ("output stack dimensions%10d%10d%10d%10d\n",    outdim[0],  outdim[1],  outdim[2],  outdim[3]);
	printf ("output stack mmppix    %10.6f%10.6f%10.6f\n", mmppixo[0], mmppixo[1], mmppixo[2]);
	printf ("output stack center    %10.4f%10.4f%10.4f\n", centero[0], centero[1], centero[2]);
/***************************/
/* create interfile header */
/***************************/
	writeifhmce (program, outfile, outdim, voxdim, orio, mmppixo, centero, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s -r%dto%d", outroot, (int) (mino - 0.5), (int) (maxo + 0.5));
	printf ("%s\n", command); status = system (command);

/**********/
/* endrec */
/**********/
	if (saverec) endrec ();

/*********************/
/* clean up and exit */
/*********************/
	free (stackspc);
	free (imgo); free (imgw); free (imgs);
	exit (status);
} 

int x4dfp2ecat (float *imag, int *dim, int orientation) {
	switch (orientation) {
		case 2:	flipx (imag, dim+0, dim+1, dim+2);	/* transverse */
			flipz (imag, dim+0, dim+1, dim+2);
			break;
		case 3:	flipx (imag, dim+0, dim+1, dim+2);	/* coronal */
			break;
		case 4: break;					/* sagittal */
		default: return -1;				/* none of the above */
	}
	return 0;
}

void spline3d (float *imgt, int *dim) {
	int		vdim;

	vdim = dim[0] * dim[1] * dim[2];
		splinex_  (imgt, dim+0, dim+1, dim+2, imgt + vdim*1);
		spliney_  (imgt, dim+0, dim+1, dim+2, imgt + vdim*2);
	if (dim[2] < 48) {
		splineza_ (imgt, dim+0, dim+1, dim+2, imgt + vdim*3);
	} else {
		splinezf_ (imgt, dim+0, dim+1, dim+2, imgt + vdim*3);
	}
}
