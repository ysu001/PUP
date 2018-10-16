/*$Header: /data/petsun4/data1/src_solaris/wrpsmg_4dfp/RCS/wrpsmg_4dfp.c,v 1.11 2010/12/23 03:37:05 avi Exp $*/
/*$Log: wrpsmg_4dfp.c,v $
 * Revision 1.11  2010/12/23  03:37:05  avi
 * long -> int in rnan()
 *
 * Revision 1.10  2010/01/15  01:56:36  avi
 * operate on multi-volume input
 *
 * Revision 1.9  2010/01/14  06:26:47  avi
 * Solaris 10 and linux compliant
 *
 * Revision 1.8  2006/03/06  05:28:11  avi
 * replace fflip.f subroutines with cflip.c subroutines
 *
 * Revision 1.7  2004/09/02  20:57:46  rsachs
 * Installed use of fns 'errm','errr','errw' & 'getroot'.
 *
 * Revision 1.6  2002/09/20  03:48:41  avi
 * eliminate wrapped lines in hard copy
 *
 * Revision 1.5  2002/09/20  03:26:10  avi
 * implement dimension consistency tests
 * eliminate call to get4dfpDimN and need for libmri
 *
 * Revision 1.4  2002/09/19  01:21:53  avi
 * xtreamlined auxio
 *
 * Revision 1.3  2002/09/19  00:29:20  avi
 * reoganize jj loops to eliminate non jj dependent code
 *
 * Revision 1.2  2002/09/18  06:39:18  avi
 * correct intensity scaling (separate scale factors for each member of pair)
 **/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <unistd.h>		/* R_OK */
#include <rec.h>
#include <endianio.h>
#include <Getifh.h>
#include <t4_sub.h>

#define MAXL	256
#define NF	2

typedef struct {
	char	imgfile[NF][MAXL];
	char	t4file[NF][MAXL];
	char	mskfile[MAXL];
	float	weight;
	int	lo_slice, hi_slice;
	int	orient;			/* ifh convention: 2 = TRA; 3 = COR; 4 = SAG; */
	int	isbig[NF], isbigm;
	int	imgdim[4];
	float	mmppix[3];
	float	center[3];
} RUN_INFO;

float rnan (void) {
	union {
		float		r;
		unsigned int	j;
        } word;
	word.j = 0x7fffffff;
	return word.r;
}

void errs (char *program) {
	fprintf (stderr, "%s: invalid input list file\n", program);
	exit (-1);
}

extern void     vrtflip_ (), ft4imgo_ ();				/* FORTRAN */
extern void	flipx (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */
extern void	flipz (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */
extern float	t4scale (char *t4file);					/* t4scale.c */

/*********/
/* below */
/*********/
extern int 	x4dfp2ecat (float *imag, int *dim, int orientation);
extern void     imgsord (float t4[NF][16], float scale[NF], float **imgt, int *nxt, float *centert, float *mmppixt,
			float *imgo, int *nxo, float *centero, float *mmppixo, float *imgd, int n1or2);  
extern void 	auxio (char *outroot, char *aux, char *title, char **argv, int argc, int odim, float *img);
extern int	dimcheck (IFH *ifh, RUN_INFO *stackspc);

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: wrpsmg_4dfp.c,v 1.11 2010/12/23 03:37:05 avi Exp $";
static char	program[MAXL];
static int	status = 0;
static char	control = '\0';

int main (int argc, char *argv[]) {
	FILE		*lstfp;			/* input image list */
	FILE		*imgfp, *mskfp;		/* input image & mask file ptrs */
	char		lstfile[MAXL], dimfile[MAXL] = "";
	char		imgroot[NF][MAXL], imgfile[MAXL], title[MAXL]; 
	char		outroot[MAXL], outfile[MAXL];

/***********/
/* utility */
/***********/
	char		*ptr, string[2*MAXL], command[MAXL];
	int		c, i, j, k, l, m, jj;

/**************/
/* stack list */
/**************/
	char		*srgv[256];		/* input string field ptrs */
	RUN_INFO	*stackspc;
	int		nstack;

/***************/
/* input image */
/***************/
	IFH		ifh;
	float		*imgt[NF], *imgm;
	float		mmppixt[3], centert[3];
	int		n1or2, n1or2g;
	int		slcdim, vdim, tdim;

/***************/
/* computation */
/***************/
	float		t4atl[NF][16];		/* affine warps */
	float		d, q, scale[NF], wgt;

/****************/
/* output image */
/****************/
	float		*imgo, *imgw, *imgs, *imgu, *imgv, *imgd;
	float		mmppixo[3], centero[3], voxdim[3];
	int		outdim[4], nxo, nyo, nzo, odim, pdim, nframe;
	int		orio = 2, yshift = 0;

/*********/
/* flags */
/*********/
	int		debug = 0, keep_wei = 0, keep_sd1 = 0;
	int		odim_flag = 2;		/* output image dimension control switch */
 	int		NaN_flag = 0;		/* enable NaN output if val= undef */

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'N': NaN_flag++;		break;
				case 'w': keep_wei++;		break;
				case 's': keep_sd1++;		break;
				case '@': control = *ptr++;	*ptr = '\0'; break;
				case 'O':
					if (!strncmp (ptr, "333", 3)) {
						odim_flag = 3;
						if (ptr[3] == '.') yshift = atoi(ptr + 4);
					} else if (!strcmp (ptr, "222")) {
						odim_flag = 2;
					} else if (!strcmp (ptr, "111")) {
						odim_flag = 1;
					} else {
						getroot (ptr, dimfile);
						odim_flag = 0;
					}
					*ptr = '\0'; break;
			}	/* end switch */
		} else switch (k) {
			case 0:	strcpy (lstfile, argv[i]);	k++; break;
			case 1:	getroot (argv[i], outroot);	k++; break;
		}	
	} 
	if (k < 2) {
		printf ("Usage:\t%s [options] <inlist> <outfile>\n", program);
		printf ("\toption\n");
		printf ("\t-N\toutput NaN (default 0.0) for undefined values\n");
		printf ("\t-w\tcreate sum of weights image\n");
		printf ("\t-s\tcreate square root variance (sd) image\n");
		printf ("\t-O111\toutput in 111 space\n");
		printf ("\t-O222\toutput in 222 space (default)\n");
		printf ("\t-O333.n\toutput in 333.n space (y shifted up by n pixels)\n");
		printf ("\t-Omy_image\tduplicate dimensions of my_image.4dfp.ifh\n");
		printf ("\t-@<b|l>\toutput big or little endian (default CPU endian)\n");
		exit (1);
	}

/*******************/
/* scan input list */
/*******************/
	if (!(lstfp = fopen (lstfile, "r"))) errr (program, lstfile);
	nstack = 0; while (fgets (string, 2*MAXL, lstfp)) nstack++; 
        rewind (lstfp);
	if (!(stackspc = (RUN_INFO *) malloc (nstack * sizeof (RUN_INFO)))) errm (program);
	nstack = 0;
	while (fgets (string, 2*MAXL, lstfp)) {
		fprintf (stderr, "%s", string);
		if (ptr = strchr (string, '#')) *ptr = '\0';
		i = m = 0;
		while (m < 256) {
			while (!isgraph ((int) string[i]) && string[i]) i++;
			if (!string[i]) break;
			srgv[m] = string + i; m++;
			while (isgraph ((int) string[i])) i++;
			if (!string[i]) break;
			string[i] = '\0'; i++;
		}
		if (m < 1) continue;		/* blank line */

                for (jj = 0; jj < NF; jj++) strcpy (stackspc[nstack].t4file[jj], "");
		strcpy (stackspc[nstack].mskfile, "");
		stackspc[nstack].weight = 1.0;
		stackspc[nstack].lo_slice = stackspc[nstack].hi_slice = 0;

		n1or2 = 0;
		for (i = 0; i < m; i++) {
			if (strstr (srgv[i], "t4=")) {
				if (!n1or2) errs (program);
				strcpy (stackspc[nstack].t4file[n1or2 - 1], srgv[i] + 3);
			}
			else if (strstr (srgv[i], "mask="))
				getroot (srgv[i] + 5, stackspc[nstack].mskfile);
			else if (strstr (srgv[i], "weight="))
				stackspc[nstack].weight = atof (srgv[i] + 7);
			else if (ptr = strstr (srgv[i], "slices=")) {
					if (ptr = strstr (srgv[i], "to")) {
						stackspc[nstack].hi_slice = atoi (ptr + 2);
						*ptr = '\0';
					}
					stackspc[nstack].lo_slice = atof (srgv[i] + 7);
			 		if (!stackspc[nstack].hi_slice) errs (program);
			} else {
				getroot (srgv[i], imgroot[n1or2++]);
			}
			if (n1or2 > 2) errs (program);
		}

		for (jj = 0; jj < n1or2; jj++) {
			sprintf (stackspc[nstack].imgfile[jj], "%s.4dfp.img", imgroot[jj]);
			if (Getifh (stackspc[nstack].imgfile[jj], &ifh)) errr (program, imgroot[jj]);
			stackspc[nstack].isbig[jj] = strcmp (ifh.imagedata_byte_order, "littleendian");

			switch (jj) {
			case 0:	stackspc[nstack].orient = ifh.orientation;
				for (k = 0; k < 3; k++) {
					stackspc[nstack].imgdim[k] = ifh.matrix_size[k];
					stackspc[nstack].mmppix[k] = ifh.mmppix[k];
					stackspc[nstack].center[k] = ifh.center[k];
				}
				stackspc[nstack].imgdim[3] = ifh.matrix_size[3];
				break;
			case 1:	if (dimcheck (&ifh, stackspc + nstack)
				||  stackspc[nstack].imgdim[3] != ifh.matrix_size[3]) {
					fprintf (stderr, "%s: %s %s dimension mismatch\n",
						program, imgroot[0], imgroot[1]);
					exit (-1);
				}
				break;
			}
			if (strlen (stackspc[nstack].t4file[jj])) {
				if (access (stackspc[nstack].t4file[jj], R_OK))
					errr (program, stackspc[nstack].t4file[jj]);
			}
		}

/**********************************************/
/* check nframe consistency across input list */
/**********************************************/
		if (nstack) {
			if (stackspc[nstack].imgdim[3] != nframe) {
				fprintf (stderr, "%s: %s %s dimension mismatch\n", program,
					stackspc[nstack].imgfile[0], stackspc[0].imgfile[0]);
				exit (-1);
			}
		} else {
			nframe = stackspc[nstack].imgdim[3];
		}

		if (stackspc[nstack].hi_slice > stackspc[nstack].imgdim[2]
		||  stackspc[nstack].hi_slice < stackspc[nstack].lo_slice
		||  stackspc[nstack].lo_slice < 0) errs (program);

		if (strlen (stackspc[nstack].mskfile)) {
			if (Getifh (stackspc[nstack].mskfile, &ifh)) errr (program, stackspc[nstack].mskfile);
			stackspc[nstack].isbigm = strcmp (ifh.imagedata_byte_order, "littleendian");
			if (dimcheck (&ifh, stackspc + nstack)) {
				fprintf (stderr, "%s: %s %s dimension mismatch\n",
					program, imgroot[0], stackspc[nstack].mskfile);
				exit (-1);
			}
			strcat (stackspc[nstack].mskfile, ".4dfp.img");
			if (access (stackspc[nstack].mskfile, R_OK)) errr (program, stackspc[nstack].mskfile);
		}
		if (nstack) {
			if (n1or2 != n1or2g) errs (program);
		} else {
			n1or2g = n1or2;
		}	
		nstack++;
	}			/* end fgets loop */
	fclose (lstfp);

	if (debug) for (i = 0; i < nstack; i++) {
		printf ("%-32s\tt4=%-32s\t%-32s\tt4=%-32s\n",
		stackspc[i].imgfile[0], stackspc[i].t4file[0],
		stackspc[i].imgfile[1], stackspc[i].t4file[1]);
	}

	switch (odim_flag) {
	case 0:	strcat (dimfile, ".4dfp");
		if (Getifh (dimfile, &ifh)) errr (program, dimfile);
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
		centero[0] =  73.5;				/* 24.5 * 3.0 */
		centero[1] = (29 - yshift) * mmppixo[2];	/*  -29 * 3.0 */
		centero[2] = -84.0;				/*  -28 * 3.0 */
		break;
	}
	pdim = nxo * nyo * nzo;
	odim = pdim * nframe;		/* output voxel count */
	outdim[0] = nxo;
	outdim[1] = nyo;
	outdim[2] = nzo;
	outdim[3] = nframe;
	if (!control) control = (CPU_is_bigendian()) ? 'b' : 'l';

/******************************/
/* construct output file name */
/******************************/
	sprintf   (outfile, "%s.4dfp.img", outroot);
	startrece (outfile, argc, argv, rcsid, control);
	printrec  ("resampling by 3D linear interpolation\n");
 	printrec  ("sub\n"); catrec (lstfile); printrec ("endsub\n");

/*******************************/
/* allocate accumulator memory */
/*******************************/
	imgo = (float *) calloc (odim, sizeof (float));
        imgd = (float *) calloc (odim, sizeof (float));
	imgw = (float *) calloc (odim, sizeof (float));
	imgs = (float *) calloc (odim, sizeof (float));
	imgv = (float *) calloc (odim, sizeof (float));
	imgu = (float *) calloc (odim, sizeof (float));
	if (!imgo || !imgw || !imgs || !imgd || !imgv || !imgu) errm (program);

	for (i = 0; i < nstack; i++) {         /* begin stack loop */
		slcdim = stackspc[i].imgdim[0] * stackspc[i].imgdim[1];
		vdim = slcdim * stackspc[i].imgdim[2];
		if (!(imgm = (float *) malloc (vdim * sizeof (float)))) errm (program);
		if (strlen (stackspc[i].mskfile)) {
			fprintf (stdout, "Reading: %s\n", stackspc[i].mskfile);
			if (!(mskfp = fopen (stackspc[i].mskfile, "r"))
			|| eread (imgm, vdim, stackspc[i].isbigm, mskfp)
			|| fclose (mskfp)) errr (program, stackspc[i].mskfile);
		} else {
			for (k = 0; k < vdim; k++) imgm[k] = 1.0;
		}
		if (stackspc[i].lo_slice) {
			for (k = 0; k < slcdim * (stackspc[i].lo_slice - 1); k++) imgm[k] = 0.;
		}
		if (stackspc[i].hi_slice) {
			for (k = slcdim * stackspc[i].hi_slice; k < vdim; k++) imgm[k] = 0.;
		}
		printf ("input orient=%d\n", stackspc[i].orient);
		printf ("input stack dimensions %10d%10d%10d%10d\n",
			stackspc[i].imgdim[0], stackspc[i].imgdim[1], stackspc[i].imgdim[2], stackspc[i].imgdim[3]);
		printf ("input stack mmppix     %10.6f%10.6f%10.6f\n",
			stackspc[i].mmppix[0], stackspc[i].mmppix[1], stackspc[i].mmppix[2]);
		printf ("input stack center     %10.4f%10.4f%10.4f\n",
			stackspc[i].center[0], stackspc[i].center[1], stackspc[i].center[2]);
		if (strlen (stackspc[i].mskfile))	
			printf ("maskfile=%s\n", stackspc[i].mskfile);
		if (stackspc[i].lo_slice)
			printf ("slices=%dto%d\n", stackspc[i].lo_slice, stackspc[i].hi_slice);

		for (jj = 0; jj < n1or2; jj++) {
			if (strlen (stackspc[i].t4file[jj])) {
				printf ("Reading: %s\n", stackspc[i].t4file[jj]);
				t4_read_ (stackspc[i].t4file[jj], t4atl[jj]);
				scale[jj] = t4scale (stackspc[i].t4file[jj]);
			} else {
				t4_init_ (t4atl[jj]);
				scale[jj] = 1.0;
			}
			printf ("scale=%f\n", scale[jj]);
			tdim = nframe * vdim;
			if (!(imgt[jj] = (float *) malloc (tdim * sizeof (float)))) errm (program);
			fprintf (stdout, "Reading: %s\n", stackspc[i].imgfile[jj]);
			if (!(imgfp = fopen (stackspc[i].imgfile[jj], "r"))
			||  eread (imgt[jj], tdim, stackspc[i].isbig[jj], imgfp)
			||  fclose (imgfp)) errr (program, stackspc[i].imgfile[jj]);
			for (m = 0; m < nframe; m++) {
				for (k = 0; k < vdim; k++) if (!(imgm[k] > 0.)) imgt[jj][m*vdim + k] = rnan ();
			}
/*****************************************/
/* virtual flip instead of x4dfp2ecat () */
/*****************************************/
			vrtflip_ (&stackspc[i].orient, stackspc[i].imgdim, stackspc[i].center, stackspc[i].mmppix,
							centert, mmppixt);
			fflush (stdout);
		}                        /* end jj loop */
		free (imgm);

/*********************/
/* execute transform */
/*********************/
		imgsord (t4atl, scale, imgt, stackspc[i].imgdim, centert, mmppixt,
							imgo, outdim, centero, mmppixo, imgd, n1or2);
		for (jj = 0; jj < n1or2; jj++) { 
			free (imgt[jj]);
			if (strlen (stackspc[i].t4file[jj])) {
				printrec ("t4\n");  
				catrec (stackspc[i].t4file[jj]);
			}
			catrec (stackspc[i].imgfile[jj]);
		}
		if (strlen (stackspc[i].mskfile)) catrec (stackspc[i].mskfile);

/*******************************/
/* restore 4dfp voxel indexing */
/*******************************/
		for (m = 0; m < nframe; m++) x4dfp2ecat (imgd + m*pdim, outdim, orio);

/***************************************/
/* compute variance-associated volumes */
/***************************************/
		printf ("weight=%f\n", stackspc[i].weight);
                wgt = stackspc[i].weight;
		for (k = 0; k < odim; k++) {
			if (!isnan (imgd[k])) {
				q = imgd[k]; 
				imgs[k] += wgt * q;
				imgw[k] += wgt;
				imgv[k] += wgt * q * q;
				imgu[k] += wgt * wgt;
	        	}
		}
	}   				/* end stack loop */

/*******************************/
/* normalize stack accumulator */
/*******************************/
	printf ("computing statistical images\n");
	for (k = 0; k < odim; k++) {
		if (imgw[k] > 0.) {
	     		imgs[k] /= imgw[k];
			d = imgw[k]*imgw[k] - imgu[k];
			if (d > 0.) {
				q = (imgw[k]*imgw[k]) / d;
				imgv[k] = sqrt (q*(imgv[k]/imgw[k] - imgs[k]*imgs[k]));
			} else {
				imgv[k] = 0.;
			}
	 	} else {
	     		imgd[k] = (NaN_flag) ? rnan () : 0.0;
	     		imgs[k] = (NaN_flag) ? rnan () : 0.0;
	 	}
    	}

/**********************/
/* write output stack */
/**********************/
	voxdim[0] = (mmppixo[0] > 0.) ? mmppixo[0]: -mmppixo[0];
	voxdim[1] = (mmppixo[1] > 0.) ? mmppixo[1]: -mmppixo[1];
	voxdim[2] = (mmppixo[2] > 0.) ? mmppixo[2]: -mmppixo[2];
	printf ("output orient=%d\n", orio);
	printf ("output stack dimensions%10d%10d%10d%10d\n", outdim[0], outdim[1], outdim[2], outdim[3]);
	printf ("output stack mmppix    %10.6f%10.6f%10.6f\n", mmppixo[0], mmppixo[1], mmppixo[2]);
	printf ("output stack center    %10.4f%10.4f%10.4f\n", centero[0], centero[1], centero[2]);
	printf ("Writing: %s\n", outfile);   
	if (!(imgfp = fopen (outfile, "wb"))
	|| ewrite (imgs, odim, control, imgfp)
	|| fclose (imgfp)) errw (program, outfile);

/***************/
/* ifh and hdr */
/***************/
	writeifhmce  (program, outfile, outdim, voxdim, orio, mmppixo, centero, control);
	sprintf (command, "ifh2hdr %s", outroot);
	printf ("%s\n", command); status = system (command);
	endrec ();

	if (keep_wei) auxio (outroot, "wei", "sum of weights", argv, argc, odim, imgw); 
	if (keep_sd1) auxio (outroot, "sd1", "variance image", argv, argc, odim, imgv); 

	free (stackspc);
	free (imgd); free (imgo); free (imgw); free (imgs); free (imgv); free (imgu);
	exit (status);
} 

int x4dfp2ecat (float *imag, int *dim, int orientation) {
	switch (orientation) {
		case 2:		flipz (imag, dim+0, dim+1, dim+2);	/* transverse */
		case 3:		flipx (imag, dim+0, dim+1, dim+2);	/* coronal */
		case 4:		break;					/* sagittal */
		default:	return -1;				/* none of the above */
	}
	return 0;
}

/****************************************************************
void imgsord - function to do image addition/subtraction
Input Arguments:
	float **imgt = *imgt[2] - the 2 input images

Output Arguments:
	imgd = The sum or difference image
	imgt[0]          , if n1or2 == 1
	imgt[0] - imgt[1], if n1or2 == 2
*****************************************************************/
void imgsord (float t4[NF][16], float scale[NF], float **imgt, int *nxt, float *centert,
		float *mmppixt, float *imgo, int *nxo, float *centero, 
		float *mmppixo, float *imgd, int n1or2) {
	int k, m, vdim, pdim;

	vdim = nxt[0] * nxt[1] * nxt[2];
	pdim = nxo[0] * nxo[1] * nxo[2];

	for (m = 0; m < nxt[3]; m++) {
		ft4imgo_ (t4[0], imgt[0] + m*vdim, nxt+0, nxt+1, nxt+2, centert, mmppixt,
				    imgd + m*pdim, nxo+0, nxo+1, nxo+2, centero, mmppixo);
	}
	for (k = 0; k < pdim * nxt[3]; k++) imgd[k] *= scale[0];
	if (n1or2 == 1) return;
	for (m = 0; m < nxt[3]; m++) {
		ft4imgo_ (t4[1], imgt[1] + m*vdim, nxt+0, nxt+1, nxt+2, centert, mmppixt,
				    imgo + m*pdim, nxo+0, nxo+1, nxo+2, centero, mmppixo);
	}
	for (k = 0; k < pdim * nxt[3]; k++) imgd[k] -= imgo[k]*scale[1];
}

void auxio (char *outroot, char *aux, char *title, char **argv, int argc, int odim, float *imgx) {
	FILE	*imgfp;
	char	command[MAXL], auxfile[MAXL], auxroot[MAXL], outfile[MAXL]; 

	sprintf (outfile, "%s.4dfp.img", outroot);
	sprintf (auxroot, "%s_%s", outroot, aux); 
	sprintf (auxfile, "%s.4dfp.img", auxroot);

	printf ("Writing: %s\n", auxfile);
	if (!(imgfp = fopen (auxfile, "wb"))
	|| ewrite (imgx, odim, control, imgfp)
	|| fclose (imgfp)) errw (program, auxfile);

	startrece (auxfile, argc, argv, rcsid, control);
	sprintf (command, "%s %s\n", outroot, title); printrec (command); 
	catrec (outfile);
	endrec ();
	sprintf (command, "/bin/cp %s.4dfp.ifh %s.4dfp.ifh", outroot, auxroot);
	status |= system (command);
	sprintf (command, "ifh2hdr %s", auxroot);
	printf ("%s\n", command); status = system (command);
}

int dimcheck (IFH *ifh, RUN_INFO *stackspc) {
	int	k, m;

	m = stackspc->orient - ifh->orientation;
	for (k = 0; k < 3; k++) {
		m |= 	   (stackspc->imgdim[k] - ifh->matrix_size[k]);
		m |= (fabs (stackspc->mmppix[k] - ifh->mmppix[k]) > 0.0001);
		m |= (fabs (stackspc->center[k] - ifh->center[k]) > 0.0001);
	}
	return m;
}
