/*$Header: /data/petsun4/data1/src_solaris/t4_actmapf_4dfp/RCS/t4_xr3d_4dfp.c,v 1.12 2007/09/20 03:16:32 avi Exp $*/
/*$Log: t4_xr3d_4dfp.c,v $
 * Revision 1.12  2007/09/20  03:16:32  avi
 * linux compliant
 *
 * Revision 1.11  2006/09/29  04:56:02  avi
 * Solaris 10
 *
 * Revision 1.10  2006/09/29  03:00:40  avi
 * before making endian invariant
 *
 * Revision 1.9  2006/09/29  02:46:30  avi
 * rsachs modernization
 *
 * Revision 1.8  2001/05/30  19:52:12  avi
 * -E (undefined voxels = 1.e-37) option
 *
 * Revision 1.7  2001/03/13  00:29:53  avi
 * -O111 now means 176 x 208 x 176 (111h convention)
 *
 * Revision 1.6  1999/03/11  02:39:37  avi
 * eliminate separate spline interpolation constant image pointers
 * better usage
 *
 * Revision 1.5  1999/02/08  23:44:04  avi
 * correct bug in -O333 parsing
 *
 * Revision 1.4  1999/02/04  02:30:44  avi
 * filter out dimfile terminal ".ifh"
 *
 * Revision 1.3  1999/02/04  02:16:26  avi
 * connect -O options to working code
 *
 * Revision 1.2  1999/01/13  00:54:15  avi
 * correct usage
 * fix t4_file R_OK error message
 *
 * Revision 1.1  1999/01/12  06:38:12  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <rec.h>
#include <endianio.h>
#include <Getifh.h>
#include <spline23dvgh.h>

#define MAXL		256
#define YSHIFT		0

/*************/
/* externals */
/*************/
extern void	pt4ixyz_ (int *imgdim, float *voxdim, float *centert, float *mmppixt, float *t4mat,
						      float *centero, float *mmppixo, float *t4atl, float *t4);	/* ft4ixyz.f */
extern void	ft4ixyz_ (int *mode, float *t4, float *imgt, int *nxt, int *nyt, int *nzt,
						float *imgo, int *nxo, int *nyo, int *nzo);			/* ft4ixyz.f */
extern void	t4_init_ (float *t4);					/* t4_sub.f */
extern void	t4_list_ (float *t4);					/* t4_sub.f */
extern void	t4_read_ (char *t4file, float *t4);			/* t4_sub.f */
extern void	ft4imgo_ (float *t4, float *imgt, int *nxt, int *nyt, int *nzt, float *centert, float *mmppixt,
				     float *imgo, int *nxo, int *nyo, int *nzo, float *centero, float *mmppixo);	/* ft4imgo.f */
extern void	vrtflip_ (int *iori, int *imgdim, float *centeri, float *mmppixi, float *centert, float *mmppixt);	/* ft4imgo.f */
extern void	spline3d (float *imgt, int *dim);			/* below */
extern int 	x4dfp2ecat (float *imag, int *dim, int orientation);	/* below */
extern void	flipx (float *imgf, int *pnx, int *pny, int *pnz);	/* cflip.c */
extern void	flipy (float *imgf, int *pnx, int *pny, int *pnz);	/* cflip.c */
extern void	flipz (float *imgf, int *pnx, int *pny, int *pnz);	/* cflip.c */

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[] = "$Id: t4_xr3d_4dfp.c,v 1.12 2007/09/20 03:16:32 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* filenames */
/*************/
	FILE		*outfp, *imgfp, *matfp;
	char		command[MAXL], program[MAXL], *ptr;
	char		imgroot[MAXL];			/* input 4dfp stack filename */
	char		outfile[MAXL], outroot[MAXL];
	char		imgfile[MAXL], matfile[MAXL];
	char		dimfile[MAXL] = "";
	char		t4file [MAXL];
	char		trailer[MAXL] = "xr3d";		/* default output filename trailer */

/***************/
/* imput image */
/***************/
	IFH		ifh;
	float		*imgt;
	float		voxdim[3];
	int		imgdim[4], vdim, isbig;
	char		control = '\0';
	
/*******************/
/* resampled image */
/*******************/
	float		*imgo;
	float		t4mat[16], t4atl[16], t4[16], s4[4], cscale = 1.0;
	float		mmppixt[3], centert[3];
	float		mmppixo[3], centero[3];
	float		fmin = +FLT_MAX;
	float		fmax = -FLT_MAX;
	int		odim, outdim[4], orient;
	int		orio = 2, yshift = YSHIFT;
	int		odim_flag = 3;		/* output image dimension control switch */
	
/***********/
/* utility */
/***********/
	int		c, i, j, k;
	float		q;
	static float	mmppix1[3] = {1., 1., 1.};	/* for using ft4imgo in index->index mode */
	static float	center0[3] = {0., 0., 0.};	/* for using ft4imgo in index->index mode */
	
/*********/
/* flags */
/*********/
	int		mode = 2048;		/* ENDSLICE enabled */
	int		fast_flag = 0;
	int		NaN_flag = 'E';		/* 'E' 1.e-37; 'Z' 0.0; 'N' NaN; */
	int		volnorm = 0;
	int		verbose = 0;
	int		debug = 0;
	int		status = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]);
			ptr = command; while (c = *ptr++) switch (c) {
			case 'e': verbose++;	break;
			case 'f': fast_flag++;	break;
			case 'N': case 'Z': case 'E': NaN_flag = c;	break;
			case 'v': volnorm = atoi (ptr);		*ptr = '\0'; break;
			case 'c': cscale = atof (ptr);		*ptr = '\0'; break;
			case 'a': strcpy (trailer, ptr);	*ptr = '\0'; break;
			case '@': control = *ptr++;		*ptr = '\0'; break;
			case 'O':
				if (!strncmp (ptr, "333", 3)) {
					if (strchr (ptr, '.')) yshift = atoi (ptr + 4);
									odim_flag = 3;
				}
				else if (!strncmp (ptr, "222", 3))	odim_flag = 2;
				else if (!strncmp (ptr, "111", 3))	odim_flag = 1;
				else {
					getroot (ptr, dimfile);		odim_flag = 0;
				}
				*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	strncpy (t4file,  argv[i], MAXL);	k++; break;
			case 1:	getroot (argv[i], imgroot);		k++; break;
		}	
	}
	if (k < 2) {
		fprintf (stderr, "Usage:\t%s [options] <t4file> <input_4dfp_stack>\n", program);
		fprintf (stderr, " e.g.,\t%s -aatl anat_ave_to_711-2B_t4 b1_rmsp_dbnd\n", program);
		fprintf (stderr, "\toption\n");
		fprintf (stderr, "\t-a<str>\tspecify outfile name trailer (default = \"xr3d\")\n");
		fprintf (stderr, "\t-c<flt>\tscale output by specified factor\n");
		fprintf (stderr, "\t-N\toutput undefined voxels as NaN\n");
		fprintf (stderr, "\t-Z\toutput undefined voxels as 0\n");
		fprintf (stderr, "\t-E\toutput undefined voxels as 1.e-37 (default)\n");
		fprintf (stderr, "\t-v[0|1]\tset per frame intensity equalization mode (default = OFF)\n");
		fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default input endian)\n");
		fprintf (stderr, "\t-f\tfast (linear interpolation resample instead of 3D cubic spline)\n");
		fprintf (stderr, "\t-e\techo mat file to stdout frame by frame (verbose mode)\n");
		fprintf (stderr, "\t-O111\toutput in 111 space\n");
		fprintf (stderr, "\t-O222\toutput in 222 space\n");
		fprintf (stderr, "\t-O333.n\toutput in 333.n space (y shifted up by n pixels)\n");
		fprintf (stderr, "\t-O<str>\toutput image dimensions according to <str>.4dfp.ifh\n");
		fprintf (stderr, "N.B.:\tdefault output format = 333.%d\n", YSHIFT);
		exit (1);
	}

	switch (odim_flag) {
	case 0:	strcat (dimfile, ".4dfp");
		if (Getifh (dimfile, &ifh)) exit (-1);
		outdim[0] = ifh.matrix_size[0];
		outdim[1] = ifh.matrix_size[1];
		outdim[2] = ifh.matrix_size[2];
		for (k = 0; k < 3; k++) {
			mmppixo[k] = ifh.mmppix[k];
			centero[k] = ifh.center[k];
		}
		orio = ifh.orientation;
		break;
	case 1:	outdim[0] = 176;
		outdim[1] = 208;
		outdim[2] = 176;
		mmppixo[0] =  1.0;
		mmppixo[1] = -1.0;
		mmppixo[2] = -1.0;
		centero[0] =   89.;
		centero[1] =  -85.;
		centero[2] = -101.;
		break;
	case 2:	outdim[0] = 128;
		outdim[1] = 128;
		outdim[2] = 75;
		mmppixo[0] =  2.0;
		mmppixo[1] = -2.0;
		mmppixo[2] = -2.0;
                centero[0] =  129.;
		centero[1] = -129.;
		centero[2] = mmppixo[2] * 41;
		break;
	case 3:	outdim[0] = 48;
		outdim[1] = 64;
		outdim[2] = 48;
		mmppixo[0] =  3.0;
		mmppixo[1] = -3.0;
		mmppixo[2] = -3.0;
		centero[0] =  73.5;				/* 24.5 * +3.0 */
		centero[1] =  (29 - yshift) * mmppixo [2];
		centero[2] = -84.0;				/* 28.0 * -3.0 */
		break;
	}
	odim = outdim[0] * outdim[1] * outdim[2];		/* output voxel count */

/****************/
/* read t4 file */
/****************/
	if (strcmp (t4file, "none")) {
		if (access (t4file, R_OK)) {
			fprintf (stderr, "%s: %s read error\n", program, t4file);
			exit (-1);
		}
		printf ("Reading: %s\n", t4file);
		t4_read_ (t4file, t4atl);
	} else {
		t4_init_ (t4atl);
	}
	
/****************************/
/* assemble image filenames */
/****************************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (ptr = strrchr (imgroot, '/')) ptr++; else ptr = imgroot;
	sprintf (matfile, "%s_xr3d.mat", imgroot);
	sprintf (outroot, "%s_%s", ptr, trailer);
	sprintf (outfile, "%s.4dfp.img", outroot);

/*****************************/
/* get 4dfp stack dimensions */
/*****************************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0
	|| Getifh (imgfile, &ifh)) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';
	vdim = imgdim[0] * imgdim[1] * imgdim[2];
	outdim[3] = imgdim[3];

/***************************************/
/* prepare virtual flip of input frame */
/***************************************/
	vrtflip_ (&ifh.orientation, imgdim, ifh.center, ifh.mmppix, centert, mmppixt);

	printf ("input orient=%d\n", ifh.orientation);
	printf ("input stack dimensions %10d%10d%10d%10d\n",    imgdim[0],  imgdim[1],  imgdim[2], imgdim[3]);
	printf ("input stack mmppix     %10.6f%10.6f%10.6f\n", mmppixt[0], mmppixt[1], mmppixt[2]);
	printf ("input stack center     %10.4f%10.4f%10.4f\n", centert[0], centert[1], centert[2]);
	printf ("output orient=%d\n", orio);
	printf ("output stack dimensions%10d%10d%10d%10d\n",    outdim[0],  outdim[1],  outdim[2], outdim[3]);
	printf ("output stack mmppix    %10.6f%10.6f%10.6f\n", mmppixo[0], mmppixo[1], mmppixo[2]);
	printf ("output stack center    %10.4f%10.4f%10.4f\n", centero[0], centero[1], centero[2]);

/**************************/
/* allocate image  memory */
/**************************/
	if (!(imgt = (float *) malloc (4 * vdim * sizeof (float))))	errm (program);
	if (!(imgo = (float *) calloc (odim, sizeof (float))))		errm (program);

	if (!(matfp = fopen (matfile, "r")))  errr (program, matfile);
	if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!(outfp = fopen (outfile, "wb"))) errw (program, outfile);

	printf ("Reading: %s\n", imgfile);
	printf ("Writing: %s frame", outfile);
	for (i = 0; i < imgdim[3]; i++) {
		fscanf (matfp, "%*s %*s %d", &k); if (k != i + 1) errr (program, matfile);
		for (k = 0; k < 4; k++)	fscanf (matfp, "%f %f %f %f", t4mat + k + 0, t4mat + k + 4, t4mat + k + 8, t4mat + k + 12);
		fscanf (matfp, "%*s %*s %d", &k); if (k != i + 1) errr (program, matfile);
		fscanf (matfp, "%f %f %f %f", s4 + 0, s4 + 1, s4 + 2, s4 + 3);
		printf (" %d", i + 1); fflush (stdout);
		if (verbose) {
			t4_list_ (t4mat);
			printf ("s4\n%10.6f%10.6f%10.6f%10.6f\n", s4[0], s4[1], s4[2], s4[3]);
		}
		if (eread (imgt, vdim, isbig, imgfp)) errr (program, imgfile);

/*********************/
/* execute transform */
/*********************/
		pt4ixyz_ (imgdim, voxdim, centert, mmppixt, t4mat, centero, mmppixo, t4atl, t4);
		if (fast_flag) {
			ft4imgo_ (t4, imgt, imgdim+0, imgdim+1, imgdim+2, mmppix1, mmppix1,
				      imgo, outdim+0, outdim+1, outdim+2, center0, mmppix1);
		} else {
			spline3d (imgt, imgdim);
			ft4ixyz_ (&mode, t4, imgt, imgdim+0, imgdim+1, imgdim+2, imgo, outdim+0, outdim+1, outdim+2);
		}

/*******************************/
/* restore 4dfp voxel indexing */
/*******************************/
		x4dfp2ecat (imgo, outdim, orio);

		q = (volnorm) ? s4[0] : 1.0;
		for (k = 0; k < odim; k++) {
			if (!isnan (imgo[k])) {
				imgo[k] *= cscale * q;
				if (fmin > imgo[k]) fmin = imgo[k];
				if (fmax < imgo[k]) fmax = imgo[k];
			} else switch (NaN_flag) {
				case 'Z': imgo[k] = 0.0;		break;
				case 'E': imgo[k] = (float) 1.e-37;	break;
				case 'N': default:			break;
			}
		}
		if (ewrite (imgo, odim, control, outfp)) errw (program, outfile);
	}
	printf ("\n");
	if (fclose (imgfp)) errr (program, imgfile);
	if (fclose (matfp)) errr (program, matfile);
	if (fclose (outfp)) errw (program, outfile);

/*******************/
/* create rec file */
/*******************/
	startrecle (outfile, argc, argv, rcsid, control);
	if (strcmp (t4file, "none")) catrec (t4file);
	catrec (imgfile);
	endrec ();

	voxdim[0] = (mmppixo[0] > 0.) ? mmppixo[0]: -mmppixo[0];
	voxdim[1] = (mmppixo[1] > 0.) ? mmppixo[1]: -mmppixo[1];
	voxdim[2] = (mmppixo[2] > 0.) ? mmppixo[2]: -mmppixo[2];
/*******/
/* ifh */
/*******/
	writeifhmce (program, outfile, outdim, voxdim, orio, mmppixo, centero, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s -r%dto%d", outroot, (int) (fmin - 0.5), (int) (fmax + 0.5));
	printf ("%s\n", command); status = system (command);

	free (imgo); free (imgt);
	exit (0);
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
