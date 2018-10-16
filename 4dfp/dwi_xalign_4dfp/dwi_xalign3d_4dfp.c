/****************************************************************************************/
/* Copyright 2003, 2004, 2005, 2006, 2007, 2008, 2009, 1010, 2011			*/
/* Washington University, Mallinckrodt Institute of Radiology.				*/
/* All Rights Reserved.									*/
/* This software may not be reproduced, copied, or distributed without written		*/
/* permission of Washington University. For further information contact A. Z. Snyder.	*/
/****************************************************************************************/
/*$Header: /data/petsun4/data1/src_solaris/dwi_xalign_4dfp/RCS/dwi_xalign3d_4dfp.c,v 1.23 2011/03/09 05:33:56 avi Exp $*/
/*$Log: dwi_xalign3d_4dfp.c,v $
 * Revision 1.23  2011/03/09  05:33:56  avi
 * *** empty log message ***
 *
 * Revision 1.22  2011/03/09  05:29:29  avi
 * correct usage
 *
 * Revision 1.21  2009/07/14  18:59:24  avi
 * correct absence of "-m" functionality description in Usage
 *
 * Revision 1.20  2009/02/20  01:17:55  avi
 * option -I
 *
 * Revision 1.19  2008/12/30  05:29:02  avi
 * more missing prototypes
 *
 * Revision 1.18  2008/12/30  05:07:56  avi
 * always compute group mean (eliminate option -m)
 *
 * Revision 1.17  2008/08/15  03:49:03  avi
 * proper prototypes
 *
 * Revision 1.16  2008/08/15  02:06:22  avi
 * linux compliant (remove f_init () and f_exit ())
 *
 * Revision 1.15  2008/07/01  03:00:46  avi
 * pass 2D vs. 3D mode using modes[] instead of &mo3d in calls to resample_()
 *
 * Revision 1.14  2008/03/13  00:46:48  avi
 * option -p (2D alignment mode)
 *
 * Revision 1.13  2007/06/21  01:03:36  avi
 * MAXG -> 48; fwrite -> gwrite (endian compliant)
 *
 * Revision 1.12  2007/03/13  19:31:59  avi
 * option -a (compute arithmetic as opposed to geometric mean group average)
 * rec file for _geom (output group average image)
 *
 * Revision 1.11  2007/03/12  02:45:51  avi
 * Solaris 10
 *
 * Revision 1.10  2005/12/29  07:34:09  avi
 * externalize to command line control over parameter search radius (fimgeta_() -> fimgetae_())
 * get_4d_dimensions() -> get_4dfp_dimo()
 *
 * Revision 1.9  2004/12/23  00:26:28  avi
 * accommodate possibility of > 99 volumes in ASCII output
 *
 * Revision 1.8  2004/12/21  05:51:35  avi
 * MAXV -> 364 (=7*52)
 *
 * Revision 1.7  2004/12/14  01:27:43  avi
 * msk0 computation accepts 1.0 as in mask
 *
 * Revision 1.6  2003/10/03  21:29:59  avi
 * -n option (zero negative voxels in output)
 *
 * Revision 1.5  2003/05/04  00:09:28  avi
 * two computations of group geometric mean in iter loop
 *
 * Revision 1.4  2003/05/01  02:14:40  avi
 * correct matrix composition operations
 *
 * Revision 1.3  2003/04/30  22:47:52  avi
 * eliminate temp file
 *
 * Revision 1.2  2003/04/30  03:14:23  avi
 * -d option
 * correct computation of msk2
 *
 * Revision 1.1  2003/04/29  23:57:45  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <rec.h>
#include <endianio.h>
#include <Getifh.h>
#include <spline23dvgh.h>
#include <librms.h>

#define MAXL		256
#define MAXV		364		/* maximum volumes in dwi_data image */
#define MAXG		48		/* maximum number of alignment groups */
#define NCYCLEG		3		/* group alignment cycles */
#define FSMOOTH		0.0		/* pre-filter */
#define DEL		5.0		/* sampling interval in mm */
#define EPSMM		3.0		/* fimgetae_() translation parameter search radius in mm */
#define EPSRAD		40.0		/* fimgetae_() eps head radius in mm */

extern void	t4_init_ (float *t4);					/* fomega.f */
extern void	t4_omega_ (long *modes, float *t4, float *omega);	/* fomega.f */
extern void	fimgetae_ (float *img1, float *d2x1, float *d2y1, float *d2z1, int *nx1, int *ny1, int *nz1,
		float *mmppix1, float *center1, short *msk1, long *mode, float *del, float *eps,
			   float *img2, float *d2x2, float *d2y2, float *d2z2, int *nx2, int *ny2, int *nz2,
		float *mmppix2, float *center2, short *msk2, float *t4);
extern void	etagen_   (float *img1, float *d2x1, float *d2y1, float *d2z1,
			   int *nx1, int *ny1, int *nz1, float *mmppix1, float *center1, short *msk1,
			   long *mode, float *omega, float *del,
			   float *img2, float *d2x2, float *d2y2, float *d2z2,
			   int *nx2, int *ny2, int *nz2, float *mmppix2, float *center2, short *msk2,
			   float *eta, float *q);							/* etadeta.f */
extern void	resample_ (float *img1, float *d2x1, float *d2y1, float *d2z1, int *nx1, int *ny1, int *nz1,
		float *mmppix1, float *center1, short *msk1, long *mode, float *t4,
			   float *img2, int *nx2, int *ny2, int *nz2, float *mmppix2, float *center2);	/* etadeta.f */
extern void	gauss3d_ (float *image, int *pnx, int *pny, int *pnz, float *cmppix, float *fhalf);	/* gauss3d.f */
extern void	imag2mask ();						/* in librms */	
void	filter3d (char *program, float *imgt, int *imgdim, float *voxdim, float f0, int wrap);	/* below */
void	spline3d (char *program, float *imgf, int *imgdim);					/* below */
void	volalign (float *tar, short *tarm, char *tarstr, float *src, short *srcm, char *srcstr,	/* below */
		  long *modes, float *t4, float *peta, float *pq);
void	load_filter_spline (char *srcroot, int isbig, int ivol, float *imag);			/* below */
void	group_mean (char *srcroot, int isbig, int group[], char *grpstr, int flag);		/* below */
int	ingroup (int iquery, int group[MAXV]);							/* below */

typedef struct {
	float	t4[16];
	float	eta;
	float	q;
	float	det;
} VOLDAT;

/********************/
/* global variables */
/********************/
char		program[MAXL];
long		mo3d = 259;			/* enable 3D cross_modal */
long		modes[8];			/* this program uses only modes[0] */
int		wrap = 0;			/* flag used by filter3d */
float		fsmooth = FSMOOTH;		/* 3D pre-filter half freq (1/mm) */
float		del = DEL;			/* sampling interval */
float		eps[12];			/* fimgetae_() parameter search radii */
float		mmppixs[3], centers[3];
int		srcsiz, srcdim[4];
float		*imgs, *img0, *img1, *img2;
short		*msk0, *msk1, *msk2, *all1;
VOLDAT		*voldat;
static char rcsid[] = "$Id: dwi_xalign3d_4dfp.c,v 1.23 2011/03/09 05:33:56 avi Exp $";

static void gparse (char *string, int group[MAXV]) {
	char	*ptr, *ptr1;
	int	i, j, k, i0, i1;

	strcat (string, ",");
	j = i = 0; while (ptr = strtok (string + i, ",")) {
		i += strlen (ptr) + 1;
		if (0) printf ("%5d\t%s\t", i, ptr);
		if (ptr1 = strchr (ptr, '-')) {
			i1 = atoi (ptr1 + 1); *ptr1 = '\0';
			i0 = atoi (ptr);
		} else {
			i1 = i0 = atoi (ptr);
		}
		for (k = i0; k <= i1; k++) {
			printf ("%d ", k);
			group[j++] = k - 1;
		}
	}
	printf ("\n");
}

static void usage () {
	printf ("Usage:\t%s <(4dfp) dwi> <(4dfp) mask>\n", program);
	printf (" e.g.:\t%s hbo08a_dwi1 hbo08a_dwi1_mskt -s -g2-4 -g5,13,18,23\n", program);
	printf ("\toption\n");
	printf ("\t-p\tplanar (2D; disable cross-slice) alignment\n");
	printf ("\t-w\tenable wrap addressing\n");
	printf ("\t-s\tenable cross DWI voxel size adjust (principal axis stretch)\n");
	printf ("\t-a\tcompute group arithmeric mean volume (default geometric mean)\n");
	printf ("\t-n\tzero negative values in output image\n");
	printf ("\t-I<int>\tspecify volume number of I0 counting from 1 (default 1)\n");
	printf ("\t-f<flt>\tspecify pre-blur filter half freq (1/mm) (default none)\n");
	printf ("\t-d<flt>\tspecify sampling interval in mm (default=%.4f)\n", DEL);
	printf ("\t-i<flt>\tspecify displacment search radius in mm (default=%.4f)\n", EPSMM);
	printf ("\t-j<flt>\tspecify parameter search object radius in mm (default=%.4f)\n", EPSRAD);
	printf ("\t-c<int>\tspecify number of within-group cycles (default=%d)\n", NCYCLEG);
	printf ("\t-g<int>[-<int>][,<int>[-<int>]][,...]\tprogram alignment group\n");
	printf ("N.B.:\t<(4dfp) mask> may be \"none\"\n");
	printf ("N.B.:\tI0 should not be named in any programmed alignent group\n");
	exit (1);
}

int main (int argc, char *argv[]) {
/************/
/* 4dfp I/O */
/************/
	FILE		*fp_out, *fp_mat, *fp_dat, *fp_rec, *fp_grp;
	char		mskroot[MAXL], mskfile[MAXL], srcroot[MAXL], srcfile[MAXL];
	char		outfile[MAXL], grpfile[MAXL], matfile[MAXL], datfile[MAXL], recfile[MAXL];

/****************/
/* image arrays */
/****************/
	float		mmppixm[3];
	int		mskdim[4], grpdim[4], orientm, orients, isbigs, isbigm;
	char		control = '\0';

/***********/
/* utility */
/***********/
	char		*str, string[MAXL], volstr[MAXL], grpstr[MAXL];
	int		c, i, j, k, m, iter, ivol, four = 4;
	int		nzero = 0;
	double		q1, q2;

/**************/
/* processing */
/**************/
	float		t4[16];
	float		q, eta;
	float		epsmm = EPSMM, epsrad = EPSRAD;
	VOLDAT		grpdat[MAXG];

/********************/
/* alignment groups */
/********************/
	int		groups[MAXG][MAXV];
	int		igroup, ngroup = 0;
	int		ncycleg = NCYCLEG;
	int		I0vol = 0;

/*********/
/* flags */
/*********/
	int		planar = 0;
	int		stretch = 0;
	int		status = 0;
	int		debug = 0;
	int		save_grp = 1;
	int		nflag = 0;	/* zero negative output */

	printf ("%s\n", rcsid);
	if (!(str = strrchr (argv[0], '/'))) str = argv[0]; else str++;
	strcpy (program, str);
	for (i = 0; i < MAXG; i++) for (j = 0; j < MAXV; j++) groups[i][j] = -1;

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (string, argv[i]); str = string;
			while (c = *str++) switch (c) {
				case 'a': save_grp = 2;			break;
				case 'p': planar++;			break;
				case 's': stretch++;			break;
				case 'w': mo3d |= 1024; wrap++; 	break;
				case 'n': nflag++;			break;
				case 'I': I0vol = atoi (str) - 1;	*str = '\0'; break;
				case 'c': ncycleg = atoi (str);		*str = '\0'; break;
				case 'd': del = atof (str);		*str = '\0'; break;
				case 'f': fsmooth = atof (str);		*str = '\0'; break;
				case 'i': epsmm = atof (str);		*str = '\0'; break;
				case 'j': epsrad = atof (str);		*str = '\0'; break;
				case 'g': printf ("%s: volumes in alignment group %d\n", program, ngroup + 1);
					gparse (str, groups[ngroup++]);
					if (ngroup == MAXG) {
						fprintf (stderr, "%s: not configured for more than %d alignment groups\n",
								program, MAXG);
						exit (-1);
					}
									*str = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	getroot (argv[i], srcroot); k++; break;
			case 1:	getroot (argv[i], mskroot); k++; break;
		}	
	}
	if (k < 2) usage ();
	for (igroup = 0; igroup < ngroup; igroup++) for (j = 0; j < MAXV; j++) if (groups[igroup][j] == I0vol) usage ();

/***********************/
/* set fimgetae_() eps */
/***********************/
	for (i = 0; i < 3;  i++) eps[i] = epsmm;
	for (i = 3; i < 12; i++) eps[i] = epsmm/epsrad;

	if (get_4dfp_dimoe (srcroot, srcdim, mmppixs, &orients, &isbigs)) errr (srcroot, program);
	if (!control) control = (isbigs) ? 'b' : 'l';
	if (srcdim[3] > MAXV) {
		fprintf (stderr, "%s: not configured for more than %d DWI volumes\n", program, MAXV);
		exit (-1);
	}
	for (srcsiz = 1, k = 0; k < 3; k++) {
		srcsiz *= srcdim[k];
		centers[k] = mmppixs[k] * (srcdim[k] / 2);
		grpdim[k] = srcdim[k];
	}
	if (!(img2 = (float *) malloc (srcsiz * sizeof (float) * 4)))	errm (program);
	if (!(img1 = (float *) malloc (srcsiz * sizeof (float) * 4)))	errm (program);
	if (!(img0 = (float *) malloc (srcsiz * sizeof (float) * 4)))	errm (program);	/* I0 */
	if (!(imgs = (float *) malloc (srcsiz * sizeof (float))))	errm (program);	/* DWI volume */
	if (!(msk0 = (short *) malloc (srcsiz * sizeof (short))))	errm (program);
	if (!(msk1 = (short *) malloc (srcsiz * sizeof (short))))	errm (program);
	if (!(msk2 = (short *) malloc (srcsiz * sizeof (short))))	errm (program);
	if (!(all1 = (short *) malloc (srcsiz * sizeof (short))))	errm (program);
	for (i = 0; i < srcsiz; i++) all1[i] = 1;

/************/
/* set mask */
/************/
        if (strstr (mskroot, "none")) {
		for (i = 0; i < srcsiz; i++) msk0[i] = 1;
	} else {
		if (get_4dfp_dimoe (mskroot, mskdim, mmppixm, &orientm, &isbigm)) errr (mskroot, program);
		status = orientm - orients;
		for (k = 0; k < 3; k++) status |= (mskdim[k] != srcdim[k] || fabs (mmppixs[k] - mmppixm[k]) > 1.e-4);
		if (status) {
			fprintf (stderr, "%s: %s %s dimension mismatch\n", program, srcroot, mskroot);
			exit (-1);
		}
		load_4dfp_frame (mskroot, srcdim, 0, isbigm, imgs);
		for (i = 0; i < srcsiz; i++) msk0[i] = (imgs[i] >= 1.0) ? 1 : 0;
	}

/**************************************/
/* allocate and initialize transforms */
/**************************************/
	if (!(voldat = (VOLDAT *) calloc (srcdim[3], sizeof (VOLDAT)))) errm (program);
	for (ivol = 0; ivol < srcdim[3]; ivol++) t4_init_ (voldat[ivol].t4);

/************************************************/
/* compute first pass cross-encoding alignments */
/************************************************/
	i = 0; modes[i++] = mo3d;				/* no stretch in first call to fimgeta */
	if (stretch) modes[i++] = mo3d | (16 + 32 + 64);	/* enable xyz stretch */
	modes[i++] = 0;
	if (planar) for (i = 0; i < 8; i++) modes[i] &= ~(2 + 64);	/* turn off bits 2 and 64 */
	printf ("modes"); for (i = 0; i < 8; i++) printf (" %x", modes[i]); printf ("\n");

	load_filter_spline (srcroot, isbigs, I0vol, img0);
	for (ivol = 0; ivol < srcdim[3]; ivol++) {
		if (ivol == I0vol) continue;
		sprintf (volstr, "volume %3d", ivol + 1);
		load_filter_spline (srcroot, isbigs, ivol, img1);
		volalign (img0, msk0, "I0", img1, all1, volstr, modes, voldat[ivol].t4, &voldat[ivol].eta, &voldat[ivol].q);
	}

/*************************/
/* group-wise alignments */
/*************************/
	if (save_grp && ngroup && ncycleg > 0) {
		sprintf (grpfile, "%s_geom.4dfp.img", srcroot);
		if (!(fp_grp = fopen (grpfile, "w+b"))) errw (program, grpfile);
		printf ("Writing: %s\n", grpfile);
	}

	i = 0; modes[i++] = (stretch) ? mo3d | (16 + 32 + 64) : mo3d;
	modes[i++] = 0;
	if (planar) for (i = 0; i < 8; i++) modes[i] &= ~(2 + 64);	/* turn off bits 2 and 64 */
	printf ("modes"); for (i = 0; i < 8; i++) printf (" %x", modes[i]); printf ("\n");

for (igroup = 0; igroup < ngroup; igroup++) {t4_init_ (grpdat[igroup].t4);
for (iter = 0; iter < ncycleg; iter++) {
	printf ("%s: aligning group %d cycle %d\n", program, igroup + 1, iter + 1);

/*****************************/
/* create group mean in img2 */
/*****************************/
	sprintf (grpstr, "group %d mean cycle %d", igroup + 1, iter + 1);
	group_mean (srcroot, isbigs, groups[igroup], grpstr, save_grp);
	if ((int) fsmooth) filter3d (program, img2, srcdim, mmppixs, fsmooth, wrap);
	spline3d (program, img2, srcdim);

/*********************************/
/* realign volumes to group mean */
/*********************************/
	for (ivol = 0; ivol < srcdim[3]; ivol++) if (ingroup (ivol, groups[igroup])) {
		sprintf (volstr, "volume %3d", ivol + 1);
		load_filter_spline (srcroot, isbigs, ivol, img1);
		volalign (img2, msk2, grpstr, img1, all1, volstr, modes,
			voldat[ivol].t4, &voldat[ivol].eta, &voldat[ivol].q);
	}

/********************************/
/* re-create group mean in img2 */
/********************************/
	group_mean (srcroot, isbigs, groups[igroup], grpstr, save_grp);
	if (save_grp) {
		if (fseek (fp_grp, (long) igroup * srcsiz * sizeof (float), SEEK_SET)
		||  gwrite ((char *) img2, sizeof (float), srcsiz, fp_grp, control)) errw (program, grpfile);
	}
	if ((int) fsmooth) filter3d (program, img2, srcdim, mmppixs, fsmooth, wrap);
	spline3d (program, img2, srcdim);
		
/*********************************/
/* align group mean (img2) to I0 */
/*********************************/
	volalign (img0, msk0, "I0", img2, msk2, grpstr, modes, grpdat[igroup].t4, &grpdat[igroup].eta, &grpdat[igroup].q);

/*************************************************/
/* compose volume->I0 and grpmean->I0 transforms */
/*************************************************/
	for (ivol = 0; ivol < srcdim[3]; ivol++) if (ingroup (ivol, groups[igroup])) {
		matmul_ (voldat[ivol].t4, grpdat[igroup].t4, t4, &four);
		for (i = 0; i < 16; i++) voldat[ivol].t4[i] = t4[i];
	}
	fflush (stdout);
}	/* iter */
}	/* igroup */

/*************************/
/* close group mean file */
/*************************/
	if (save_grp && (grpdim[3] = ngroup) && ncycleg > 0) {
		fclose (fp_grp);
		writeifhe (program, grpfile, grpdim, mmppixs, orients, control);
		sprintf (string, "ifh2hdr %s -r2000", grpfile); system (string);
	}
	startrecle (grpfile, argc, argv, rcsid, control);
	sprintf (string, "Group %s mean image[s]\n", (save_grp == 1) ? "geometric" : "arithmetic"); printrec (string);
	for (igroup = 0; igroup < ngroup; igroup++) {
		sprintf (string, "Alignment group %d volumes (counting from 1):", igroup + 1); printrec (string);
		for (ivol = 0; ivol < srcdim[3]; ivol++) if (ingroup (ivol, groups[igroup])) {
			sprintf (string, " %d", ivol + 1); printrec (string);
		}
		printrec ("\n");
	}
	sprintf (srcfile, "%s.4dfp.img", srcroot); catrec (srcfile);
	endrec ();

/********************************/
/* start cross-encoding matfile */
/********************************/
	sprintf (matfile, "%s_xenc.mat", srcroot);
	if (!(fp_mat = fopen (matfile, "w"))) errw (program, matfile);

/*********************************************************/
/* apply compound transforms and write realigned volumes */
/*********************************************************/
	sprintf (outfile, "%s_xenc.4dfp.img", srcroot);
	printf ("Writing: %s\n", outfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);

	for (ivol = 0; ivol < srcdim[3]; ivol++) {
		load_4dfp_frame (srcroot, srcdim, ivol, isbigs, img1);
		spline3d (program, img1, srcdim);
		j = srcsiz;
		resample_ (img1,img1+1*j,img1+2*j,img1+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,all1,modes,
	                              voldat[ivol].t4,imgs,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers);
		determ12_ (voldat[ivol].t4, &voldat[ivol].det);
		for (q1 = q2 = m = i = 0; i < srcsiz; i++) {
			q1 += img1[i];
			imgs[i] *= voldat[ivol].det;
			if (imgs[i] == 0.0) {
				m++;
			} else {
				q2 += imgs[i];
			}
			if (nflag && imgs[i] < 0.0) {imgs[i] = 0.0; nzero++;}
		}
		printf ("image mean %10.2f -> %10.2f 0 voxel count=%d\n", q1/srcsiz, q2/(srcsiz-m), m);
		fprintf (fp_mat, "t4 contrast %3d\n", ivol + 1);
		for (j = 0; j < 4; j++) {
			fprintf (fp_mat, "%10.6f%10.6f%10.6f%10.4f\n",
			voldat[ivol].t4[0+j], voldat[ivol].t4[4+j], voldat[ivol].t4[8+j], voldat[ivol].t4[12+j]);
		}
		if (gwrite ((char *) imgs, sizeof (float), srcsiz, fp_out, control)) errw (program, outfile);
	}
	fclose (fp_mat);
	fclose (fp_out);

/***********************/
/* outfile ifh and hdr */
/***********************/
	sprintf (string, "/bin/cp %s.4dfp.ifh %s_xenc.4dfp.ifh", srcroot, srcroot); system (string);
	sprintf (string, "/bin/cp %s.4dfp.hdr %s_xenc.4dfp.hdr", srcroot, srcroot); system (string);

/************/
/* dat file */
/************/
	sprintf (datfile, "%s_xenc.4dfp.dat", srcroot);
        if (!(fp_dat = fopen (datfile, "w"))) errw (program, datfile);
	printf ("Writing: %s\n", datfile);
	fprintf (fp_dat, "#   volume       eta         q       det\n");
	for (ivol = 0; ivol < srcdim[3]; ivol++) {
		fprintf (fp_dat, "%10d%10.6f%10.6f%10.6f\n", ivol + 1,
		voldat[ivol].eta, voldat[ivol].q, voldat[ivol].det);
	}
	fclose (fp_dat);

/***********/
/* recfile */
/***********/
	startrecle (outfile, argc, argv, rcsid, control);
	sprintf (recfile, "%s.rec", outfile);
	if (!(fp_rec = fopen (recfile, "a"))) errw (program, recfile);
	fprintf (fp_rec, "I0 volume = %d (counting from 1)\n", I0vol);
	fprintf (fp_rec, "Alignment pre-filter fhalf = %.4f (1/mm)\n", fsmooth);
	fprintf (fp_rec, "Sampling interval = %.4f mm\n", del);
	fprintf (fp_rec, "Translation parameter search radius = %.4f mm\n", epsmm);
	fprintf (fp_rec, "Parameter search object radius = %.4f mm\n", epsrad);
	for (igroup = 0; igroup < ngroup; igroup++) {
		fprintf (fp_rec, "Alignment group %d volumes (counting from 1):", igroup + 1);
		for (i = 0; i < srcdim[3]; i++) if (ingroup (i, groups[igroup])) fprintf (fp_rec, " %d", i + 1);
		fprintf (fp_rec, "\n");
	}
	if (nflag) fprintf (fp_rec, "%d output voxels zeroed to eliminate negative values\n", nzero);
	fclose (fp_rec);
	catrec (datfile);
	sprintf (srcfile, "%s.4dfp.img", srcroot); catrec (srcfile);
	catrec (mskfile);
	endrec ();

/************/
/* clean up */
/************/
	free (voldat);
	free (imgs); free (img2); free (img1); free (img0);
	free (msk1); free (msk2); free (msk0); free (all1);
	exit (status);
}

void filter3d (char *program, float *imgt, int *imgdim, float *mmppix, float f0, int wrap) {
	float		val, *imgp;
	int		k, margin, nxp, nyp, nzp;

	if (wrap) {
		nxp = imgdim[0];
		nyp = imgdim[1];
	} else {
		val = (0.5 + (2.0 * 0.1874 / (mmppix[0] * f0))); margin = val; nxp = npad_ (imgdim + 0, &margin);
		val = (0.5 + (2.0 * 0.1874 / (mmppix[1] * f0))); margin = val; nyp = npad_ (imgdim + 1, &margin);
	}
		val = (0.5 + (4.0 * 0.1874 / (mmppix[2] * f0))); margin = val; nzp = npad_ (imgdim + 2, &margin);
	printf ("image dimensions %d %d %d padded to %d %d %d \n", imgdim[0], imgdim[1], imgdim[2], nxp, nyp, nzp);
        if (!(imgp = (float *) calloc (nxp * nyp * nzp, sizeof (float)))) errm (program);
	imgpad_ (imgt, imgdim + 0, imgdim + 1, imgdim + 2, imgp, &nxp, &nyp, &nzp);
	gauss3d_ (imgp, &nxp, &nyp, &nzp, mmppix, &f0);
	imgdap_ (imgt, imgdim + 0, imgdim + 1, imgdim + 2, imgp, &nxp, &nyp, &nzp);
	free (imgp);
}

void spline3d (char *program, float *imgf, int *imgdim) {
	float		*imgp;
	int		k, nzt, four = 4, vdim;

	vdim = 1; for (k = 0; k < 3; k++) vdim *= imgdim[k];
	splinex_  (imgf, imgdim + 0, imgdim + 1, imgdim + 2, imgf + vdim*1);
	spliney_  (imgf, imgdim + 0, imgdim + 1, imgdim + 2, imgf + vdim*2);
	if (imgdim[2] < 48) {
		splineza_ (imgf,  imgdim + 0, imgdim + 1, imgdim + 2, imgf + vdim*3);
	} else {
		nzt = npad_ (imgdim + 2, &four);
		k = imgdim[0] * imgdim[1] * nzt;
		if (!(imgp = (float *) malloc (k * sizeof (float) * 2))) errm (program);
		imgpad_   (imgf,          imgdim + 0, imgdim + 1, imgdim + 2, imgp,     imgdim + 0, imgdim + 1, &nzt);
		splinezf_ (imgp,          imgdim + 0, imgdim + 1, &nzt,       imgp + k);
		imgdap_   (imgf + vdim*3, imgdim + 0, imgdim + 1, imgdim + 2, imgp + k, imgdim + 0, imgdim + 1, &nzt);
		free (imgp);
	}
}

void volalign (float *tar, short *tarm, char *tarstr, float *src, short *srcm, char *srcstr,
		long *modes, float *t4, float *peta, float *pq) {
	int		i, j, k, msrc, mtar;
	float		omega[12];

	printf ("aligning %s to %s\n", srcstr, tarstr);
	for (msrc = mtar = i = 0; i < srcsiz; i++) {
		if (tarm[i]) mtar++;
		if (srcm[i]) msrc++;
	} 
	printf ("target voxels=%d\tsource voxels=%d\n", mtar, msrc);

	j = srcsiz;
	i = 0; while (modes[i]) {
		fimgetae_ (tar,tar+1*j,tar+2*j,tar+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,tarm,modes+i,&del,eps,
			   src,src+1*j,src+2*j,src+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,srcm,t4);
		i++;
	} if (i) i--;
	t4_omega_(modes+i, t4, omega);
	etagen_  (tar,tar+1*j,tar+2*j,tar+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,tarm,modes+i,omega,&del,
		  src,src+1*j,src+2*j,src+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,srcm,peta,pq);
}

void load_filter_spline (char *srcroot, int isbig, int ivol, float *imag) {
	load_4dfp_frame (srcroot, srcdim, ivol, isbig, imag);
	if ((int) fsmooth) filter3d (program, imag, srcdim, mmppixs, fsmooth, wrap);
	spline3d (program, imag, srcdim);
}

/*****************************/
/* create group mean in img2 */
/*****************************/
void group_mean (char *srcroot, int isbig, int group[], char *grpstr, int flag) {
	int	ivol, i, j;
	double	q;

	printf ("computing: %s\n", grpstr);
	for (i = 0; i < srcsiz; i++) msk2[i] = img2[i] = 0.;
	for (ivol = 0; ivol < srcdim[3]; ivol++) if (ingroup (ivol, group)) {
		load_4dfp_frame (srcroot, srcdim, ivol, isbig, img1);
		spline3d (program, img1, srcdim);
		j = srcsiz;
		resample_ (img1,img1+1*j,img1+2*j,img1+3*j,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers,all1,modes,
	                              voldat[ivol].t4,imgs,srcdim+0,srcdim+1,srcdim+2,mmppixs,centers);
		for (i = 0; i < srcsiz; i++) if (imgs[i] > 1.0) {
			switch (flag) {
				case 1: q = log (imgs[i]);	break;
				case 2: q = imgs[i];		break;
			}
			img2[i] += q;
			msk2[i]++;
		}
	}
	for (i = 0; i < srcsiz; i++) {
		if (!msk2[i]) {
			q = 0;
		} else {
			switch (flag) {
				case 1: q = exp (img2[i]/msk2[i]);	break;
				case 2: q =      img2[i]/msk2[i];	break;
			}
		}
		img2[i] = q;
 		msk2[i] *= msk0[i];	/* msk0[i] is 1 or 0 */
	}
}

int ingroup (int iquery, int group[MAXV]) {
	int	i;

	for (i = 0; i < MAXV; i++) if (iquery == group[i]) return 1;
	return 0;
}
