/*$Header: /data/petsun4/data1/src_solaris/normalize_4dfp/RCS/normalize_4dfp.c,v 1.9 2011/01/28 03:59:25 avi Exp $*/
/*$Log: normalize_4dfp.c,v $
 * Revision 1.9  2011/01/28  03:59:25  avi
 * histogram output (option -h) now refelcts mode 1000 scaling
 *
 * Revision 1.8  2010/12/24  03:49:59  avi
 * option -m
 *
 * Revision 1.7  2009/02/22  01:57:11  avi
 * option -h
 *
 * Revision 1.6  2007/09/19  03:39:30  avi
 * linux compliant
 *
 * Revision 1.5  2006/09/29  17:07:31  avi
 * Solaris 10
 *
 * Revision 1.4  2005/12/19  02:21:06  avi
 * remove links to libmri
 *
 * Revision 1.3  1999/03/24  00:27:46  avi
 * cleaning
 * new img2hist
 *
 * Revision 1.2  1997/09/17  22:44:47  avi
 * new options -v0 -v1 -v2 (variable intensity normalization volume)
 *
 * Revision 1.1  1997/05/23  02:01:12  yang
 * Initial revision
 *
 * Revision 1.2  1996/07/15  05:48:43  avi
 * -z
 *
 * Revision 1.1  1996/07/10  17:54:19  avi
 * Initial revision
 **/
/*_________________________________________________________________________
Module:		normalize_4dfp.c
Date:		April 12, 1996
Authors:	Bob McKinstry, Avi Snyder, Tom Yang
Description:	Remove time dependent scaling inhomogenieties from a 4dfp 
		stack and write normalized 4dfp stack as "_norm" file.
_________________________________________________________________________*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>				/* getpid () */
#include <math.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>
#include <librms.h>

#define MAXL	256
#define NBIN	105

/*************/
/* externals */
/*************/
extern void	fsquash_	(float *stack, int *nvox, int *nvol, float *imag);				/* fnormalize.f */
extern void	frmo_		(float *imgt,  int *nvox, int *nf_func);					/* fnormalize.f */
extern void	fts_var_	(float *stack, int *np, int *nv, short *mask, float *vol, float *var);		/* fnormalize.f */
extern void	fslice_means_	(float *imag,  int *np, int *nz, short *mask, int *n_z, float *v_z);		/* fnormalize.f */
extern void	fnorm_volume_	(float *stack, int *np, int *nv, short *mask, float *vol, float *amean);	/* fnormalize.f */
extern void	fnorm_slice_	(float *stack, int *np, int *nz, int *nv, short *mask, float *v_z);		/* fnormalize.f */
extern void	ffind_mode_	(int *hist, int *nbin, float *binval, float *amode);				/* fnormalize.f */
extern void	img2hist (float *imgt, int nvox, float *minval, float *maxval, int *hist, int nbin);		/* img2hist.c */

/***********/
/* globals */
/***********/
static char rcsid[] = "$Id: normalize_4dfp.c,v 1.9 2011/01/28 03:59:25 avi Exp $";
static char program[MAXL];

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void write_command_line (FILE *outfp, int argc, char **argv) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
}

int main (int argc, char *argv[]) {
/*************/
/* image i/o */
/*************/
	FILE		*imgfp, *mskfp;			/* 4dfp stacks */
	FILE		*datfp;				/* data file for timeseries rms */
	FILE		*tmpfp;				/* temp process log for rec file */
	char		imgfile[MAXL], imgroot[MAXL];
	char		mskfile[MAXL], mskroot[MAXL] = "";
	char		outfile[MAXL], outroot[MAXL];
	char		trailer[MAXL] = "norm", datfile[MAXL], tmpfile[MAXL];

/********************/
/* image dimensions */
/********************/
	IFH		ifh;
	long		total_dim;
	int		nf_func, nf_anat = 0;
	int		volume_dim, slice_dim;
	int		imgdim[4], orient,  isbig;
	int		mskdim[4], orientm, isbigm;
	float		voxdim[3], voxdimm[3];		/* voxel dimensions in mm */
	float		fhalf = 0.02;			/* imag2mask half frequency in reciprocal mm */
	float		crit = 0.3;			/* mask inclusion criterion */
	char		control = '\0';

/****************/
/* image memory */
/****************/
	float		*v_z;				/* slice means */
	float		*imgt;				/* 4dfp stack */
	float		*imgf, *img1;			/* one frame */
	float		*vart;				/* timeseries variance */
	short		*mask;				/* region of evaulation */
	int		*n_z;				/* pixels per slice in mask */

/*************************/
/* voxel value histogram */
/*************************/
	float		amean, amode, factor;
	float		minval, maxval;
	float		binval[NBIN];			/* histogram bin values */
	int		hist[NBIN];			/* voxel value histogram */
	int		nbin = NBIN;
	int		nbin1 = NBIN - NBIN / 3;

/***********/
/* utility */
/***********/
	int		c, i, j, k;
	char		command[MAXL], *ptr;

/*********/
/* flags */
/*********/
	int		status = 0;
	int		hist_flag = 0;
	int		self_flag = 0;			/* normalize to self */
	int		rmo_flag = 0;			/* remove timeseries mean from each voxel */
	int		norm_type = 1;			/* 1 = by volume (frame); 2 = by slice; 0 = nop */

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 's': self_flag++;			break;
				case 'z': rmo_flag++;			break;
				case 'h': hist_flag++;			break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
				case 'v': norm_type = atoi (ptr);	*ptr = '\0'; break;
				case 'n': nf_anat = atoi (ptr);		*ptr = '\0'; break;
				case 'a': strcpy (trailer, ptr);	*ptr = '\0'; break;
				case 'm': getroot (ptr, mskroot);	*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0:	getroot (argv[i], imgroot); 	k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:\t%s <(4dfp) image>\n", program);
		printf (" e.g.,\t%s -n3 my_run_4dfp\n", program);
		printf ("      \t%s -n3 -v2 my_run_4dfp\n", program);
		printf ("\toption\n");
		printf ("\t-n\tspecify number of pre-functional frames\n");
		printf ("\t-v0\tno           frame to frame intensity stabilization\n");
		printf ("\t-v1\tvolume based frame to frame intensity stabilization (default)\n");
		printf ("\t-v2\tslice  based frame to frame intensity stabilization\n");
		printf ("\t-s\tdisable mode=1000 normalization\n");
		printf ("\t-z\tsubtract mean volume from functional frames\n");
		printf ("\t-h\tcreate <image>.hist file suitable for plotting, e.g., with xmgr\n");
		printf ("\t-a<str>\tspecify trailer (default=\"norm\")\n");
		printf ("\t-m<str>\tread specified 4dfp mask (default blur & threshold input image)\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		exit (1);
	}

	sprintf (datfile, "%s_norm.dat", imgroot);
	sprintf (imgfile, "%s.4dfp.img", imgroot);
/***********************************/
/* get input 4dfp image dimensions */
/***********************************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	if (Getifh (imgfile, &ifh)) errr (program, imgroot);
	if (!control) control = (isbig) ? 'b' : 'l';
	nf_func		= imgdim[3]  - nf_anat;
	slice_dim	= imgdim[0]  * imgdim[1];
	volume_dim	= slice_dim  * imgdim[2];
	total_dim	= volume_dim * imgdim[3];

	if (!(imgt = (float *) malloc (total_dim  * sizeof (float)))
	||  !(img1 = (float *) malloc (volume_dim * sizeof (float)))
	||  !(imgf = (float *) malloc (volume_dim * sizeof (float)))
	||  !(mask = (short *) calloc (volume_dim,  sizeof (short)))
	||  !(vart = (float *) malloc (imgdim[3]  * sizeof (float)))
	||  !(v_z =  (float *) malloc (imgdim[2]  * sizeof (float)))
	||  !(n_z =  (int *)   malloc (imgdim[2]  * sizeof (int)))) errm (program);

	if (!(imgfp = fopen (imgfile, "rb")) || eread (imgt, total_dim, isbig, imgfp)
	|| fclose (imgfp)) errr (program, imgfile);

/********************************************/
/* compute average functional frames volume */
/********************************************/
	fsquash_ (imgt + volume_dim * nf_anat, &volume_dim, &nf_func, imgf);

/******************************/
/* open temp process log file */
/******************************/
	sprintf (tmpfile, "norm%d", getpid ());
	if (!(tmpfp = fopen (tmpfile, "w"))) errw (program, tmpfile);

	if (strlen (mskroot)) {
/*************/
/* read mask */
/*************/
		sprintf (mskfile, "%s.4dfp.img", mskroot);
		if (get_4dfp_dimoe (mskroot, mskdim, voxdimm, &orientm, &isbigm) < 0) errr (program, mskfile);
		status = orientm != orient;
		for (k = 0; k < 3; k++) {
			status |= mskdim[k] - imgdim[k];
			status |= (fabs ((double) (voxdimm[k] - voxdim[k])) > 0.0001);
		}
		if (status) {
			fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot, mskroot);
			remove (tmpfile);
			exit (-1);
		}
		printf ("Reading: %s\n", mskfile);
		if (!(mskfp = fopen (mskfile, "rb"))
		|| eread (img1, volume_dim, isbigm, mskfp)
		|| fclose (mskfp)) errr (program, mskfile);
		for (i = 0; i < volume_dim; i++) if (img1[i] > (float) 1.e-37) mask[i] = 1;
		fprintf (tmpfp, "Mask read from %s\n", mskroot);
	} else {
/***********************************/
/* create mask from squashed image */
/***********************************/
		imag2mask (imgdim + 0, imgdim + 1, imgdim + 2, imgf, mask, voxdim, &fhalf, &crit);
	}
/*********************************/
/* apply mask to averaged volume */
/*********************************/
	for (i = 0; i < volume_dim; i++) if (!mask[i]) imgf[i] = 0.;

/*********************************/
/* compute voxel value histogram */
/*********************************/
	if (self_flag) {
		fprintf (tmpfp, "4dfp stack normalized to self\n");
	} else {
		img2hist (imgf, volume_dim, &minval, &maxval, hist, nbin);
		for (k = 0; k < NBIN; k++) binval[k] =  minval + k * (maxval - minval) / NBIN;
		ffind_mode_ (hist + NBIN / 3, &nbin1, binval + NBIN / 3, &amode);
		fprintf (tmpfp,  "4dfp stack original      mode intensity= %8.2f\n", amode);
		fprintf (tmpfp,  "4dfp stack normalized to mode intensity= %8.2f\n", 1000.);
		fprintf (stdout, "4dfp stack original      mode intensity= %8.2f\n", amode);
		fprintf (stdout, "4dfp stack normalized to mode intensity= %8.2f\n", 1000.);
		factor = 1000. / amode;
		for (k = 0; k < total_dim;  k++) imgt[k] *= factor;
		for (k = 0; k < volume_dim; k++) imgf[k] *= factor;
	}

	if (hist_flag) {
/************************************************************/
/* remake histogram after mode 1000 scaling on range 0-1470 */
/************************************************************/
		for (j = 0; j < nbin; j++) hist[j] = 0.;
		maxval = 1470.;		/* 14*nbin */
		for (i = 0; i < volume_dim; i++) {
			j = (int) (imgf[i]*nbin)/maxval;
			if (j > 10 && j < nbin) hist[j]++;
		}
		sprintf (outfile, "%s.hist", imgroot);
		if (!(datfp = fopen (outfile, "w"))) errw (program, outfile);
		printf ("Writing: %s\n", outfile);
		write_command_line (datfp, argc, argv);
		fprintf (datfp,  "#4dfp original mode intensity= %8.2f\n", amode);
		fprintf (datfp, "%10.2f%10d\n", 0., 0);
		for (j = 0; j < nbin; j++) {
			fprintf (datfp, "%14.6f %10d\n", (j+0)*maxval/nbin, hist[j]);
			fprintf (datfp, "%14.6f %10d\n", (j+1)*maxval/nbin, hist[j]);
		}
		fprintf (datfp, "%10.2f%10d\n", maxval, 0);
		fclose (datfp);
	}

	switch (norm_type) {
	case 2:
		fslice_means_ (imgf, &slice_dim, imgdim + 2, mask, n_z, v_z);
		fnorm_slice_  (imgt, &slice_dim, imgdim + 2, imgdim + 3, mask, v_z);
		for (k = 0; k < imgdim[2]; k++) fprintf (tmpfp, "Masked slice %2d mean=%10.4f\n", k + 1, v_z[k]);
		break;
	case 1:
		fnorm_volume_ (imgt, &volume_dim, imgdim + 3, mask, imgf, &amean);
		fprintf (tmpfp, "Masked image mean=%10.4f\n", amean);
		break;
	default:
		fprintf (tmpfp, "No frame to frame intensity stabilization\n");
		break;
	}

/*******************************/
/* compute timeseries variance */
/*******************************/
	if (imgdim[3] > 1) {
		fts_var_ (imgt, &volume_dim, imgdim + 3, mask, imgf, vart);
		fprintf (stdout, "Writing: %s\n", datfile);
		if (!(datfp = fopen (datfile, "w"))) errw (program, datfile);
		for (i = 0; i < nf_anat; i++)		fprintf (datfp, "#%9d%10.4f\n", i + 1, sqrt (vart[i]));
		for (i = nf_anat; i < imgdim[3]; i++)	fprintf (datfp, "%10d%10.4f\n", i + 1, sqrt (vart[i]));
		fclose (datfp);
	}

	if (rmo_flag) {
		printf ("Removing voxel value offsets from timeseries\n");
		frmo_ (imgt + volume_dim * nf_anat, &volume_dim, &nf_func);
		fprintf (tmpfp, "Functional frames voxel value offset removed\n");
	}
	fclose (tmpfp);

	sprintf (outroot, "%s_%s", imgroot, trailer);
	sprintf (outfile, "%s.4dfp.img", outroot);
/*******************************/
/* write normalized 4dfp stack */
/*******************************/
	fprintf (stdout, "Writing: %s\n", outfile);
	if (!(imgfp = fopen (outfile, "wb")) || ewrite (imgt, total_dim, control, imgfp)
	|| fclose (imgfp)) errw (program, imgfile);

/*******/
/* rec */
/*******/
	startrece (outfile, argc, argv, rcsid, control);
	catrec    (tmpfile); remove (tmpfile);
	catrec    (imgfile);
	endrec    ();

/*******/
/* ifh */
/*******/
	if (writeifhe (program, outfile, imgdim, voxdim, orient, control)) errw (program, outfile);
	
/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s -r%d", outroot, 2000);
	printf ("%s\n", command); status = system (command);
	
	free (imgt); free (imgf); free (img1); free (mask); free (vart);
	free (v_z); free (n_z);
	exit (0);
}
