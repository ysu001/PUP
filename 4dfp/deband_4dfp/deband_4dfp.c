/*$Id: deband_4dfp.c,v 1.13 2009/02/22 02:06:48 avi Exp $*/
/*$Log: deband_4dfp.c,v $
 * Revision 1.13  2009/02/22  02:06:48  avi
 * mask short -> int
 *
 * Revision 1.12  2006/09/27  22:18:01  avi
 * Solaris 10
 *
 * Revision 1.11  2004/10/12  20:10:43  rsachs
 * Removed 'writeifh' from this module. Replaced it with the one that's in ../nonlin.
 *
 * Revision 1.10  2004/10/11  22:02:54  rsachs
 * Installed 'errm','errr','errw','getroot','get_4dfp_dimo'. Removed 'Get4dfpDimN'.
 *
 * Revision 1.9  2001/06/29  18:52:03  avi
 * MAXF -> 2048
 *
 * Revision 1.8  1999/11/20  00:55:49  avi
 * print input image orientation as a check on Getifh ()
 *
 * Revision 1.7  1999/03/05  05:42:07  avi
 * operate on one frame at a time
 * accept run formats (-F option)
 *
 * Revision 1.6  1998/06/28  01:13:24  avi
 * Revision 1.5  1998/06/20  01:17:29  avi
 *  debandg.f rec.c
 *
 * Revision 1.4  1997/05/23  00:54:37  yang
 * new rec macros
 *
 * Revision 1.3  1997/04/28  20:56:37  yang
 * Working Solaris version.
 *
 * Revision 1.2  1997/04/28  00:43:47  yang
 * Working Solaris version.
 *
 * Revision 1.3  1996/07/10  06:09:30  avi
 * rec file saving
 * Revision 1.2  1996/05/17  23:43:28  avi
 * nf_anat logic
 * Revision 1.1  1996/05/17  21:47:37  avi
 * Initial revision
 **/
/*___________________________________________________________________________________
Module:		deband_4dfp.c
Date:		April 12, 1996
Authors:	Avi Snyder with Tom Yang support
Description:	Remove slice dependent scaling inhomogenieties
		from a 4dfp stack and write debanded 4dfp stack to "_dbnd" file.
____________________________________________________________________________________*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>
#include <librms.h>

#define MAXL 	256
#define MAXF	2048
#define MAXS 	256

/*************/
/* externals */
/*************/
extern void	fband_     (float *imag, int *np, int *nz, int *mask, int *n_z, float *v_z);		/* fdeband.f */
extern void	fbandr_    (float *imag, int *np, int *nz, int *mask, float *alpha);			/* fdeband.f */
extern void	fdebandr_  (float *imag, int *np, int *nz, int *nv, float *alpha);			/* fdeband.f */
extern void	fdebandg_  (float *imag, int *np, int *nz, int *nv, float *beta);			/* debandg.f */
extern void	get_flipz_ (float *imag, int *np, int *nz, int *n_z, float *v_z, int *flipz);		/* debandg.f */
extern void	fbandg_    (float *imag, int *np, int *nz, int *mask, float *beta);			/* debandg.f */
extern void	fbande_   (int *flipz, float *imag, int *np, int *nz, int *mask, float *gamma);		/* debandg.f */
extern void	fdebande_ (int *flipz, float *imag, int *np, int *nz, int *nv, float *gamma);		/* debandg.f */
extern int	expandf (char *string, int len);							/* expandf.c */

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[] = "$Id: deband_4dfp.c,v 1.13 2009/02/22 02:06:48 avi Exp $";
int main (int argc, char *argv[]) {
	FILE		*outfp, *datfp;			/* output */
	FILE		*imgfp;				/* input */
	char		imgroot[MAXL], imgfile[MAXL];	/* input 4dfp file name */
	char		outroot[MAXL], outfile[MAXL];	/* debanded 4dfp file name */
	char		datfile[MAXL];			/* debanding parameter listing */


/**************/
/* brain mask */
/**************/
	float		fhalf = 0.02;			/* img2lmask half frequency in reciprocal mm */
	float		crit = 0.3;			/* mask inclusion criterion */

/***************/
/* 4dfp arrays */
/***************/
	int		*mask;				/* region of evaulation */
	float		*imgt;				/* frame buffer */
	float		*imgv, *imgd;			/* stack averages before & after debanding */
	float		voxdim[3];			/* voxel dimensions in mm */
	int		imgdim[4], slice_dim, vdim, zdim, orient, isbig, osbig;
	int		nf_anat = 0, nf_func;
	char		control = '\0';			/* output endian */
	
/***********************/
/* debanding algorithm */
/***********************/
	char		format[MAXF] = "";
	float		alpha, alphaf;			/* debanding parameter */
	float		v_z[MAXS], v_z1[MAXS];		/* before and after slice means */
	int		n_z[MAXS];			/* pixels per slice in mask */

/***********/
/* utility */
/***********/
	int		c, i, j, k;
	int		one = 1;
	char		*ptr, command[MAXL], program[MAXL];

/*********/
/* flags */
/*********/
	int		status;
	int		flipz;
	int		method = 0;	/* 0 = flat even/odd; 1 = linear gradient; 2 = exponential gradient */

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'e': method = 2;		break;
				case 'g': method = 1;		break;
				case 'n': nf_anat = atoi (ptr);	*ptr = '\0'; break;
				case 'F': strcpy (format, ptr);	*ptr = '\0'; break;
				case '@': control = *ptr++;	*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	getroot (argv[i], imgroot); k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:	%s <(4dfp) image>\n", program);
		printf (" e.g.:	%s -n3 mybold\n", program);
		printf (" e.g.:	%s -F\"3x125+\" mybold\n", program);
		printf ("\toption\n");
		printf ("\t-e\tdeband by exponential gradient model (default flat model)\n");
		printf ("\t-g\tdeband by linear gradient model (default flat model)\n");
		printf ("\t-n<int>\tspecify number of pre-functional frames\n");
		printf ("\t-F<str>\tspecify complete functional/non-functional format\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		exit (1);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
	sprintf (outroot, "%s_dbnd",     imgroot);
	sprintf (datfile, "%s.dat",      outroot);
	sprintf (outfile, "%s.4dfp.img", outroot);

/***********************************/
/* get input 4dfp image dimensions */
/***********************************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig)) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';
	slice_dim	= imgdim[0]*imgdim[1];
	zdim		= imgdim[2];
	vdim		= slice_dim*zdim;

/*****************/
/* expand format */
/*****************/
	if (!strlen (format)) {
		nf_func = imgdim[3] - nf_anat;
		sprintf (format, "%dx%d+", nf_anat, nf_func);
	}
	if (expandf (format, MAXF)) exit (-1);
	for (nf_anat = i =  0;  i < imgdim[3]; i++) if (format[i] == 'x') nf_anat++;
	if ((nf_func = imgdim[3] - nf_anat) <= 0) {
		fprintf (stderr, "%s: no functional frames\n", imgfile);
		exit (-1);
	}

	if (!(imgt = (float *) calloc (vdim, sizeof (float)))
	||  !(imgd = (float *) calloc (vdim, sizeof (float)))
	||  !(imgv = (float *) calloc (vdim, sizeof (float)))
	||  !(mask = (int *) calloc (vdim, sizeof (int)))) errm (program);

	fprintf (stdout, "Reading: %s\n", imgfile);
	if (!(imgfp = fopen (imgfile, "rb")))  errr (program, imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "w+b"))) errw (program, outfile);

/*******************************************/
/* compute average functional frame volume */
/*******************************************/
	for (i = 0; i < imgdim[3]; i++) {
		if (eread (imgt, vdim, isbig, imgfp)) errr (program, imgfile);
		if (format[i] != 'x') for (k = 0; k < vdim; k++) imgv[k] += imgt[k];
	}
	for (k = 0; k < vdim; k++) imgv[k] /= nf_func;

/***************/
/* create mask */
/***************/
	img2lmask (imgdim + 0, imgdim + 1, imgdim + 2, imgv, mask, voxdim, &fhalf, &crit);

/*********************************/
/* apply mask to averaged volume */
/*********************************/
	for (i = 0; i < vdim; i++) if (!mask[i]) imgv[i] = 0.0;

/***********************/
/* compute slice means */
/***********************/
	fband_ (imgv, &slice_dim, &zdim, mask, n_z, v_z);

/************************************************************/
/* compute debanding parameter and deband functional frames */
/************************************************************/
	get_flipz_ (imgv, &slice_dim, &zdim, n_z, v_z, &flipz);
	switch (method) {
		case 0:	fbandr_ (imgv, &slice_dim, &zdim, mask, &alphaf);		break;
		case 1:	fbandg_ (imgv, &slice_dim, &zdim, mask, &alphaf);		break;
		case 2:	fbande_ (&flipz, imgv, &slice_dim, &zdim, mask, &alphaf);	break;
	}
	for (i = 0; i < imgdim[3]; i++) {
		if (format[i] == 'x') continue;
		printf ("Frame %3d ", i + 1);
		if (fseek (imgfp, (long) i * vdim * sizeof (float), SEEK_SET)
		||  eread (imgt, vdim, isbig, imgfp)) errr (program, imgfile);
		switch (method) {
		case 0:	fdebandr_ (imgt, &slice_dim, &zdim, &one, &alphaf);		break;
		case 1:	fdebandg_ (imgt, &slice_dim, &zdim, &one, &alphaf);		break;
		case 2:	fdebande_ (&flipz, imgt, &slice_dim, &zdim, &one, &alphaf);	break;
		}
		if (fseek (outfp, (long) i * vdim * sizeof (float), SEEK_SET)
		|| ewrite (imgt, vdim, control, outfp)) errw (program, outfile);
	}

/***************************************/
/* recompute average functional volume */
/***************************************/
/*	osbig = (CPU_is_bigendian ()) ? !(control == 'l' || control == 'L') : (control == 'b' || control == 'B');	*/
	switch (control) {
		case 'b': case 'B': osbig = 1; break;
		case 'l': case 'L': osbig = 0; break;
		default: osbig = CPU_is_bigendian();
	}
	for (i = 0; i < imgdim[3]; i++) {
		if (format[i] == 'x') continue;
		if (fseek (outfp, (long) i * vdim * sizeof (float), SEEK_SET)
		||  eread (imgt, vdim, osbig, outfp)) errw (program, outfile);
		for (k = 0; k < vdim; k++) imgd[k] += imgt[k];
	}
	for (k = 0; k < vdim; k++) imgd[k] /= nf_func;

/*********************************/
/* compute and write slice means */
/*********************************/
	fband_ (imgd, &slice_dim, &zdim, mask, n_z, v_z1);
	datfp = fopen (datfile, "w");
	for (i = 0; i < zdim; i++) if (n_z[i]) fprintf (datfp, "%10d%10.4f%10.4f%10d\n", i + 1, v_z[i], v_z1[i], n_z[i]);
	fclose (datfp);

/*******/
/* rec */
/*******/
	startrece (outfile, argc, argv, rcsid, control);

/*********************************************/
/* deband non-functional frames individually */
/*********************************************/
	for (i = 0; i < imgdim[3]; i++) {
		if (format[i] != 'x') continue;
		printf ("Frame %3d ", i + 1);
		if (fseek (imgfp, (long) i * vdim * sizeof (float), SEEK_SET)
		||  gread ((char *) imgt, sizeof (float), vdim, imgfp, isbig)) errr (program, imgfile);
		switch (method) {
		case 0:	fbandr_   (imgt, &slice_dim, &zdim, mask, &alpha);
			fdebandr_ (imgt, &slice_dim, &zdim, &one, &alpha);
			break;
		case 1:	fbandg_   (imgt, &slice_dim, &zdim, mask, &alpha);
			fdebandg_ (imgt, &slice_dim, &zdim, &one, &alpha);
			break;
		case 2:	fbande_   (&flipz, imgt, &slice_dim, &zdim, mask, &alpha);
			fdebande_ (&flipz, imgt, &slice_dim, &zdim, &one, &alpha);
			break;
		}
		if (fseek (outfp, (long) i * vdim * sizeof (float), SEEK_SET)
		|| ewrite (imgt, vdim, control, outfp)) errw (program, outfile);
		sprintf (command, "Frame %10d slice multipliers: even=%8.6f odd=%8.6f\n",
			i + 1, 1.0 + alpha, 1.0 - alpha);
		printrec (command);
	}
	sprintf  (command, "Functional frame slice multipliers: even=%8.6f odd=%8.6f\n", 1.0 + alphaf, 1.0 - alphaf);
	printrec (command);
	catrec   (imgfile);
	endrec   ();

/***************/
/* ifh and hdr */
/***************/
	if (writeifhe (program, outroot, imgdim, voxdim, orient, control)) errw (program, outroot);
	sprintf (command, "ifh2hdr -r2000 %s", outroot); status = system (command);

	free (imgt); free (imgd); free (imgv); free (mask);
	exit (status);
}

