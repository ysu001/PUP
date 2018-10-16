/*$Header: /data/petsun4/data1/src_solaris/cross_realign3d_4dfp/RCS/mat2dat.c,v 1.10 2010/02/10 06:59:09 avi Exp avi $*/
/*$Log: mat2dat.c,v $
 * Revision 1.10  2010/02/10  06:59:09  avi
 * options -f and -L
 * differentiated output (*.ddat) slightly altered
 * command line included in output
 *
 * Revision 1.9  2007/08/08  02:40:15  avi
 * gcc compliant
 *
 * Revision 1.8  2006/09/28  21:59:18  avi
 * Solaris 10
 *
 * Revision 1.7  2004/08/02  21:31:12  avi
 * -R option (save per-run quantities)
 * -D option (save differentiated)
 *
 * Revision 1.6  2001/06/29  21:34:04  avi
 * MAXF -> 2048 (with new limit check)
 *
 * Revision 1.5  2000/04/20  03:38:05  avi
 * correct usage
 *
 * Revision 1.4  1999/10/06  21:20:13  avi
 * RADIUS -> 50.
 *
 * Revision 1.3  1999/10/06  05:25:21  avi
 * compute trajectory statistics
 *
 * Revision 1.2  1998/07/17  03:41:30  avi
 * fix Usage
 *
 * Revision 1.1  1998/07/17  03:37:55  avi
 * Initial revision
 **/
/*_________________________________________________________________________
Module:		mat2dat.c
Date:		03-Aug-96
Authors:	Avi Snyder
Description:	Convert cross_realign3d_4dfp mat file to dat.
_________________________________________________________________________*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <librms.h>		/* warp2param() in frmsmri.f */

#define MAXL	256
#define NF_ANAT	3
#define MAXF	2048
#define RADIUS	50.

/*************/
/* externals */
/*************/
int	expandf (char *string, int len);		/* expandf.c */

/***********/
/* globals */
/***********/
static char rcsid[] = "$Id: mat2dat.c,v 1.10 2010/02/10 06:59:09 avi Exp avi $";
static char program[MAXL];

void errr (char* program, char* filespc) {
	fprintf (stderr, "%s: %s read error\n", program, filespc);
	exit (-1);
}

void errw (char* program, char* filespc) {
	fprintf (stderr, "%s: %s write error\n", program, filespc);
	exit (-1);
}

void errF (char* program, char* filespc) {
	fprintf (stderr, "%s: frame count limit (%d) exceeded in %s\n", program, MAXF, filespc);
	exit (-1);
}

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n#%s\n", rcsid);
}

int main (int argc, char *argv[]) {
	FILE	*matfp;			/* mat file */
	FILE	*datfp;			/* dat file */
	char	matroot[MAXL], matfile[MAXL], datfile[MAXL];
	float	dpr;			/* degrees per radian */
	float	s4[4];			/* intensity scaling factors */
	float	t4[16];			/* 4 x 4 transform */
	float	param[12];
	int	mode = 3;		/* assume 3D no stretch */

/***********************************/
/* trajectory variance computation */
/***********************************/
	char	format[MAXF + 1] = "";		/* pattern of frames to count */
	int	frame_count = 0, nframe, nframed;
	int	nf_anat = NF_ANAT;	/* prefunctional frames */
	float	u1[6], u2[6], du1[6], du2[6];
	float	s1 = 0.0, s2 = 0.0, ds1 = 0.0, ds2 = 0.0;
        float	trajectory[MAXF][6], dtraj[MAXF][6];
        float	scale[MAXF], dscale[MAXF];
	float	radius = RADIUS;

/***********/
/* utility */
/***********/
	char	*ptr, command[MAXL];
	int	c, i, j, k;
	float	t;

/*********/
/* flags */
/*********/
	int	debug = 0;
	int	save_dif = 0;
	int	save_rel = 0;
	int	write_cwd = 0;

	fprintf (stdout, "%s\n", rcsid);
	dpr = 45. / atan (1.);
	if (ptr = strrchr (argv[0], '/')) ptr++; else ptr = argv[0];
	strcpy (program, ptr);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'L': write_cwd++;				break;
				case 'R': save_rel++;				break;
				case 'D': save_dif++;				break;
				case 'n': nf_anat = atoi (ptr);			*ptr = '\0'; break;
				case 'r': radius = atof (ptr);			*ptr = '\0'; break;
				case 'f': strncpy (format, ptr, MAXF);		*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0:	strcpy (matroot, argv[i]); k++; break;
		}	
	}
	if (k < 1) {
		fprintf (stderr, "Usage:\t%s <mat_file>\n", program);
		fprintf (stderr, "e.g.,\t%s atten5_b1_rms4_dbnd_xr3d[.mat]\n", program);
		fprintf (stderr, "\toption\n");
		fprintf (stderr, "\t-R\tsave trajectory relative to run mean (remove accumulated movememnt)\n");
		fprintf (stderr, "\t-D\tsave differentiated trajectory\n");
		fprintf (stderr, "\t-L\twrite local (in $cwd) (default write parallel to <mat_file>)\n");
		fprintf (stderr, "\t-n<int>\tspecify number of pre steady state frames (default=%d)\n", NF_ANAT);
		fprintf (stderr, "\t-r<flt>\tspecify head radius in mm for total motion computation (default=%.0fmm)\n", RADIUS);
		fprintf (stderr, "\t-f<str>\tspecify frames to count format, e.g., \"4x120+4x76+\"\n");
		fprintf (stderr, "N.B.:\t-f option overrides -n\n");
		exit (1);
	}

/***********************************************/
/* parse out mat filename extension if present */
/***********************************************/
	while (ptr = strrchr (matroot, '.')) {
		if (!strcmp (ptr, ".mat")) *ptr = '\0';
		else break;
	}
	sprintf (matfile, "%s.mat", matroot);

	if (!(matfp = fopen (matfile, "r"))) errr (program, matfile);
	printf ("Reading: %s\n", matfile);
	while (fscanf (matfp, "%*s %*s %d", &i) != EOF) {
		if (frame_count >= MAXF) errF (program, matfile);
		if (debug) printf ("t4 frame %d\n", i);
		for (k = 0; k < 4; k++) {
			fscanf (matfp, "%f %f %f %f", t4 + 0 + k, t4 + 4 + k, t4 + 8 + k, t4 + 12 + k);
		}
		if (debug) for (k = 0; k < 4; k++) {
			printf ("%10.6f%10.6f%10.6f%10.4f\n", t4[k + 0], t4[k + 4], t4[k + 8], t4[k + 12]);
		}
		warp2param_ (&mode, t4, param);
		fscanf (matfp, "%*s %*s %d", &i);
		if (debug) printf ("s4 frame %d\n", i);
		fscanf (matfp, "%f %f %f %f", s4 + 0, s4 + 1, s4 + 2, s4 + 3);
		if (debug) printf ("%10.6f%10.6f%10.6f%10.4f\n", s4[0], s4[1], s4[2], s4[3]);
		for (k = 0; k < 6; k++) trajectory[frame_count][k] = param[k];
		scale[frame_count] = s4[0];
		frame_count++;
	}
	fclose (matfp);

/*********************************/
/* set up frames-to-count format */
/*********************************/
	if (strlen (format)) {
		if (k = expandf (format, MAXF)) exit (k);
		if ((k = strlen (format)) != frame_count) {
			fprintf (stderr, "format codes frame count (%d) not equal to data length (%d)\n", k, frame_count);
			exit (-1);
		}
	} else {
		for (i = 0; i < nf_anat; i++) format[i] = 'x';
		for (i = nf_anat + 1; i < frame_count; i++) format[i] = '+';
	}

/*********************************/
/* compute trajectory statistics */
/*********************************/
	for (i = 0; i < frame_count; i++) {
		for (k = 0; k < 6; k++) {
			dtraj[i][k] = (i) ? trajectory[i][k] - trajectory[i-1][k] : 0.0;
		}
		dscale[i] = (i) ? scale[i] - scale[i-1] : 0.0;
	}

	for (k = 0; k < 6; k++) du1[k] = du2[k] = u1[k] = u2[k] = 0.0;
	for (nframe = i = 0; i < frame_count; i++) if (format[i] != 'x') {
		for (k = 0; k < 6; k++) {
			u1[k] += trajectory[i][k];
			u2[k] += trajectory[i][k]*trajectory[i][k];
		}
		s1 += scale[i];
		s2 += scale[i]*scale[i];
		nframe++;
	}
	for (nframed = 0, i = 1; i < frame_count; i++) if (format[i-1] != 'x' && format[i] != 'x') {
		for (k = 0; k < 6; k++) {
			du1[k] += dtraj[i][k];
			du2[k] += dtraj[i][k]*dtraj[i][k];
		}
		ds1 += dscale[i];
		ds2 += dscale[i]*dscale[i];
		nframed++;
	}

	for (k = 0; k < 6; k++) {
		u1[k]  /= nframe;
		u2[k]  -= nframe * u1[k]*u1[k];
		u2[k]  /= nframe - 1;
	}
	s1 /= nframe;
	s2 -= nframe * s1*s1;
	s2 /= nframe - 1;

	for (k = 0; k < 6; k++) {
		du1[k] /= nframed;
		du2[k] -= nframed * du1[k]*du1[k];
		du2[k] /= nframed - 1;
	}
	ds1 /= nframed;
	ds2 -= nframed * ds1*ds1;
	ds2 /= nframed - 1;

/*******************/
/* write dat files */
/*******************/
	if ((ptr = strrchr (matroot, '/')) && write_cwd) ptr++; else ptr = matroot;
	sprintf (datfile, "%s.dat", ptr);
	printf ("Writing: %s\n", datfile);
	if (!(datfp = fopen (datfile, "w"))) errw (program, datfile);
	write_command_line (datfp, argc, argv);
	fprintf (datfp, "#frame    dx(mm)    dy(mm)    dz(mm)    X(deg)    Y(deg)    Z(deg)     scale\n");
	for (i = 0; i < frame_count; i++) {
		fprintf (datfp, "%6d", i + 1);
		for (k = 0; k < 3; k++) fprintf (datfp, "%10.4f", trajectory[i][k]);
		for (k = 3; k < 6; k++) fprintf (datfp, "%10.4f", dpr * trajectory[i][k]);
		fprintf (datfp, "%10.4f\n", scale[i]);
	}
	fprintf (datfp, "#counting %d out of %d frames\n", nframe, frame_count);
	fprintf (datfp, "#%-5s", "mean");
	for (k = 0; k < 3; k++) fprintf (datfp, "%10.4f", u1[k]);
	for (k = 3; k < 6; k++) fprintf (datfp, "%10.4f", dpr * u1[k]);
	fprintf (datfp, "%10.4f\n", s1);
	fprintf (datfp, "#%-5s", "s.d.");
	for (k = 0; k < 3; k++) fprintf (datfp, "%10.4f", sqrt (u2[k]));
	for (k = 3; k < 6; k++) fprintf (datfp, "%10.4f", dpr * sqrt (u2[k]));
	fprintf (datfp, "%10.4f\n", sqrt (s2));
	fprintf (datfp, "#rms translation (mm)         %.4f\n",       sqrt (u2[0] + u2[1] + u2[2]));
	fprintf (datfp, "#rms rotation (deg)           %.4f\n", dpr * sqrt (u2[3] + u2[4] + u2[5]));
	t = u2[0] + u2[1] + u2[2] + radius*radius * (u2[3] + u2[4] + u2[5]);
	fprintf (datfp, "#total rms movement at radius=%.0f mm     %.4f\n", radius, sqrt (t));
	fclose (datfp);
        sprintf (command, "cat %s", datfile); system (command);

if (save_rel) {
	if ((ptr = strrchr (matroot, '/')) && write_cwd) ptr++; else ptr = matroot;
	sprintf (datfile, "%s.rdat", ptr);
	printf ("Writing: %s\n", datfile);
	if (!(datfp = fopen (datfile, "w"))) errw (program, datfile);
	write_command_line (datfp, argc, argv);
	fprintf (datfp, "#frame    dx(mm)    dy(mm)    dz(mm)    X(deg)    Y(deg)    Z(deg)     scale\n");
	for (i = 0; i < frame_count; i++) {
		fprintf (datfp, "%6d", i + 1);
		for (k = 0; k < 3; k++) fprintf (datfp, "%10.4f", trajectory[i][k] - u1[k]);
		for (k = 3; k < 6; k++) fprintf (datfp, "%10.4f", dpr * (trajectory[i][k] - u1[k]));
		fprintf (datfp, "%10.4f\n", scale[i]/s1);
	}
	fprintf (datfp, "#counting %d out of %d frames\n", nframe, frame_count);
	fprintf (datfp, "#%-5s", "mean");
	for (k = 0; k < 3; k++) fprintf (datfp, "%10.4f", u1[k]);
	for (k = 3; k < 6; k++) fprintf (datfp, "%10.4f", dpr * u1[k]);
	fprintf (datfp, "%10.4f\n", s1);
	fprintf (datfp, "#%-5s", "s.d.");
	for (k = 0; k < 3; k++) fprintf (datfp, "%10.4f", sqrt (u2[k]));
	for (k = 3; k < 6; k++) fprintf (datfp, "%10.4f", dpr * sqrt (u2[k]));
	fprintf (datfp, "%10.4f\n", sqrt (s2)/s1);
	fprintf (datfp, "#rms translation (mm)         %.4f\n",       sqrt (u2[0] + u2[1] + u2[2]));
	fprintf (datfp, "#rms rotation (deg)           %.4f\n", dpr * sqrt (u2[3] + u2[4] + u2[5]));
	t = u2[0] + u2[1] + u2[2] + radius*radius * (u2[3] + u2[4] + u2[5]);
	fprintf (datfp, "#total rms movement at radius=%.0f mm     %.4f\n", radius, sqrt (t));
	fclose (datfp);
}

if (save_dif) {
	if ((ptr = strrchr (matroot, '/')) && write_cwd) ptr++; else ptr = matroot;
	sprintf (datfile, "%s.ddat", ptr);
	printf ("Writing: %s\n", datfile);
	if (!(datfp = fopen (datfile, "w"))) errw (program, datfile);
	write_command_line (datfp, argc, argv);
	fprintf (datfp, "#frame   ddx(mm)   ddy(mm)   ddz(mm)   dX(deg)   dY(deg)   dZ(deg)100*dscale\n");
	for (i = 0; i < frame_count; i++) {
		fprintf (datfp, "%6d", i + 1);
		for (k = 0; k < 3; k++) fprintf (datfp, "%10.4f", dtraj[i][k]);
		for (k = 3; k < 6; k++) fprintf (datfp, "%10.4f", dpr * dtraj[i][k]);
		fprintf (datfp, "%10.4f\n", 100*dscale[i]/s1);
	}
	fprintf (datfp, "#counting %d out of %d frames\n", nframed, frame_count);
	fprintf (datfp, "#%-5s", "s.d.");
	for (k = 0; k < 3; k++) fprintf (datfp, "%10.4f", sqrt (du2[k]));
	for (k = 3; k < 6; k++) fprintf (datfp, "%10.4f", dpr * sqrt (du2[k]));
	fprintf (datfp, "%10.4f\n", 100*sqrt (ds2)/s1);
	fprintf (datfp, "#rms dtranslation (mm/frame)  %.4f\n",       sqrt (du2[0] + du2[1] + du2[2]));
	fprintf (datfp, "#rms drotation (deg/frame)    %.4f\n", dpr * sqrt (du2[3] + du2[4] + du2[5]));
	t = du2[0] + du2[1] + du2[2] + radius*radius * (du2[3] + du2[4] + du2[5]);
	fprintf (datfp, "#total rms velocity at radius=%.0f mm     %.4f mm/frame\n", radius, sqrt (t));
	fclose (datfp);
}

	exit (0);
}
