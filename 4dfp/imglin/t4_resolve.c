/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4_resolve.c,v 1.5 2009/09/09 05:53:14 avi Exp $*/
/*$Log: t4_resolve.c,v $
 * Revision 1.5  2009/09/09  05:53:14  avi
 * Solaris10/linux compliant
 *
 * Revision 1.4  2002/03/22  21:02:25  avi
 * more stdout messages
 *
 * Revision 1.3  1999/10/18  04:58:23  avi
 * iterate on fabs (err - err0)
 *
 * Revision 1.2  1999/10/14  03:43:01  avi
 * write resolved t4 sub and mat files
 *
 * Revision 1.1  1999/10/10  07:20:28  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <t4_sub.h>

#define	MAXF	32
#define MAXL	256
#define RADIUS	50.

/*****************/
/* ft4_resolve.f */
/*****************/
extern void	symerr_ (int *n, float *ck6);
extern void	resolve_rigid_ (float *radius, int *n, float *jac, float *err, float *domega, float *errsum);
extern void	resolve_scale_ (int *n, float *err, float *dscale);

/***********/
/* CKsub.f */
/***********/
extern void	ck6_to_t4_ (float *ed, float *t4);
extern void	ck6_to_param6_ (float *ed, float *param6);
extern void	t4_to_ck6_ (float *t4, float *ed);

/************/
/* getjxy.f */
/************/
extern void	getjxy_ (float *txa, float *tya, float *txy, float *that, float *terr, float *j_xy_xa, float *j_xy_ya);

void errm (char *program) {
	fprintf (stderr, "%s: memory allocation error\n", program);
	exit (-1);
}

void errr (char *program, char *file) {
	fprintf (stderr, "%s: %s read error\n", program, file);
	exit (-1);
}

void errw (char *program, char *file) {
	fprintf (stderr, "%s: %s write error\n", program, file);
	exit (-1);
}

void printT (float* T, float* s) {
	printf ("%10.6f%10.6f%10.6f%10.4f%10.4f%10.4f%10.6f\n", T[0], T[1], T[2], T[3], T[4], T[5], s[0]);
}

void usage (char *program) {
	fprintf (stderr, "Usage:\t%s <image1> <image2> ...\n", program);
	fprintf (stderr, "\toption\n");
	fprintf (stderr, "\t-v	verbose mode\n");
	fprintf (stderr, "\t-m	generate mat file output\n");
	fprintf (stderr, "\t-s	include intensity scale factor in t4 file output\n");
	fprintf (stderr, "\t-w	weight inversely in proportion to scale in sub file output (sum counts mode)\n");
	fprintf (stderr, "\t-o<str>	write resolved output with specified fileroot\n");
	fprintf (stderr, "\t-r<flt>	set VOI rms radius in mm (default=%.0f)\n", RADIUS);
	fprintf (stderr, "N.B.:	%s looks for t4 files <image1>_to_<image2>_t4, <image1>_to_<image3>_t4, ...\n", program);
	fprintf (stderr, "N.B.:	%s automatically strips filename extensions when constructing t4 filenames\n", program);
	exit (1);
}

static char rcsid[] = "$Id: t4_resolve.c,v 1.5 2009/09/09 05:53:14 avi Exp $";
int main (int argc, char *argv[]) {
	FILE 		*fp_t4, *fpsub;
        char		scans[MAXF][MAXL];

	char		*str, command[MAXL], program[MAXL];
	char		outroot[MAXL] = "", subfile[MAXL], t4file[MAXL];
	float		t4[16], scale, Serr[MAXF*(MAXF-1)], dS[MAXF*(MAXF-1)];
	float		Sbar[MAXF*(MAXF-1)], Shat[MAXF*(MAXF-1)];
	float		param6[6];
	float		Tzer[6] = {0., 0., 0., 0., 0., 0.};
	float		Tbar[MAXF*(MAXF-1)][6], T[6], *fptr;
	float		That[MAXF*(MAXF-1)][6], Terr[MAXF*(MAXF-1)][6], dT[MAXF*(MAXF-1)][6];
	float		J_xy_xa[36], J_xy_ya[36];
	float		*J, *Jptr;
	float		radius = RADIUS, degperad;
	float		err0, err = 0.0;
	int		niter = 0;
	int		m, n = 0;
	int		c, i, j, k, l;

/*********/
/* flags */
/*********/
	int		mat_flag = 0;
	int		scale_flag = 0;
	int		weight_flag = 0;
	int		status = 0;
	int		verbose = 0;
	int		zero_test = 0;
	int		debug = 0;

	printf ("%s\n", rcsid);
	if (!(str = strrchr (argv[0], '/'))) str = argv[0]; else str++;
	strcpy (program, str);
	degperad = 45. / atan (1.);

/************************/
/* process command line */
/************************/
	for (n = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case 'm': mat_flag++;		break;
				case 's': scale_flag++;		break;
				case 'v': verbose++;		break;
				case 'w': weight_flag++;	break;
				case 'r': radius = atof (str);		*str = '\0'; break;
				case 'o': strcpy (outroot, str);	*str = '\0'; break;
			}
		} else {
			if (n < MAXF) strcpy (scans[n++], argv[i]);
		}
	}
	printf ("n=%d\n", n);
	if (n < 2) usage (program);

/*******************************/
/* pre-process image filenames */
/*******************************/
	for (i = 0; i < n; i++) {
		while (str = strrchr (scans[i], '.')) {
				if (!strcmp (str, ".img"))  *str = '\0';
			else	if (!strcmp (str, ".ifh"))  *str = '\0';
			else	if (!strcmp (str, ".4dfp")) *str = '\0';
			else	break;
		}
	}

/*****************************************************/
/* read t4files and convert to rigid body parameters */
/*****************************************************/
	for (m = i = 0; i < n; i++) for (j = 0; j < n; j++) {
		if (i == j) continue;
		sprintf (t4file, "%s_to_%s_t4", scans[j], scans[i]);
		if (access (t4file, R_OK)) {
			fprintf (stderr, "%s: %s read error\n", program, t4file);
			exit (-1);
		}
		t4_read_ (t4file, t4);			/* t4_sub.f */
                t4_to_ck6_ (t4, That[m]);		/* CKsub.f */
		iget_t4_scale_ (t4file, &scale);	/* t4_sub.f */
		Shat[m] = log (scale);
                if (verbose) {
			printf ("%s\n", t4file);
			sprintf (command, "cat %s", t4file); system (command);
			printT (That[m], Shat + m);
		}
		if (!i) {
			Sbar[m] = Shat[m];
			for (k = 0; k < 6; k++) Tbar[m][k] = (zero_test) ?  Tzer[k] : That[m][k];
		}
		m++;
	}

/*******************************************/
/* compute symmetric transform discrepancy */
/*******************************************/
	symerr_ (&n, That[0]);

/************/
/* allocate */
/************/
	if (!(J = (float *) calloc (n*(n-1)*(n-1) * 36, sizeof (float)))) errm (program);

	printf ("%s: VOI rms radius=%.4f\n", program, radius);
	printf ("%s: begin Gauss-Newton trajectory estimation\n", program);
/******************************/
/* compute error and jacobian */
/******************************/
	do {	err0 = err;
		for (m = i = 0; i < n; i++) for (j = 0; j < n; j++) {
			if (i == j) continue;
			if (verbose) printf ("%3d%3d%3d %s_to_%s_t4\n", i, j, m, scans[j], scans[i]);
			if (!i) {		/* Txa */
				getjxy_ (Tbar[j-1], Tzer,      T,       That[m], Terr[m], J_xy_xa, J_xy_ya);
			} else if (!j) {	/* Tax */
				getjxy_ (Tzer,      Tbar[i-1], Tbar[m], That[m], Terr[m], J_xy_xa, J_xy_ya);
				Sbar[m] = -Sbar[i-1];
			} else {		/* Txy */
				getjxy_ (Tbar[j-1], Tbar[i-1], Tbar[m], That[m], Terr[m], J_xy_xa, J_xy_ya);
				Sbar[m] = Sbar[j-1] - Sbar[i-1];
			}
			fptr = (i) ? Tbar[m] : T;
			Serr[m] = Sbar[m] - Shat[m];
			if (verbose) {
				printT (fptr,    Sbar + m);
				printT (That[m], Shat + m);
				printT (Terr[m], Serr + m);
			}

			if (j) {
				Jptr = J + 6*m + 36*n*(n-1)*(j-1);
				fptr = J_xy_xa;
				for (l = 0; l < 6; l++) {
					for (k = 0; k < 6; k++) Jptr[k] = *fptr++;
					Jptr += 6*n*(n-1);
				}
			}
			if (i) {
				Jptr = J + 6*m + 36*n*(n-1)*(i-1);
				fptr = J_xy_ya;
				for (l = 0; l < 6; l++) {
					for (k = 0; k < 6; k++) Jptr[k] = *fptr++;
					Jptr += 6*n*(n-1);
				}
			}
			m++;
		}
		resolve_rigid_ (&radius, &n, J, Terr[0], dT[0], &err);
		for (i = 0; i < n - 1; i++) for (k = 0; k < 6; k++) Tbar[i][k] -= dT[i][k];
		resolve_scale_ (&n, Serr, dS);
		for (i = 0; i < n - 1; i++) Sbar[i] -= dS[i];
	printf ("err=%f err0=%f\n", err, err0);
	} while (fabs (err - err0) > 1.e-6*err0 && ++niter < 10);
	free (J);

/****************************************/
/* write resolved t4 files and sub file */
/****************************************/
	if (strlen (outroot)) {
		sprintf (subfile, "%s.sub", outroot);
		if (!(fpsub = fopen (subfile, "w"))) errw (program, subfile);
		printf ("Writing: %s\n", subfile);
		for (i = 0; i < n; i++) {
			if (i) {
				ck6_to_t4_ (Tbar[i-1], t4);
			} else {
				t4_init_ (t4);
			}
			sprintf (t4file, "%s_to_%s_t4", scans[i], outroot);
			t4_write_ (t4, t4file);
			scale = (i) ? exp (Sbar[i-1]) : 1.0;
			if (scale_flag) {
				if (!(fp_t4 = fopen (t4file, "a"))) errw (program, t4file);
				fprintf (fp_t4, "scale:    %10.6f\n", scale);
				fclose (fp_t4);
			}
			fprintf (fpsub, "%s\tt4=%s", scans[i], t4file);
			if (weight_flag) fprintf (fpsub, "\tweight=%.6f", 1.0 / scale);
			fprintf (fpsub, "\n");
		}
		fclose (fpsub);
		sprintf (command, "cat %s", subfile); system (command);
	}

/*****************************/
/* write mat file trajectory */
/*****************************/
	if (strlen (outroot) && mat_flag) {
/*****************************/
/* find maxlen scan filename */
/*****************************/
		l = 11; for (i = 0; i < n; i++) if ((k = strlen (scans[i])) > l) l = k; l++;
		sprintf (subfile, "%s.mat0", outroot);
		if (!(fpsub = fopen (subfile, "w"))) errw (program, subfile);
		printf ("Writing: %s\n", subfile);
		for (i = 0; i < n; i++) {
			if (i) {
				ck6_to_param6_ (Tbar[i-1], param6);
			} else {
				ck6_to_param6_ (Tzer, param6);
			}
			scale = (i) ? exp (Sbar[i-1]) : 1.0;
			k = strlen (scans[i]);
			for (j = 0; j < k; j++) fprintf (fpsub, "%.1s", scans[i] + j);
			for (; j < l; j++) fprintf (fpsub, "%.1s", " ");
			fprintf (fpsub, "%9.5f", scale);
			for (k = 0; k < 3; k++) fprintf (fpsub, "%9.4f", param6[k]);
			for (k = 3; k < 6; k++) fprintf (fpsub, "%9.4f", param6[k] * degperad);
			fprintf (fpsub, "\n");
		}
		fclose (fpsub);
		sprintf (command, "cat %s", subfile); system (command);
	}

	exit (status);
}
