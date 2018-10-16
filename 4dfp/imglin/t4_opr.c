/*$Header*/
/*$Log: t4_opr.c,v $
 * Revision 1.1  2010/02/19  02:06:06  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <librms.h>
#include <t4_io.h>

#define MAXL	256	/* maximum string length */
#define MAXR	4096	/* maximum number of input t4file */

/********************/
/* global variables */
/********************/
char		program[MAXL];
static char	rcsid[] = "$Id: t4_opr.c,v 1.1 2010/02/19 02:06:06 avi Exp $";

void usage (char *program) {
	fprintf (stderr, "Usage:	%s <t4file_1> <t4file_2> ...\n", program);
	fprintf (stderr, "	option\n");
	fprintf (stderr, "	-e<out>	save arithmetic  mean\n");
	fprintf (stderr, "	-g<out>	save logarithmic mean\n");
	fprintf (stderr, "	-f	factor each input into rigid-body and strech components\n");
	fprintf (stderr, "	-l<lst>	read input file names from specified list file\n");
	fprintf (stderr, "N.B.:	options -e and -g should used only with isolated stretch t4 files\n");
	exit (1);
}

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

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;

	if (ptr = strchr (string, '#')) *ptr = '\0';
	i = m = 0;
	while (m < maxp) {
		while (!isgraph ((int) string[i]) && string[i]) i++;
		if (!string[i]) break;
		srgv[m++] = string + i;
		while (isgraph ((int) string[i])) i++;
		if (!string[i]) break;
		string[i++] = '\0';
	}
	return m;
}

void expS (float *S, float *E) {
	float	w[9], q[9], wt[9], a[9], b[9];
	int	i, j, three = 3;

	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		q[i + 3*j] = S[i + 4*j];	/* q <- S */
	}
	eigen_		(q, w,  &three);	/* w <- eigenvectors of S */
	transpos_	(w, wt, &three);

	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		q[i + 3*j] = (i == j) ? exp (q[i + 3*j]) : 0.;	/* q <- exp eigenvalues of (S) */
	}
	matmul_		(w, q,  a, &three);
	matmul_		(a, wt, b, &three);	/* b <- exp(S) */
	t4_init (E);
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		E[i + 4*j] = b[i + 3*j];
	}
}

void logS (float *S, float *L) {
	float	w[9], q[9], wt[9], a[9], b[9];
	int	i, j, three = 3;

	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		q[i + 3*j] = S[i + 4*j];	/* q <- S */
	}
	eigen_		(q, w,  &three);	/* w <- eigenvectors of S */
	transpos_	(w, wt, &three);

	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		q[i + 3*j] = (i == j) ? log (q[i + 3*j]) : 0.;	/* q <- ln eigenvalues of (S) */
	}
	matmul_		(w, q,  a, &three);
	matmul_		(a, wt, b, &three);	/* b <- log(S) */
	t4_init (L);
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		L[i + 4*j] = b[i + 3*j];
	}
}

void RSfactor (float *RS, float *R, float *S) {
	float	g[9], w[9], q[9], wt[9], a[9], b[9], det;
	int	i, j, three = 3;

	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		g[i + 3*j] = RS[i + 4*j];	/* g <- RS */
	}
	transpos_	(g, b,  &three);	/* b <- SRt */
	matmul_		(b, g,  q, &three);	/* q <- SS */
	eigen_		(q, w,  &three);	/* w <- eigenvectors of SS */
	transpos_	(w, wt, &three);

	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		q[i + 3*j] = (i == j) ? 1./sqrt (q[i + 3*j]) : 0.;	/* q <- eigenvalues of inv(S) */
	}
	matmul_		(w, q,  a, &three);
	matmul_		(a, wt, b, &three);	/* b <- inv(S) */
	matmul_		(g, b,  q, &three);	/* q <- R */
	matinv_		(b, &three, &det);	/* b <- S */
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		R[i + 4*j] = q[i + 3*j];
		S[i + 4*j] = b[i + 3*j];
	}
	for (i = 0; i < 3; i++) {
		R[i + 12] = RS[i + 12];		/* preserve displacement in returned R */
		S[i + 12] = 0.;
		R[3 + 4*i] = S[3 + 4*i] = 0.;	/* fourth row */
	}
	R[15] = S[15] = 1.;			/* always in t4 */
}

int main (int argc, char *argv[]) {
	FILE		*fp;
	char		t4file[MAXR][MAXL], lstfile[MAXL] = "";
	char		outfile[MAXL] = "", outroot[MAXL];
	float		t4[MAXR][16], scale[MAXR];
	float		RD[16], S[16], L[16], T[16], E[16];
	int		nt4 = 0;

/***********/
/* utility */
/***********/
	char		*ptr, command[MAXL];
	char		*srgv[MAXL];				/* list file field pointers */
	int		c, i, j, k, m, opr = 0;
	float		t;

/*********/
/* flags */
/*********/
	int		status = 0;
	int		debug = 0;

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);
/************************/
/* process command line */
/************************/
	for (nt4 = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'd': debug++;				break;
				case 'f': opr = c;				break;
				case 'e':
				case 'g': opr = c; strcpy (outfile, ptr);	*ptr = '\0'; break;
				case 'l': strcpy (lstfile, ptr);		*ptr = '\0'; break;
			}
		} else {
			if (nt4 < MAXR) strcpy (t4file[nt4++], argv[i]);
		}
	}

/*******************/
/* parse list file */
/*******************/
	if (strlen (lstfile)) {
		if (!(fp = fopen (lstfile, "r"))) errr (program, lstfile);
		while (fgets (command, MAXL, fp)) {
			m = split (command, srgv, MAXL);
			if (m < 1) continue;
			if (nt4 >= MAXR) {
				fprintf (stderr, "%s: maximum number of input t4file (%d) exceeded\n", program, MAXR);
				exit (-1);
			}
			strcpy (t4file[nt4++], srgv[0]);
		}
		fclose (fp);
	}
	if (!nt4) usage (program);
	printf ("nt4=%d\n", nt4);

/****************/
/* read t4files */
/****************/
	for (i = 0; i < nt4; i++) {
		t4_read (t4file[i], t4[i]);		/* t4_io.c */
		k = iget_t4_scale (t4file[i], &t);	/* t4_io.c */
		scale[i] = (k) ? 0. : t;
                if (debug) t4_write (stdout, t4[i], scale[i]);
	}

	for (i = 0; i < 16; i++) T[i] = 0.;
/***********/
/* process */
/***********/
	for (m = 0; m < nt4; m++) switch (opr) {
		case 'f':
			if (!(ptr = strrchr (t4file[m], '/'))) ptr = t4file[m]; else ptr++;
			strcpy 	(outroot, ptr);
			RSfactor (t4[m], RD, S);
 			sprintf (outfile, "%s_rigidbody", outroot);
			if (!(fp = fopen (outfile, "w"))) errw (program, outfile);
			printf ("Writing: %s\n", outfile);
			t4_write (fp, RD, scale[m]);
			if (fclose (fp)) errw (program, outfile);
		 	sprintf (outfile, "%s_stretch", outroot);
			if (!(fp = fopen (outfile, "w"))) errw (program, outfile);
			printf ("Writing: %s\n", outfile);
			t4_write (fp, S, 1.);
			if (fclose (fp)) errw (program, outfile);
			break;
		case 'e':
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) T[i + 4*j] += t4[m][i + 4*j];
			break;
		case 'g':
			logS (t4[m], L);
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) T[i + 4*j] += L[i + 4*j];
			break;
	}

/**********/
/* finish */
/**********/
	switch (opr) {
		case 'e': case 'g':
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) T[i + 4*j] /= nt4;
			break;
	}
	if (opr == 'g') {
			expS (T, E);
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) T[i + 4*j] = E[i + 4*j];
	}
	switch (opr) {
		case 'e': case 'g':
			if (!strlen (outfile)) usage (program);
			if (!(fp = fopen (outfile, "w"))) errw (program, outfile);
			printf ("Writing: %s\n", outfile);
			T[15] = 1.;
			t4_write (fp, T, 1.);
			if (fclose (fp)) errw (program, outfile);
			break;
	}

	exit (status);
}
