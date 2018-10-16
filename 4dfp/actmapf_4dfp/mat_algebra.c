/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/mat_algebra.c,v 1.4 2010/04/25 04:06:58 avi Exp $*/
/*$Log: mat_algebra.c,v $
 * Revision 1.4  2010/04/25  04:06:58  avi
 * options -J and -c
 *
 * Revision 1.3  2010/02/15  05:35:51  avi
 * options -I and -m
 *
 * Revision 1.2  2009/12/01  05:32:40  avi
 * option -p
 *
 * Revision 1.1  2009/03/05  03:48:43  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>		/* R_OK */
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#ifndef _FEATURES_H		/* defined by math.h in linux */
	#include <sunmath.h>	/* needed for isnormal() in Solaris but not linux */
#endif
#include <librms.h>

#define MAXL		256	/* maximum string length */
#define MAXM		128	/* maximum number of input matrices */
#define MAXS		4096	/* maximum length of matfile input line */

/********************/
/* global variables */
/********************/
static char	program[MAXL];
static char	rcsid[] = "$Id: mat_algebra.c,v 1.4 2010/04/25 04:06:58 avi Exp $";

void errm (char* program) {
	fprintf (stderr, "%s: memory allocation error\n", program);
	exit (-1);
}

void errr (char* program, char* filespc) {
	fprintf (stderr, "%s: %s read error\n", program, filespc);
	exit (-1);
}

void errw (char* program, char* filespc) {
	fprintf (stderr, "%s: %s write error\n", program, filespc);
	exit (-1);
}

double **calloc_double2 (int n1, int n2) {
	int	i;
	double	**a;

	if (!(a = (double **) malloc (n1 * sizeof (double *)))) errm (program);
	if (!(a[0] = (double *) calloc (n1 * n2, sizeof (double)))) errm (program);
	for (i = 1; i < n1; i++) a[i] = a[0] + i*n2;
	return a;
}

void free_double2 (double **a) {
	free (a[0]);
	free (a);
}

double ***calloc_double3 (int n1, int n2, int n3) {
	unsigned int	i, j;
	double		***a;

	if (!(a = (double ***) malloc (n1 * sizeof (double **)))) errm (program);
	if (!(a[0] = (double **) malloc (n1 * n2 * sizeof (double *)))) errm (program);
	if (!(a[0][0] = (double *) calloc (n1 * n2 * n3, sizeof (double)))) errm (program);
	for (i = 0; i < n1; i++) {
		a[i] = a[0] + n2*i;
		for (j = 0; j < n2; j++) {
			a[i][j] = a[0][0] + n3*(n2*i + j);
		}
	}
	return a;
}

void free_double3 (double ***a) {
	free (a[0][0]);
	free (a[0]);
	free (a);
}

static int split (char *string, char *srgv[], int maxp) {
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

int matdim (char* matfile) {
	FILE		*fp;
	int		m, n;
	char		string[MAXS];
	char		*srgv[MAXL];				/* string field pointers */

	n = 0;
	if (!(fp = fopen (matfile, "r"))) errr (program, matfile);
	while (fgets (string, MAXS, fp)) {
		m = split (string, srgv, MAXL);
		if (!m) continue;	/* skip blank lines */
		if (n) {
			if (n != m) {
				fprintf (stderr, "%: %s format error\n", program, matfile);
				exit (-1);
			}
		} else {
			n = m;
		}
	}
	return n;
}

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n#%s\n", rcsid);
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void usage () {
	fprintf (stderr, "Usage:	%s -<operation><outmat> <mat1> <mat2> ...\n", program);
	fprintf (stderr, "	operation\n");
	fprintf (stderr, "	a	add\n");
	fprintf (stderr, "	m	multiply (compute product)\n");
	fprintf (stderr, "	s	subtract (mat1 - mat2)\n");
	fprintf (stderr, "	p	solve simple eigenvalue problem |mat1 - lambda*I| = 0\n");
	fprintf (stderr, "	g	solve generalized eigenvalue problem |mat1 - lambda*mat2| = 0\n");
	fprintf (stderr, "	e	mean (expectation)\n");
	fprintf (stderr, "	v	variance\n");
	fprintf (stderr, "	option\n");
	fprintf (stderr, "	-F	Fisher r-to-z transform input\n");
	fprintf (stderr, "	-H	Fisher z-to-r transform output\n");
	fprintf (stderr, "	-I	invert output\n");
	fprintf (stderr, "	-J	report Kullback's J = tr[mat1*inv(mat2)] + tr[mat2*inv(mat1)] -2k\n");
	fprintf (stderr, "	-l<lst>	read input file names from specified list file\n");
	fprintf (stderr, "	-c<flt>	multiply first output matrix by specified scalar\n");
	fprintf (stderr, "N.B.:	input matrix dimensions must match\n");
	fprintf (stderr, "N.B.:	options -I and -H are mutually exclusive\n");
	fprintf (stderr, "N.B.:	the default operation is output <- input\n");
	exit (1);
}

int main (int argc, char *argv[]) {
	FILE		*fp, *fpt;
	char		matfile[MAXM][MAXL], lstfile[MAXL] = "", outfile[MAXL] = "", tmpfile[MAXL];
	int		nmat = 0, opr = 0, dim, nout;
	double		***mati, ***matt, scale = 1.;

/***********/
/* utility */
/***********/
	double		q, u, det;
	int		c, i, j, k, l, m;
	char		*ptr, command[MAXL], string[MAXS];
	char		*srgv[MAXL];				/* string field pointers */

/*********/
/* flags */
/*********/
	int		debug = 0;
	int		status = 0;
	int		Fisherin = 0;
	int		reportJ = 0;
	int		tanhout = 0;
	int		invertout = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);
if (0) {
	u = 1.; q = atanh(u); printf("q=%f\n", q);
	k = isnormal(q); printf("isnormal(q)=%d\n", k);
}
/************************/
/* process command line */
/************************/
	for (i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'd': debug++;				break;
				case 'F': Fisherin++;				break;
				case 'I': invertout++;				break;
				case 'H': tanhout++;				break;
				case 'a':
				case 'm':
				case 'e':
				case 'p':
				case 'g':
				case 's':
				case 'v': opr = c; strcpy (outfile, ptr);	*ptr =  '\0'; break;
				case 'J': reportJ++;				break;
				case 'l': strcpy (lstfile, ptr);		*ptr =  '\0'; break;
				case 'c': scale = atof (ptr);			*ptr =  '\0'; break;
			}
		} else {
			strcpy (matfile[nmat++], argv[i]);
		}	
	}

/*******************/
/* parse list file */
/*******************/
	if (strlen (lstfile)) {
		if (!(fp = fopen (lstfile, "r"))) errr (program, lstfile);
		while (fgets (command, MAXL, fp)) {
			if (nmat >= MAXM) {
				fprintf (stderr, "%s: maximum number of input matrices (%d) exceeded\n", program, MAXM);
				exit (-1);
			}
			m = split (command, srgv, MAXL);
			if (m != 1) continue;		/* valid list entries are one field */
			strcpy (matfile[nmat++], srgv[0]);
		}
		fclose (fp);
	}
	if (!nmat) usage ();
	printf ("nmat=%d\n", nmat);
	if (opr == 'p' && nmat > 1) {
		fprintf (stderr, "%s: simple diagonalization must be done one matrix at a time\n");
		usage ();
	}
	if (opr == 'g' && nmat != 2) {
		fprintf (stderr, "%s: generalized eigenvalue problem requires exactly 2 input matrices\n");
		usage ();
	}
	if (invertout && tanhout) usage ();

/*********************************/
/* check dimensional consistency */
/*********************************/
	for (i = 0; i < nmat; i++) {
		printf ("matfile[%d]=%s\n", i, matfile[i]);
		if (!strcmp (matfile[i], outfile)) {
			fprintf (stderr, "%s: output %s matches one or more inputs\n", program, outfile);
			exit (-1);
		}
		if (!i) {
			dim = matdim (matfile[i]);
			printf ("dim=%d\n", dim);
		} else {
			status |= (dim != matdim (matfile[i]));
		}
		if (status) {
			fprintf (stderr, "%s: %s %s dimension mismatch\n", program, matfile[0], matfile[i]);
			exit (-1);
		}
	}

/************/
/* allocate */
/************/
	mati = calloc_double3 (nmat, dim, dim);
	matt = calloc_double3 (4,    dim, dim);

/**************************************/
/* open file to collect comment lines */
/**************************************/
	sprintf (tmpfile, "%s.tmp", outfile);
	if (!(fpt = fopen (tmpfile, "w"))) errw (program, tmpfile);

/********/
/* read */
/********/
	for (i = 0; i < nmat; i++) {
		printf ("Reading: %s\n", matfile[i]); fflush (stdout);
		if (!(fp = fopen (matfile[i], "r"))) errr (program, matfile[i]);
		j = 0; while (fgets (string, MAXS, fp)) {
			if (string[0] == '#') fprintf (fpt, "#%s", string);
			m = split (string, srgv, MAXL);
			if (!m) continue;	/* skip blank lines */
			assert (m == dim);
			for (k = 0; k < dim; k++) {
				mati[i][j][k] = atof (srgv[k]);
				if (Fisherin) {
					q = mati[i][j][k];
					if (fabs(q) > 1.0) {
						fprintf (stderr, "%s: illegal input value (%.6f) in %s\n",
							program, q, matfile[i]);
						exit (-1);
					}
					mati[i][j][k] = atanh (q);
				}
			}
			j++;
		}
		assert (j == dim);
		if (fclose (fp)) errr (program, matfile[i]);
	}
	if (fclose (fpt)) errw (program, tmpfile);

/*******************/
/* execute algebra */
/*******************/
	for (l = 0; l < nmat; l++) {
		if (opr == 'm') {
			if (!l) for (j = 0; j < dim; j++) matt[0][j][j] = 1.0;
			dmatmul_ (&matt[0][0][0], &mati[l][0][0], &matt[1][0][0], &dim);
			dmatcop_ (&matt[1][0][0], &matt[0][0][0], &dim);
		}
		for (j = 0; j < dim; j++) for (i = 0; i < dim; i++) {
			switch (opr) {
				case 'a':
				case 'e':
					matt[0][j][i] += mati[l][j][i];			break;
				case 's':
					switch (l) {
						case 0: matt[0][j][i]  = mati[l][j][i]; break;
						case 1: matt[0][j][i] -= mati[l][j][i];
					}						break;
				case 'v':
					q = mati[l][j][i];
					matt[0][j][i] += q;
					matt[1][j][i] += q*q;				break;
				default: matt[0][j][i] = mati[0][j][i]; 		break;
			}
		}
	}

	q  = 1.0/nmat;
	for (j = 0; j < dim; j++) for (i = 0; i < dim; i++) {
		switch (opr) {
			case 'e':
				matt[0][j][i] *= q;
				break;
			case 'v':
				u = matt[0][j][i] * q;
				matt[0][j][i] = (matt[1][j][i] - u*u/q) / (nmat - 1);
				break;
			default:
				break;
		}
	}

	nout = 1;
	switch (opr) {
		case 'p':
			dmatcop_ (&mati[0][0][0], &matt[3][0][0], &dim);
			deigen_  (&mati[0][0][0], &matt[0][0][0], &dim);
			if (0) {	/* output is input if this code is executed */
/***************************************/
/* eigenvectors are columns of matt[0] */
/***************************************/
				dtranspos_ (&matt[0][0][0], &matt[1][0][0], &dim);
				dmatmul_   (&matt[0][0][0], &mati[0][0][0], &matt[2][0][0], &dim);
				dmatmul_   (&matt[2][0][0], &matt[1][0][0], &matt[0][0][0], &dim);
				printf ("factored output should equal input\n");
				for (j = 0; j < dim; j++) for (i = 0; i < dim; i++) {
					printf ("i=%2d j=%2d %9.6f %9.6f %9.6f\n",
						i, j, matt[3][j][i], matt[0][j][i], matt[3][j][i] - matt[0][j][i]);
				}
			}
			break;
		case 'g':
			dgeigen (&mati[0][0][0], &mati[1][0][0], &matt[0][0][0], &matt[1][0][0], &dim);
			nout = 2;
			if (0) {
/***************************************/
/* eigenvectors are columns of matt[1] */
/***************************************/
				dtranspos_ (&matt[1][0][0], &matt[2][0][0], &dim);
				dmatmul_   (&matt[1][0][0], &matt[2][0][0], &matt[3][0][0], &dim);
				printf ("eigenvector product should equal mat2\n");
				for (j = 0; j < dim; j++) for (i = 0; i < dim; i++) {
					printf ("i=%2d j=%2d %9.6f %9.6f %9.6f\n",
						i, j, matt[3][j][i], mati[1][j][i], matt[3][j][i] - mati[1][j][i]);
				}
			}
			break;
	}

	if (reportJ) {
		nout = 0;
		dmatcop_ (&mati[0][0][0], &matt[0][0][0], &dim);
		dmatcop_ (&mati[1][0][0], &matt[1][0][0], &dim);
		dmatinv_ (&matt[0][0][0], &dim, &det);
		dmatinv_ (&matt[1][0][0], &dim, &det);
		dmatmul_ (&mati[0][0][0], &matt[1][0][0], &matt[2][0][0], &dim);
		dmatmul_ (&mati[1][0][0], &matt[0][0][0], &matt[3][0][0], &dim);
		q = -2*dim;
		for (i = 0; i < dim; i++) q += matt[2][i][i] + matt[3][i][i];
		printf ("J= %.6f\n", q);
	}

if (scale != 1.) for (i = 0; i < dim; i++) for (j = 0; j < dim; j++) matt[0][j][i] *= scale;
/****************/
/* write result */
/****************/
if (strlen (outfile)) {
	printf ("Writing: %s\n", outfile);
	if (!(fp = fopen (outfile, "w"))) errw (program, outfile);
	write_command_line (fp, argc, argv);
	if (fclose (fp)) errw (program, outfile);
	sprintf (command, "cat %s >> %s", tmpfile, outfile);
	status = system (command);
	remove (tmpfile);
	if (!(fp = fopen (outfile, "a"))) errw (program, outfile);
} else {
	fp = stdout;
}

	if (opr == 'p') {
		fprintf (fp, "#");
		for (i = 0; i < dim; i++) fprintf (fp, "%9.4f ", mati[0][i][i]);
		fprintf (fp, "\n");
	}
	for (l = 0; l < nout; l++) {
		if (invertout) {
			dmatinv_ (&matt[l][0][0], &dim, &det);
			fprintf (fp, "#det= %.6g\n", det);
		}
		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++) {
				if (tanhout) {
					q = matt[0][j][i];
					matt[0][j][i] = (isnormal(q)) ? tanh (q) : 1.0;
				}
				fprintf (fp, " %9.6f", matt[l][j][i]);
			}
			fprintf (fp, "\n");
		}
	}
	if (strlen (outfile)) if (fclose (fp)) errw (program, outfile);

	free_double3 (mati); free_double3 (matt);
	exit (status);
}
