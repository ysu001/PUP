/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/LDA.c,v 1.12 2011/09/15 02:28:20 avi Exp $*/
/*$Log: LDA.c,v $
 * Revision 1.12  2011/09/15  02:28:20  avi
 * trap nonsensical diemsnionality conditions
 *
 * Revision 1.11  2011/09/11  03:00:59  avi
 * options -c and -E
 *
 * Revision 1.10  2011/09/10  04:46:59  avi
 * minor corrections to shrinkage code
 *
 * Revision 1.9  2011/09/09  04:34:20  avi
 * option -s (shrinkage)
 *
 * Revision 1.8  2011/08/29  19:39:53  avi
 * increase dimensionality of feature space to 4096
 *
 * Revision 1.7  2011/08/25  00:18:59  avi
 * include classification results in output LDA files
 *
 * Revision 1.6  2011/08/24  04:21:41  avi
 * classification algebra working
 *
 * Revision 1.5  2011/08/23  23:05:35  avi
 * read/write option test dataset
 *
 * Revision 1.4  2011/08/22  04:31:51  avi
 * correct usage
 *
 * Revision 1.3  2011/08/22  04:02:21  avi
 * minor typos
 *
 * Revision 1.2  2011/08/22  03:35:39  avi
 * complete rewrite to read and write multiple files
 *
 * Revision 1.1  2011/08/21  23:24:54  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>		/* R_OK */
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <librms.h>

#define MAXL		256	/* maximum filename length */
#define MAXF		4096	/* maximum fnumber of input string fields */
#define MAXS		MAXF*16	/* maximum length of Xfile input line */
#define MAXC		16	/* maximum number of input files (classes) */

/********************/
/* global variables */
/********************/
static char	program[MAXL];
static char	rcsid[] = "$Id: LDA.c,v 1.12 2011/09/15 02:28:20 avi Exp $";

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

void getirange (char *string, int *minval, int *maxval) {
        char	*str;

	str = strstr (string, "to");
	if (str) {
		*str = '\0';
		*minval = atoi (string);
		*maxval = atoi (str + 2);
	} else {
		*minval = 0;
		*maxval = atoi (string);
	}
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
	fprintf (stderr, "Usage:	%s <(ASCII) input1> <(ASCII) input2> ...\n", program);
	fprintf (stderr, " e.g.:	%s setosa.txt versicolor.txt virginica.txt\n", program);
	fprintf (stderr, "	option\n");
	fprintf (stderr, "	-d	debug mode\n");
	fprintf (stderr, "	-s	apply shrinkage correction to within intertia matrix\n");
	fprintf (stderr,"	-c<int>[to[<int>] specify input columns to use (counting from 1; default all)\n");
	fprintf (stderr,"	-E<str> write discriminating eigenvetor[s] to named file\n");
	fprintf (stderr, "	-T<str>	classify test data in specified file\n");
	fprintf (stderr, "N.B.:	%s expects input with one line per item and dimen fields per line\n", program);
	fprintf (stderr, "N.B.:	each input file generates an output file named <input>.LDA\n", program);
	fprintf (stderr, "N.B.:	-s shrinkage follows Schafer & Strimmer 2005 Target B\n");
	exit (1);
}

int main (int argc, char *argv[]) {
	FILE		*fp, *fpt;
	char		Xfile[MAXC + 1][MAXL], tstfile[MAXL] = "", outfile[MAXL], Efile[MAXL] = "";
	double		**xbar, **B, **V, **lambda, **E, **X, **Y;
	double		***C, ***Cinv, **ybar, det[MAXC], **logp;
	double		nu, gammastar, vark;		/* shrinkage parameters */
	int		nclass, ntot, ntst, iclass, *jclass, dimen, dimenp, nitem[MAXC + 1], neig;
	int		col0 = 0, col1 = 0;

/***********/
/* utility */
/***********/
	double		q, u, t;
	int		c, i, j, k, l, m;
	char		*ptr, command[MAXL], string[MAXS];
	char		*srgv[MAXF];				/* string field pointers */

/*********/
/* flags */
/*********/
	int		debug = 0;
	int		shrink = 0;
	int		status = 0;
	int		Etest = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (nclass = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'd': debug++;				break;
				case 's': shrink++;				break;
				case 'E': strcpy (Efile, ptr);			*ptr = '\0'; break;
				case 'c': getirange (ptr, &col0, &col1);	*ptr = '\0'; break;
				case 'T': strcpy (tstfile, ptr);		*ptr = '\0'; break;
			}
		} else {
			if (nclass == MAXC) {
				fprintf (stderr, "%s: maximum number of input files is %d\n", program, MAXC);
				usage ();
			}
			strcpy (Xfile[nclass++], argv[i]);
		}	
	}
	if (nclass < 1) usage ();
	if (col0) col0--; if (col0 < 0 || (col1 && (col1 < 1))) usage ();
	if (col0 > col1) usage ();

/*************************************/
/* append test data to training data */
/*************************************/
	if (strlen (tstfile)) {
		strcpy (Xfile[nclass], tstfile);
		Etest++;
	}

	ntot = ntst = dimen = 0;
	for (iclass = 0; iclass < nclass + Etest; iclass++) {
		if (!(fp = fopen (Xfile[iclass], "r"))) errr (program, Xfile[iclass]);
		printf ("Reading: %s\n", Xfile[iclass]);
		while (fgets (string, MAXS, fp)) {
			m = split (string, srgv, MAXF);
			if (!m) continue;	/* skip blank lines */
			nitem[iclass]++;
			if (dimen || iclass) {
				if (dimen != m) {
					fprintf (stderr, "%s: %s line %d has %d fields\n",
							program, Xfile[iclass], nitem[iclass], m);
					status++;
				}
			} else {
				dimen = m;
				if (m < col1) {
					fprintf (stderr, "%s: input column count (%d) less than requested high column (%d)\n",
						program, m, col1);
					usage ();
				}
			}
		}
		printf ("%s is %d x %d\n", Xfile[iclass], nitem[iclass], dimen);
		if (iclass < nclass) {
			ntot += nitem[iclass];
		} else {
			ntst  = nitem[iclass];
		}
		if (fclose (fp)) errr (program, Xfile[iclass]);
	}
	if (status) exit (-1);
	printf ("ntot=%d\n", ntot);
	printf ("ntst=%d\n", ntst);

	if (col1) {
		dimen = col1 - col0;
	} else {
		col1 = dimen;
	}
	if (shrink && dimen < 2) {
		fprintf (stderr, "%s: shrinkage cannot be used with feature space of dimensionality 1\n", program);
		exit (-1);
	}
	dimenp = (dimen < 10) ? dimen : 10;	/* dimen used for screen printing */
	printf ("%s: using columns %d to %d\n", program, col0 + 1, col1);

/*********************/
/* allocate and read */
/*********************/
	X = calloc_double2 (ntot + ntst, dimen);
	for (j = iclass = 0; iclass < nclass + Etest; iclass++) {
		if (!(fp = fopen (Xfile[iclass], "r"))) errr (program, Xfile[iclass]);
		for (i = 0; i < nitem[iclass]; i++, j++) {
			fgets (string, MAXS, fp);
			m = split (string, srgv, MAXF);
			if (!m) continue;	/* skip blank lines */
			for (k = 0; k < dimen; k++) {
				X[j][k] = atof (srgv[k + col0]);
			}
		}
		if (fclose (fp)) errr (program, Xfile[iclass]);
	}
	if (debug) for (i = 0; i < dimenp; i++) {
		for (j = 0; j < ntot + ntst; j++) printf (" %9.4f", X[j][i]);
		printf ("\n");
	}

/************/
/* center X */
/************/
	xbar = calloc_double2 (nclass + 1, dimen);	/* xbar[nclass] will hold global mean */
	for (i = 0; i < dimen; i++) {
		for (q = j = 0; j < ntot; j++) q += X[j][i];
		xbar[nclass][i] = q/ntot;
		for (j = 0; j < ntot + ntst; j++) X[j][i] -= xbar[nclass][i];
	}
	if (0)  for (i = 0; i < dimenp; i++) {
		for (j = 0; j < ntot + ntst; j++) printf (" %9.4f", X[j][i]);
		printf ("\n");
	}

/***********************/
/* compute class means */
/***********************/
	for (i = 0; i < dimen; i++) {
		for (j = l = 0; l < nclass + Etest; l++) {
			for (q = k = 0; k < nitem[l]; k++) q += X[j++][i];
			xbar[l][i] = q/nitem[l];
		}
	}
	printf ("class means\n");
	for (i = 0; i < dimenp; i++) {
		for (l = 0; l < nclass; l++) printf (" %9.6f", xbar[l][i]);
		if (Etest) printf (" |%10.6f", xbar[l][i]);
		printf ("\n");
	}

/*******************/
/* between inertia */
/*******************/
	B = calloc_double2 (dimen, dimen);
	for (i = 0; i < dimen; i++) for (j = 0; j < dimen; j++) {
		for (q = iclass = 0; iclass < nclass; iclass++) {
			q += xbar[iclass][i]*xbar[iclass][j];
			B[i][j] = nitem[iclass]*q/ntot;
		}
	}
	printf ("between inertia\n");
	for (i = 0; i < dimenp; i++) {
		for (j = 0; j < dimenp; j++) printf (" %9.6f", B[j][i]);
		printf ("\n");
	}

/******************/
/* within inertia */
/******************/
	V = calloc_double2 (dimen, dimen);
	for (i = 0; i < dimen; i++) for (j = 0; j < dimen; j++) {
		for (q = k = 0; k < ntot; k++) q += X[k][i]*X[k][j];
		V[i][j] = q/ntot;
	}
/******************/
/* shrinkage code */
/******************/
	if (shrink) {
		for (vark = i = 0; i < dimen; i++) for (j = 0; j < dimen; j++) {
			for (u = k = 0; k < ntot; k++) {
				t = X[k][i]*X[k][j] - V[i][j];
				u += t*t;
			}
			vark += u/(ntot - 1);
		}
		for (nu = i = 0; i < dimen; i++) nu += V[i][i];
		nu /= dimen;		/* nu now has mean eigenvalue of uncorrected within intertia */
		for (q = i = 0; i < dimen; i++) for (j = 0; j < dimen; j++) {
			t = V[i][j];
			if (i == j) t -= nu;
			q += t*t;
		}
		t = ntot - 1;
		gammastar = vark*ntot/(q*t*t);
		printf ("%s shrinkage parameters: vark=%f nu=%f gammastar=%f\n", program, vark, nu, gammastar);
		for (i = 0; i < dimen; i++) for (j = 0; j < dimen; j++) V[i][j] *= 1. - gammastar;
		for (i = 0; i < dimen; i++) V[i][i] += nu*gammastar;
	}
	printf ("within inertia\n");
	for (i = 0; i < dimenp; i++) {
		for (j = 0; j < dimenp; j++) printf (" %9.6f", V[j][i]);
		printf ("\n");
	}

	neig = nclass - 1;
	if (!neig) goto DONE1;
	if (dimen < neig) {
		fprintf (stderr, "%s: dimensionality of feature space (%d) must exceed nclass (%d)\n",
				program, dimen, nclass);
		status = 1;
		goto DONE1; 
	}

/****************************************/
/* solve generalized eigenvalue problem */
/****************************************/
	lambda	= calloc_double2 (dimen, dimen);
	E	= calloc_double2 (dimen, dimen);
	dgeigen (&B[0][0], &V[0][0], &lambda[0][0], &E[0][0], &dimen);
	printf ("lambda\n");
	for (i = 0; i < dimenp; i++) {
		for (j = 0; j < dimenp; j++) printf (" %9.4f", lambda[j][i]);
		printf ("\n");
	}
	printf ("E\n");
	for (i = 0; i < dimenp; i++) {
		for (j = 0; j < neig; j++) printf (" %9.6f", E[j][i]);
		printf ("\n");
	}

/***************************************/
/* output discriminating eigenvetor[s] */
/***************************************/
	if (strlen (Efile)) {
		printf ("Writing: %s\n", Efile);
		if (!(fp = fopen (Efile, "w"))) errw (program, Efile);
		write_command_line (fp, argc, argv);
		fprintf (fp, "#discriminating eigenvetor");
		if (neig > 1) fprintf (fp, "s");
		fprintf (fp, "\n");
		for (i = 0; i < dimen; i++) {
			for (j = 0; j < neig; j++) fprintf (fp, " %9.6f", E[j][i]);
			fprintf (fp, "\n");
		}
		if (fclose (fp)) errr (program, outfile);
	}

/***********************/
/* compute Y = [Et][X] */
/***********************/
	Y = calloc_double2 (ntot + ntst, dimen);
	for (i = 0; i < ntot + ntst; i++) for (j = 0; j < dimen; j++) {
		for (q = k = 0; k < dimen; k++) q += E[j][k]*X[i][k];
		Y[i][j] = q;
	}

/*****************************/
/* compute ybar = [Et][xbar] */
/*****************************/
	ybar = calloc_double2 (nclass, neig);
	for (i = 0; i < nclass; i++) {
		for (j = 0; j < neig; j++) {
			for (q = k = 0; k < dimen; k++) q += E[j][k]*xbar[i][k];
			ybar[i][j] = q;
		}
		printf ("ybar[%d]=", i);
		for (k = 0; k < neig; k++) printf (" %9.4f", ybar[i][k]);
		printf ("\n");
	}

/*********************/
/* build classifiers */
/*********************/
	C	= calloc_double3 (nclass, neig, neig);
	Cinv	= calloc_double3 (nclass, neig, neig);
	for (i = 0; i < neig; i++) for (j = 0; j < neig; j++) {
		for (l = iclass = 0; iclass < nclass; iclass++) {
			for (q = k = 0; k < nitem[iclass]; k++, l++) {
				q += (Y[l][i] - ybar[iclass][i])*(Y[l][j] - ybar[iclass][j]);
			}
			Cinv[iclass][i][j] = C[iclass][i][j] = q/nitem[iclass];
		}
		assert (l == ntot);
	}
	for (iclass = 0; iclass < nclass; iclass++) {
		dmatinv_ (&Cinv[iclass][0][0], &neig, det + iclass);
		printf ("iclass=%d det=%f\n", iclass, det[iclass]);
	}

/*************************************************************/
/* evaluate class membership modeled as Gaussian distributed */
/*************************************************************/
	logp = calloc_double2 (ntot + ntst, nclass);
	if (!(jclass = (int *) calloc (ntot + ntst, sizeof (int)))) errm (program);
	for (i = 0; i < ntot + ntst; i++) for (iclass = 0; iclass < nclass; iclass++) {
		logp[i][iclass] = log ((double) nitem[iclass]/ntot) - 0.5*(neig*log (2.*M_PI) + log (det[iclass]));
		for (q = k = 0; k < neig; k++) for (j = 0; j < neig; j++) {
			q += (Y[i][k] - ybar[iclass][k])*Cinv[iclass][k][j]*(Y[i][j] - ybar[iclass][j]);
		}
		logp[i][iclass] -= 0.5*q;
		if (logp[i][iclass] > logp[i][jclass[i]]) jclass[i] = iclass;
	}

/************************/
/* write LDA components */
/************************/
	if (debug)  for (i = 0; i < ntot + ntst; i++) {
		printf ("%5d", i);
		for (k = 0; k < neig; k++) printf (" %9.4f", Y[i][k]);
		for (iclass = 0; iclass < nclass; iclass++) printf (" %9.4f", logp[i][iclass]);
		printf ("  class=%d\n", jclass[i]);
	}
	for (j = iclass = 0; iclass < nclass + Etest; iclass++) {
		sprintf (outfile, "%s.LDA", Xfile[iclass]);
		printf ("Writing: %s\n", outfile);
		if (!(fp = fopen (outfile, "w"))) errw (program, outfile);
		write_command_line (fp, argc, argv);
		for (i = 0; i < nitem[iclass]; i++, j++) {
			for (k = 0; k < neig; k++)   fprintf (fp, " %9.4f", Y[j][k]);
			for (l = 0; l < nclass; l++) fprintf (fp, " %9.4f", logp[j][l]);
			fprintf (fp, "  class=%d\n", jclass[j]);
		}
		if (fclose (fp)) errr (program, outfile);
	}
	assert (j == ntot + ntst);

	free_double3 (C); free_double3 (Cinv); free_double2 (ybar); free_double2 (logp); free (jclass);
	free_double2 (lambda); free_double2 (E); free_double2 (Y);
DONE1:	free_double2 (X); free_double2 (V); free_double2 (B); free_double2 (xbar);
	exit (status);
}
