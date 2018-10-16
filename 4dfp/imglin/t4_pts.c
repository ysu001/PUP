/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4_pts.c,v 1.1 2009/10/27 05:44:20 avi Exp $*/
/*$Log: t4_pts.c,v $
 * Revision 1.1  2009/10/27  05:44:20  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <t4_io.h>

#define MAXL	256

/*************/
/* externals */
/*************/
void	t4inv_ (float *t4, float *t4inv);	/* param12opr.f */

/***********/
/* globals */
/***********/
static char	program[MAXL];
static char	rcsid[] = "$Id: t4_pts.c,v 1.1 2009/10/27 05:44:20 avi Exp $";

void errr (char* program, char* filespc) {
	fprintf (stderr, "%s: %s read error\n", program, filespc);
	exit (-1);
}

void errw (char* program, char* filespc) {
	fprintf (stderr, "%s: %s write error\n", program, filespc);
	exit (-1);
}

void errx (char* program, char* filespc) {
	fprintf (stderr, "%s: %s format error\n", program, filespc);
	exit (-1);
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
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

int main (int argc, char *argv[]) {
	FILE	*outfp, *pntfp;
	char	t4file[MAXL], pntfile[MAXL], outfile[MAXL] = "";
	float	t4[16], t4inv[16], x0[3], x1[3];

/***********/
/* utility */
/***********/
	char		*str, string[MAXL], command[MAXL];
	int 		c, i, k, m;

/*********/
/* flags */
/*********/
	int		debug = 0;
	int		verbose = 0;

	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case 'd': debug++;		break;
				case 'v': verbose++;		break;
			}
		} else switch (k) {
			case 0: strcpy (t4file,  argv[i]);	k++; break;
			case 1: strcpy (pntfile, argv[i]);	k++; break;
			case 2: strcpy (outfile, argv[i]);	k++; break;
		}
	}
	if (k < 2) {
		fprintf (stderr, "%s\n", rcsid);
		fprintf (stderr, "Usage:	%s <t4file> <pts.lst> [new pts.lst]\n", program);
		fprintf (stderr, " e.g.,	%s 711-2B_to_MNI152lin_T1_t4 711-2B_coords MNI152_coords\n", program);
		fprintf (stderr, "\toption\n");
		exit (1);
	}
	if (k > 2) fprintf (stdout, "%s\n", rcsid);

	if (t4_read (t4file, t4)) exit (-1);
	t4inv_ (t4, t4inv);

	if (!(pntfp = fopen (pntfile, "r"))) errr (program, pntfile);
	printf ("Reading: %s\n", pntfile);
	if (strlen (outfile)) {
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		printf ("Writing: %s\n", outfile);
	} else {
		outfp = stdout;
	}

	while (fgets (string, MAXL, pntfp) != NULL) {
		m = sscanf (string, "%f %f %f", x0 + 0, x0 + 1, x0 + 2);
		if (m != 3) errx (program, pntfile);

		for (k = 0; k < 3; k++) {
			x1[k] = t4inv[12 + k];
			for (i = 0; i < 3; i++) x1[k] += t4inv[k + 4*i]*x0[i];
		}
		fprintf (outfp, "%8.2f%8.2f%8.2f\n", x1[0], x1[1], x1[2]);
	}
	fclose (pntfp);
	fclose (outfp);

	exit (0);
}
