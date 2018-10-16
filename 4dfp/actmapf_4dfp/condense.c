/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/condense.c,v 1.4 2008/01/29 01:33:34 avi Exp $*/
/*$Log: condense.c,v $
 * Revision 1.4  2008/01/29  01:33:34  avi
 * option -f (accept input format from file)
 * safe strncpy from command line input field
 * define exit status
 *
 * Revision 1.3  2006/10/02  01:56:05  avi
 * Solaris 10
 *
 * Revision 1.2  2005/08/30  07:05:38  avi
 * fix several bugs
 *
 * Revision 1.1  2005/08/30  03:34:18  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>

#define MAXL		256
#define MAXS		4096

void errm (char* program) {
	fprintf (stderr, "%s: memory allocation error\n", program);
	exit (-1);
}

void errr (char* program, char* filespc) {
	fprintf (stderr, "%s: %s read error\n", program, filespc);
	exit (-1);
}

static char rcsid[] = "$Id: condense.c,v 1.4 2008/01/29 01:33:34 avi Exp $";
void usage (char *program) {
	fprintf (stderr, "Usage:\t%s <format>\n", program);
	fprintf (stderr, " e.g.,\t%s \"4x86+4x86+4x86+4x86+4x86+4x86+4x86+4x86+4x86+\"\n", program);
	fprintf (stderr, "	options\n", program);
	fprintf (stderr, "	-v	verbose mode\n", program);
	fprintf (stderr, "	-f<str>	read input format string from specified file (default command line)\n", program);
	exit (1);
}

int main (int argc, char *argv[]) {
	FILE	*strfp;
	int	c, i, ii, j, k, l, m;
	int	len;
	int	*A, *B, maxs = MAXS;
	int	verbose = 0, alter, test;
	char	command[MAXL], *ptr, program[MAXL], strfile[MAXL] = "";
	char	*string, *head, *dupl, *tail, temp[64], template[32];

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);

	if (!(string =	(char *) malloc (MAXS*sizeof (char)))) errm (program);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
			case 'f': strcpy (strfile, ptr); k++;	*ptr = '\0'; break;
			case 'v': verbose++;			break;
			default:				break;
			}
		}
		else switch (k) {
			case 0:	strncpy (string, argv[i], MAXS);	k++; break;
			default:					k++; break;
		}	
	}
	if (k != 1) usage (program);

	if (strlen (strfile)) {
		if (!(strfp = fopen (strfile, "r"))) errr (program, strfile);
		fseek (strfp, 0l, SEEK_END);
		maxs = ftell (strfp) + 1;
		if (!(string = (char *) realloc ((void *) string, maxs*sizeof (char)))) errm (program);
		rewind (strfp);
		if (!(fgets (string, maxs, strfp))) errr (program, strfile);
		fclose (strfp);
/************************/
/* convert '\n' to NULL */
/************************/
		if (ptr = strpbrk (string, "\n")) *ptr = '\0';
	} else {
		if (string[MAXS - 1]) {
			fprintf (stderr, "%s: input string exceeds buffer limit (%d chars)\n", program, MAXS - 1);
			exit (-1);
		}
	}
	if (verbose) printf ("%s\n", rcsid);

	if (!(head	= (char *) malloc (maxs*sizeof (char)))) errm (program);
	if (!(dupl	= (char *) malloc (maxs*sizeof (char)))) errm (program);
	if (!(tail	= (char *) malloc (maxs*sizeof (char)))) errm (program);
	if (!(A		= (int *)  malloc (maxs*sizeof (int))))  errm (program);
	if (!(B		= (int *)  malloc (maxs*sizeof (int))))  errm (program);
	if (verbose) printf ("maxs=%d\n", maxs);
ITER:	if (verbose) printf ("\n%s\n", string);
	len = strlen (string);
	if (verbose) printf ("len=%d\n", len);
	do {
		for (l = 1; l <= 15; l++) {
			for (ii = 0; ii < len; ii++) {
				if (ii > 0 && isdigit (string[ii - 1])) continue;
				if (verbose) printf ("looking for templates of length %d\n", l);
				strncpy (template, string + ii, l); template[l] = '\0';
				if (isdigit (template[l - 1])) continue;
				for (test = k = i = 0; i < l; i++) {
					if (template[i] == '(') k++;
					if (template[i] == ')') k--;
					if (k < 0) test++;
				}
				if (k || test) continue;
				if (verbose) printf ("template=%s\n", template);
				for (k = 0; k < len; k++) A[k] = B[k] = 0;
				k = 0; for (i = ii; i < len; i++) {
					strncpy (temp, string + i, l); temp[l] = '\0';
					if (!strcmp (temp, template)) A[k++] = i;
				}
				m = k = 0; B[m++] = A[k++];
				while ((A[k] - A[k-1]) == l) B[m++] = A[k++];
				alter = 0; if (m > 1) {
					if (verbose) {
						printf ("duplicate sequence of length %d\n", m);
						for (j = 0; j < m; j++) printf (" %d", B[j]);
						printf ("\n");
					}
					head[0] = '\0'; strncpy (head, string, ii); head[ii] = '\0';
					if (verbose) printf ("head=%s\n", head);
					strncpy (dupl, string + ii, m*l); dupl[m*l] = '\0';
					if (verbose) printf ("dupl=%s\n", dupl);
					if (l == 1) {
						sprintf (dupl, "%d%s\0",   m, template);
					} else {
						sprintf (dupl, "%d(%s)\0", m, template);
					}
					if (verbose) printf ("dupl=%s\n", dupl);
					strcpy (tail, string + ii + m*l);
					if (verbose) printf ("tail=%s\n", tail);
					strcpy (string, head); strcat (string, dupl); strcat (string, tail);
					alter++;
					goto ITER;
				}
				ii += m*l - 1;
			}
		}
	} while (alter);
	printf ("%s", string);

	free (string); free (A); free (B); free (dupl); free (head); free (tail);
	exit (0);
}
