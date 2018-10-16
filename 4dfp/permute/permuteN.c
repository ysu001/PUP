/*$Header: /data/petsun4/data1/src_solaris/permute/RCS/permuteN.c,v 1.3 2010/03/17 03:52:10 avi Exp $*/
/*$Log*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXL		256

/*************/
/* externals */
/*************/
void permute (int *seq, int nseq);		/* permute.c */

static char rcsid[] = "$Id: permuteN.c,v 1.3 2010/03/17 03:52:10 avi Exp $";
int main (int argc, char *argv[]) {
	char		*ptr, program[MAXL], command[MAXL];
	int		nseq, i, j, k, c, *seq;
	unsigned int	iseed = 0;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);
/******************************/
/* get command line arguments */
/******************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 's': iseed = atoi (ptr);	*ptr = '\0'; break;
			}
		} else switch (k) {
		 	case 0: nseq = atoi (argv[i]);		k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage: %s <int n>\n", program);
		printf (" e.g., %s 31 -s10\n", program);
		printf ("\toption\n");
		printf ("\t-s<int>\tspecifcy randomization seed\n");
		printf ("N.B.:\tif a seed (option -n) is not specified or if the seed is 0 permutation is suppressed\n");
		printf ("N.B.:\t%s prints to stdout seed-dependent randomly permuted integers in the range [1 <n>]\n",
				program);
		exit (1);
	}

	seq = (int *) malloc (nseq*sizeof (int));
	if (!seq) {
		printf ("%s: memory allocation error\n", program);
		exit (-1);
	}
	if (iseed) {
		srandom (iseed);
		permute (seq, nseq);
	} else {
		for (k = 0; k < nseq; k++) seq[k] = k;
	}
	for (k = 0; k < nseq; k++) printf ("%d ", seq[k] + 1);

	free (seq);
	return 0;
}
