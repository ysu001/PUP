/*$Header: /data/petsun4/data1/src_solaris/permute/RCS/permute.c,v 1.1 2006/12/30 01:49:00 avi Exp $*/
/*$Log: permute.c,v $
 * Revision 1.1  2006/12/30  01:49:00  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXL		256

typedef struct {
	int	seq;
	long	val;
} SE;

static int secompare (const void *pi, const void *pj) {
	SE	*i, *j;

	i = (SE *) pi;
	j = (SE *) pj;
	return (i->val - j->val);
}

void permute (int *seq, int nseq) {
	SE		*list;
	int		k;
	int		debug = 0;

	list = (SE *) malloc (nseq * sizeof (SE));
	if (!list) {
		fprintf (stderr, "permute: memory allocation error\n");
		exit (-1);
	}
	for (k = 0; k < nseq; k++) {
		list [k].seq = k;
		list [k].val = random ();
	}
	if (debug) {fprintf (stdout, "call qsort...");  fflush (stdout);}
	qsort ((void *) list, nseq, sizeof (SE), secompare);
	if (debug) printf ("done\n");

	for (k = 0; k < nseq; k++) seq[k] = list[k].seq;
	free (list);
}

static char rcsid[] = "$Id: permute.c,v 1.1 2006/12/30 01:49:00 avi Exp $";
int permute_test (int argc, char *argv[]) {
	char		*ptr, program[MAXL];
	int		k, nseq = 20, *seq;
	unsigned int	iseed;

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);

	seq = (int *) malloc (nseq*sizeof (int));
	if (!seq) {
		printf ("%s: memory allocation error\n", argv[0]);
		exit (-1);
	}
	if (argc < 2) {
		printf ("Usage:	%s <iseed>\n", program);
		exit (-1);
	}
	iseed = atoi (argv[1]);
	srandom (iseed);
	permute (seq, nseq);
	for (k = 0; k < nseq; k++) printf ("%10d%10d\n", k, seq[k]);
	return 0;
}
