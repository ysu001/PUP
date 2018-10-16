/*$Header$*/
/*$Log$*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

extern void set_rnan_ (float *);

float rnan (void) {
	union {
		float		r;
		unsigned int	j;	/* works also with unsigned long at least on littleendian CPUs */ 
        } word;
	word.j = 0x7fffffff;
	return word.r;
}

int main (int argc, char *argv[]) {
	FILE		*fp;
	int		c, i, j, k, l, m;
	float		v[4] = {0., 1., -1., 0.};
	char		command[256];

	v[3] = rnan ();
	printf ("isnan(v[3])=%d\n", isnan(v[3]));
	set_rnan_ (v + 2);

	if (!(fp = fopen ("binary.out", "wb")) || fwrite (v, sizeof (float), 4, fp) != 4
	|| fclose (fp)) {
		fprintf (stderr, "write error\n");
		exit (-1);
	}
	printf  ("Writing: binary.out\n");
	sprintf (command, "od -f binary.out");
	printf  ("%s\n", command); system (command);
	exit (0);
}
