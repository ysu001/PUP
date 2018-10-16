/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/dnormal.c,v 1.1 2006/08/17 08:06:30 avi Exp $*/
/*$Log: dnormal.c,v $
 * Revision 1.1  2006/08/17  08:06:30  avi
 * Initial revision
 **/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXL		256

double dnormal () {
	static int	odd = 1;
	static double	r, theta;

	if (odd) {
		do {r = sqrt (-2.*log(drand48()));} while (isnan (r));
		theta = 2.*M_PI*(drand48() - .5);
		odd = 0;
		return r*cos(theta);
	} else {
		odd++;
		return r*sin(theta);
	}
}

static char	rcsid[] = "$Id: dnormal.c,v 1.1 2006/08/17 08:06:30 avi Exp $";
void dnormal_test (long iseed) {
	int		i;

	printf ("initilization %d\n", iseed); srand48 (iseed);
	for (i = 0; i < 10; i++) {
		printf ("%10d%10.6f\n", i, dnormal());
	}
}
