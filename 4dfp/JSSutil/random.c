/*$Header: /home/usr/shimonyj/JSSutil/RCS/random.c,v 1.1 2007/08/28 03:33:33 avi Exp $*/
/*$Log: random.c,v $
 * Revision 1.1  2007/08/28  03:33:33  avi
 * Initial revision
 **/

#include <math.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* ran1: random number generator NR page 280 */
/* must be called with a negative integer to init */
/* do not alter seed between calls */
float ran1(long *seed)
{
	int j;
	long k;
	static long	iy=0,iv[NTAB];
	float temp;
	
	if (*seed <=0 || !iy) {
		if (-(*seed) < 1) *seed=1;
		else *seed = -(*seed);
		for (j=NTAB+7; j>=0; j--) {
			k = (*seed)/IQ;
			*seed = IA*(*seed - k*IQ) - IR*k;
			if (*seed < 0) *seed += IM;
			if (j < NTAB) iv[j] = *seed;
		}
		iy = iv[0];
	}
	
	k=(*seed)/IQ;
	*seed = IA*(*seed - k*IQ) - IR*k;
	if (*seed < 0) *seed += IM;
	j = iy/NDIV;
	iy = iv[j];
	iv[j] = *seed;
	if ((temp=AM*iy) > RNMX) return (RNMX);
	else return (temp);
}

/* expdev: return deviate with exponential distribution unit mean */
/* uses ran1, NR page 287 */
float expdev(long *seed)
{
	float ran1(long *seed);
	float dum;
	
	do 
		dum = ran1(seed);
	while (dum == 0.0);
	return (-log(dum));
}

/* gasdev: return deviate with gaussian distribution zero mean, unit variance */
/* uses ran1, NR page 289 */
float gasdev(long *seed)
{
	float ran1(long *seed);
	static int iset=0;
	static float gset;
	float fac, rsq, v1, v2;
	
	if (iset == 0) {
		do {
			v1 = 2.0*ran1(seed) - 1.0;
			v2 = 2.0*ran1(seed) - 1.0;
			rsq = v1*v1 + v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		iset = 1;
		return (v2*fac);
	}
	else {
		iset=0;
		return (gset);
	}
}

