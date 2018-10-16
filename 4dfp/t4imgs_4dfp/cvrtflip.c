/*$Header: /data/petsun4/data1/src_solaris/t4imgs_4dfp/RCS/cvrtflip.c,v 1.2 2007/09/23 02:58:25 avi Exp $*/
/*$Log: cvrtflip.c,v $
 * Revision 1.2  2007/09/23  02:58:25  avi
 * #include <stdlib.h> and <stdio.h>
 *
 * Revision 1.1  2005/08/05  03:54:11  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void vrtflip (int *iori, int *imgdim, float *centeri, float *mmppixi, float *centert, float *mmppixt) {
	int	i, k;
	float flips[3][3] = {
		{-1., +1., -1.},	/* transverse */
		{-1., +1., +1.},	/* coronal */
		{+1., +1., +1.}		/* sagittal */
	};

	k = *iori - 2;		/* ifh.orientation */
	if (k < 0 || k > 2) {
		fprintf (stderr, "vrtflip: out of range orientation (%d)\n", *iori);
		exit (-1);
	}
	for (i = 0; i < 3; i++) {
		mmppixt[i] = mmppixi[i]*flips[k][i];
		centert[i] = centeri[i]*flips[k][i];
		if (flips[k][i] < 0.) centert[i] = mmppixt[i]*(imgdim[i] + 1) - centert[i];
	}
}
