/*$Header: /data/petsun4/data1/src_solaris/frame_align/RCS/slice_sequence.c,v 1.1 2011/07/26 05:24:02 avi Exp $*/
/*$Log: slice_sequence.c,v $
 * Revision 1.1  2011/07/26  05:24:02  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <math.h>

void slice_sequence (int interleaved, int force0, int reverse, int *seq, int nslice) {
	int		i, j, k, istart, inc;

	istart = (force0 || !interleaved) ? 0 : (nslice + 1) % 2;
	inc = (interleaved) ? 2 : 1;

	k = 0;
	i = istart;
	while (k < nslice) {
		seq[k++] = i;
		if ((i += inc) >= nslice) i = (istart + 1) % 2;
	}

	if (reverse) for (i = 0, j = nslice - 1; i < nslice/2; i++, j--) {
		k = seq[i];
		seq[i] = seq[j];
 		seq[j] = k;
	}
}
