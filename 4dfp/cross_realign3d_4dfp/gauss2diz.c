/*$Header: /data/petsun4/data1/src_solaris/cross_realign3d_4dfp/RCS/gauss2diz.c,v 1.2 2006/09/28 21:45:31 avi Exp $*/
/*$Log: gauss2diz.c,v $
 * Revision 1.2  2006/09/28  21:45:31  avi
 * Solaris 10
 *
 * Revision 1.1  1997/05/23  01:30:29  yang
 * Initial revision
 *
 * Revision 1.1  1996/05/17  21:23:57  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
 
/*************/
/* externals */
/*************/
extern void	gauss3d_  (float *imag, int *nx, int *ny, int *nz, float *cmppix, float *fhalf);	/* FORTRAN librms */

static char rcsid[] = "$Id: gauss2diz.c,v 1.2 2006/09/28 21:45:31 avi Exp $";
void	gauss2diz (float *imag, int *pnx, int *pny, int *pnz, float *mmppix, float *pfhalf) {
	int			k;
	int			nzp;

	nzp = 1;
	printf ("gauss2diz: nxp=%d nyp=%d nzp=%d\n", *pnx, *pny, nzp);
	for (k = 0; k < *pnz; k++) {
		gauss3d_ (imag + (*pnx) * (*pny) * k, pnx, pny, &nzp, mmppix, pfhalf);
	}
}
