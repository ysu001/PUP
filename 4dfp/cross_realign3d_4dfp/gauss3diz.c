/*$Header: /data/petsun4/data1/src_solaris/cross_realign3d_4dfp/RCS/gauss3diz.c,v 1.2 2006/09/28 21:45:10 avi Exp $*/
/*$Log: gauss3diz.c,v $
 * Revision 1.2  2006/09/28  21:45:10  avi
 * Solaris 10
 *
 * Revision 1.1  1997/05/23  01:30:29  yang
 * Initial revision
 *
 * Revision 1.3  1996/05/14  02:43:41  avi
 * add Id, Log and header fields
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*************/
/* externals */
/*************/
extern int	npad_ (int *n, int *margin);									/* FORTRAN librms */
extern void	imgpad_   (float *imag, int *nx, int *ny, int *nz, float *imgp, int *nxp, int *nyp, int *nzp);	/* FORTRAN librms */
extern void	imgdap_   (float *imag, int *nx, int *ny, int *nz, float *imgp, int *nxp, int *nyp, int *nzp);	/* FORTRAN librms */
extern void	gauss3d_  (float *imag, int *nx, int *ny, int *nz, float *cmppix, float *fhalf);		/* FORTRAN librms */

static char rcsid[] = "$Id: gauss3diz.c,v 1.2 2006/09/28 21:45:10 avi Exp $";
void gauss3diz (float *imag, int *pnx, int *pny, int *pnz, float *mmppix, float *pfhalf) {
	float			*imagp;
	int			margin;
	int			nxp, nyp, nzp;

	nxp = *pnx;
	nyp = *pny;
	margin = (int) (2.0 * 0.1874 / (mmppix [2] * (*pfhalf)));
	nzp = npad_ (pnz, &margin);
	printf ("gauss3diz: nxp=%d nyp=%d nzp=%d\n", nxp, nyp, nzp);

	imagp = (float *) malloc (nxp * nyp * nzp * sizeof (float));
	imgpad_ (imag, pnx, pny, pnz, imagp, &nxp, &nyp, &nzp);
	gauss3d_ (imagp, &nxp, &nyp, &nzp, mmppix, pfhalf);
	imgdap_ (imag, pnx, pny, pnz, imagp, &nxp, &nyp, &nzp);
	free (imagp);
}
