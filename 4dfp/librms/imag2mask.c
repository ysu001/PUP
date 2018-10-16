/*$Header: /data/petsun4/data1/src_solaris/librms/RCS/imag2mask.c,v 1.2 2007/09/18 23:23:52 avi Exp $*/
/*$Log: imag2mask.c,v $
 * Revision 1.2  2007/09/18  23:23:52  avi
 * update cast and #include librms.h
 *
 * Revision 1.1  1996/04/19  17:04:58  ty7777
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <librms.h>

static char rcsid[] = "$Id: imag2mask.c,v 1.2 2007/09/18 23:23:52 avi Exp $";
void imag2mask (int *pnx, int *pny, int *pnz, float *imag, short *mask, float *mmppix, float *pfhalf, float *pcrit) {
	float			vmax;
	float			*imgm;
	float			*imgmp;
	int			i;
	int			dimension;
	int			margin;
	int			nxp, nyp, nzp;

	dimension = (*pnx) * (*pny) * (*pnz);
	margin = (int) (2.0 * 0.1874 / (mmppix [0] * (*pfhalf))); nxp = npad_ (pnx, &margin);
	margin = (int) (2.0 * 0.1874 / (mmppix [1] * (*pfhalf))); nyp = npad_ (pny, &margin);
	margin = (int) (2.0 * 0.1874 / (mmppix [2] * (*pfhalf))); nzp = npad_ (pnz, &margin);

	imgmp = (float *) malloc (nxp * nyp * nzp * sizeof (float));
	imgpad_ (imag, pnx, pny, pnz, imgmp, &nxp, &nyp, &nzp);
	gauss3d_ (imgmp, &nxp, &nyp, &nzp, mmppix, pfhalf);
	imgm = (float *) malloc (dimension * sizeof (float));
	imgdap_ (imgm, pnx, pny, pnz, imgmp, &nxp, &nyp, &nzp);
	free (imgmp);

/*****************************************************/
/* threshold and convert floating imgm to short mask */
/*****************************************************/
	vmax = 0.;
	for (i = 0; i < dimension; i++) {if (imgm [i] > vmax) vmax = imgm [i];}
	vmax *= *pcrit;
	for (i = 0; i < dimension; i++) {
		if (imgm [i] >= vmax) mask [i] = 1; else mask [i] = 0;
	}
	free (imgm);
}
