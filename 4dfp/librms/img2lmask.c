/*$Header: /data/petsun4/data1/src_solaris/librms/RCS/img2lmask.c,v 1.1 2007/04/11 06:40:21 avi Exp $*/
/*$Log: img2lmask.c,v $
 * Revision 1.1  2007/04/11  06:40:21  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <librms.h>

static char rcsid [] = "$Id: img2lmask.c,v 1.1 2007/04/11 06:40:21 avi Exp $";
void img2lmask (int *pnx, int *pny, int *pnz, float *imag, int *mask, float *mmppix, float *pfhalf, float *pcrit) {
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

/***********************************************/
/* threshold and convert floating imgm to mask */
/***********************************************/
	for (vmax = i = 0; i < dimension; i++) if (imgm [i] > vmax) vmax = imgm [i];
	vmax *= *pcrit;
	for (i = 0; i < dimension; i++) {
		if (imgm [i] >= vmax) mask [i] = 1; else mask [i] = 0;
	}
	free (imgm);
}
