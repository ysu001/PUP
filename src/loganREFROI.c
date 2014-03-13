/**********************************************************************************
 Calculating regional binding potential using Logan Reference Tissue model
  
 USAGE:
 loganREFROI reftac roitac roiname k2 sf ef opt
 
 reftac: reference region tac file
 roitac: target region tac file
 roiname: name of the roi
 k2: eflux rate constant for reference region
 sf: start frame
 ef: end frame
 opt: output file name

 Yi Su, 01/03/2014
**********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>
#include <librms.h>
#include <nrutil.h>
#include <RSF.h>
int main(int argc, char **argv)
{
	char reftacfn[MAXL], roitacfn[MAXL], optfn[MAXL], roiname[MAXL];
	int nframes, sf, ef, i, j, k;
	float *fptr, t[MAXL], st[MAXL], frd[MAXL], reftac[MAXL], roitac[MAXL], nv, k2;
	float bp, intc, slope, rsq;
	FILE *fp;
	

	/* Get inputs */
	if (argc<8) { 
		printf("USAGE:\n");
		printf("	loganREFROI reftac roitac roiname k2 st ef opt\n");
		exit(-1);
	}
	strcpy(reftacfn, argv[1]);
	strcpy(roitacfn, argv[2]);
	strcpy(roiname, argv[3]);
	k2=(float)atof(argv[4]);
	sf=atoi(argv[5]);
	ef=atoi(argv[6]);
	strcpy(optfn, argv[7]); 
	
	/* Read TACs */
	readtac(reftacfn, reftac, frd, st, &nv, &nframes);
	readtac(roitacfn, roitac, frd, st, &nv, &nframes);

	/* Calculate BP using logan reference tissue model */
	for (i=0; i<nframes; i++)
	{
		t[i]=st[i]+frd[i]/2.-st[0];
/*		printf("%f\n",t[i]); */
	}
/*	printf("%d\t%d\t%d\n",nframes,sf,ef);*/


	bp=loganREF(reftac, roitac, k2, t, frd, nframes, sf, ef, &intc, &slope, &rsq);

	/* Write Output */
	if (!(fp =fopen(optfn, "w"))) errw("loganREFROI", optfn);
	fprintf(fp,"@-> %s %s %s %f %d %d %s\n", argv[0], reftacfn, roitacfn, k2, sf, ef, optfn);
	fprintf(fp,"%-35s %11s %11s %11s %11s\n", "Structure_Name", "NVox", "bp", "R^2", "Intc");
	fprintf(fp,"%-35s %11.0f %11.4f %11.9f %11.4f\n",roiname,nv,bp,rsq,intc);
	fclose(fp);
}
