/*$Header: /home/usr/shimonyj/diff4dfp/RCS/nonlinear.c,v 1.9 2012/12/19 00:40:19 shimonyj Exp $*/
/*$Log: nonlinear.c,v $
 * Revision 1.9  2012/12/19  00:40:19  shimonyj
 * minor fixes in nonlinear routine.
 *
 * Revision 1.8  2012/10/15  12:40:04  shimonyj
 * Change priors for constant and iso diffusion
 *
 * Revision 1.7  2012/09/02  20:39:02  shimonyj
 * minor error correction
 *
 * Revision 1.6  2012/09/02  17:51:18  shimonyj
 * correct a minor error
 *
 * Revision 1.5  2012/09/01  20:40:54  shimonyj
 * modify non linear model to include more complicated extra term
 * S = S0*exp(-b*DT) + C*exp(-bd)
 *
 * Revision 1.4  2008/09/16  02:57:23  avi
 * #include <JSSutil.h>
 *
 * Revision 1.3  2008/06/11  03:04:31  avi
 * restore carelessly commented out test on eigenvalue limits
 *
 * Revision 1.2  2008/05/29  22:46:06  avi
 * remove lprob_array from argument list
 *
 * Revision 1.1  2008/05/26  07:24:25  avi
 * Initial revision
 **/

/* perform non-linear fitting of tensor data */
/*
1. include the CO option flag
2. output the errors 
3. non convergent error bars
4. decrease number of population
5. output residuals
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <JSSutil.h>

/* LM parameters */
#define NPOP	50	/* number of population members */
#define MA	9	/* number of coefficients fits in function */
#define DELCHI	0.5	/* min change in logprob; original 0.1 */ 
#define EPS	0.01	/* tolerance for reset of itst  */
#define NITER	2	/* number of iterations; original 2 */
#define MAXITER 100
#define PRT_FLG 0	/* extra print flag */
#define LRGNEG	-1.0e9	/* large negative number */

int mrqprior(float x[], float y[],int ndata,float a[],int ia[],
		float aavg[],int ma,float **covar,float **alpha,float *logprob,
		float *alamda, float *bval, float **bvec);
void mrqpcof(float x[],float y[],int ndata,float a[],int ia[],
		float aavg[],int ma, float **alpha, float beta[],float *logprob,
		float *bval, float **bvec);
void dti_full(float x, float a[], float *y, float dyda[], 
		float *bval, float **bvec);
int gaussjordan(float **a, int n, float **b, int m);
int range_angle(float *a1, float *a2, float *a3);

/*
** range_angle - chagne angles back into standard range
**
** angles are called by reference, order alpha, beta, gamma
** normal return is 0, if error returns a 1
*/
int range_angle(float *a1, float *a2, float *a3)
{
	/* error checking */
	if (*a1 < -4.0*PI || *a1 > 4.0*PI) return(1);
	if (*a2 < -4.0*PI || *a2 > 4.0*PI) return(1);
	if (*a3 < -4.0*PI || *a3 > 4.0*PI) *a3 = 0.0;

	/* move angles between 0-2PI */
	while (*a1 < 0.0)	{ *a1 += 2.0*PI; }
	while (*a1 > 2.0*PI)	{ *a1 -= 2.0*PI; }
	while (*a2 < 0.0)	{ *a2 += 2.0*PI; }
	while (*a2 > 2.0*PI)	{ *a2 -= 2.0*PI; }
	while (*a3 < 0.0)	{ *a3 += 2.0*PI; }
	while (*a3 > 2.0*PI)	{ *a3 -= 2.0*PI; }

	/* next ensure that theta(beta) between 0-PI */
	/* simple, but causes flip in direction, usually ok for diffusion */
	if (*a2 > PI) *a2 += -PI;
	/* slightly longer, preserves direction, phi extends 0-2PI  */
	/* usually commented out 
	if (*a2 > PI) {
		*a2 = 2.0*PI - *a2;
		if (*a1 < PI) *a1 += PI;
		else *a1 += -PI;
	}
	*/

	/* next, limits phi(alpha) between 0-PI, may flip */
	if (*a1 >= PI) {
		*a1 += -PI;
		*a2 = PI - *a2;
	}
	/* next limit psi(gamma) between 0-PI */
	if (*a3 > PI) *a3 += -PI;

	return(0);
}


/*
** nonlin_dti - nonlinear fitting of diffusion tensor
** 	uses priors and returns estimated error bars for each parameter
** input: 
**	ndir - number of measurements
**	in - input data; range 1:ndir
**	bval - bvalues; range 1:ndir
**	bvec - bvectors; range [1:ndir][1:3]
**	badvox - bad voxel flag 0:off 1:on -1:debug mode
**	eval - eigenvalues from LLS [1:3]
**	ang - angles from LLS [1:3]
** output:
**	parval - returned parameter values; range [1:MA]
**	1: lambda 1
**	2: lambda 2
**	3: lambda 3
**	4: alpha/phi
**	5: beta/theta
**	6: gamma/psi
**	7: S0
**	8: C, additive constant offset, the formula: S = S0*exp(-b*DT) + C
**	9: d, additive isotropic diffusion term: S = S0*exp(-b*DT) + C*exp(-bd)
**	N.B.: 8 is activated with CO_flag=1, 8 and 9 activated with CO_flag=2
**	parstd - std of parameter values; range [1:MA]
*/
int nonlin_dti(int ndir, float *in, float *bval, float **bvec, int badvox, int CO_flag, float *eval, float *ang, 
	float *parval, float *parstd) 
{

	static long int seed = -127;
	FILE	*fp;
	int i,j,k,l,m,itst,*lista,iopt, niter, prtflg = PRT_FLG, nlprobopt;
	float a1,a2,a3,temp,dbar,inmax;
	float alamda,lprob,olprob,lprobopt,*dindx,*dmeas,**covar,**alpha, lprob_array[NPOP];
	float apar[MA+1],amin[MA+1],amax[MA+1],aavg[MA+1],asig[MA+1],aest[MA+1];

/*****************************************/
/* general set up for later calculations */
/*****************************************/
	if (prtflg) printf("Enter nonlin_dti\n");
	/* allocate for LM method */
	lista=ivector(1,MA);
	dindx = vector(1,ndir);
	dmeas = vector(1,ndir);
	covar=matrix(1,MA,1,MA);
	alpha=matrix(1,MA,1,MA);

	/* for testing problem voxels*/
	if (badvox == -1) prtflg = 2; 

/*********************************************/
/* calculate tensor using Non Linear Methods */
/*********************************************/
	inmax = 0.0;
	for (k=1; k<=ndir; k++) {
		dindx[k] = (float)k;
		dmeas[k] = in[k];
		if (in[k] > inmax) inmax = in[k];
	}

	/* set up the ranges */
	lista[1]=2; amin[1] = 0.0; amax[1] =   3.0; aavg[1] =   0.5; asig[1] =   0.3;  /* lambda1 */
	lista[2]=2; amin[2] = 0.0; amax[2] =   3.0; aavg[2] =   0.7; asig[2] =   0.3;  /* lambda2 */
	lista[3]=2; amin[3] = 0.0; amax[3] =   3.0; aavg[3] =   0.9; asig[3] =   0.3;  /* lambda3 */
	lista[4]=1; amin[4] = 0.0; amax[4] =  3.14; aavg[4] =  1.57; asig[4] =   0.5;  /* alpha */
	lista[5]=1; amin[5] = 0.0; amax[5] =  3.14; aavg[5] =  1.57; asig[5] =   0.5;  /* beta */
	lista[6]=1; amin[6] = 0.0; amax[6] =  3.14; aavg[6] =  1.57; asig[6] =   0.5;  /* gamma */
	lista[7]=2; amin[7] = 0.0; amax[7] =1.5*inmax; aavg[7] = inmax; asig[7] = 0.5*inmax;  /* S0 */
	lista[8]=0; amin[8] = 0.0; amax[8] =1.5*inmax; aavg[8] = 0.0;   asig[8] = 0.5*inmax;  /* const */
	lista[9]=0; amin[9] = 0.0; amax[9] =   4.0; aavg[9] = 0.0;   asig[9] = 1.0;  /* iso diffusion */

	/* turn on extra terms, only set aavg if in use */
	if (CO_flag==1) {
		lista[8] = 2;
		aavg[8] = 1.0*inmax;
	}
	if (CO_flag==2) {
		lista[8] = 2;
		aavg[8] = 0.50*inmax;
		lista[9] = 2;
		aavg[9] = 2.0;
	}


	/* angle transformation section */
	if (range_angle(&ang[1], &ang[2], &ang[3])) {
		printf("Warning(0): input angle out of range: %f %f %f\n",ang[1],ang[2],ang[3]);
		ang[1] = 1.57;
		ang[2] = 1.57;
		ang[3] = 1.57;
	}

	/* set up based on log linear least squares */
	aest[1] = eval[1];
	aest[2] = eval[2];
	aest[3] = eval[3];
	aest[4] = ang[1];
	aest[5] = ang[2];
	aest[6] = ang[3];
	aest[7] = aavg[7];
	aest[8] = aavg[8];
	aest[9] = aavg[9];
	if (eval[1] < amin[1] || eval[1] > amax[1] ||
	    eval[2] < amin[2] || eval[2] > amax[2] ||
	    eval[3] < amin[3] || eval[3] > amax[3] ) badvox = 1;
		
	if (badvox > 0) for(k=1;k<=MA;k++) aest[k] = aavg[k];	
	

	/* debug statement */
	if (prtflg > 1) {
		printf("**************************************************\n");
		for (k=1; k<=MA; k++) 
		printf("%d: max %f avg %f sig %f\n",k,amax[k],aest[k],asig[k]);
	}
	lprobopt = LRGNEG;
	for (l=0; l<NPOP; l++) {
		for(k=1;k<=MA;k++) { 
			if (lista[k]==0) apar[k] = aavg[k];
			if (lista[k]!=0) do {
				apar[k] = aest[k] + asig[k]*gasdev(&seed);
			} while (apar[k] >= amax[k] || apar[k] <=amin[k]);
		}

		/* sort eigenvalues */
		if (apar[1] > apar[3] && apar[1] > apar[2]) { 
			temp=apar[3]; apar[3]=apar[1]; apar[1]=temp; }
		else if (apar[2] > apar[3] && apar[2] > apar[1]) { 
			temp=apar[3]; apar[3]=apar[2]; apar[2]=temp; }
		if (apar[1] > apar[2]) { temp=apar[2]; apar[2]=apar[1]; apar[1]=temp; }

		alamda = -1;
		if ( mrqprior(dindx,dmeas,ndir,apar,lista,aavg,MA,covar,alpha,&lprob,&alamda,bval,bvec) ) continue;

		niter=1; itst=0;
		while(itst < NITER && niter < MAXITER) {
		    if (prtflg > 1) {
			printf("\nPop# %d Iter# %2d %2d lprb %13.4e lprobopt %13.4e %10s %9.2e\n",
				l,niter,itst,lprob,lprobopt,"alamda:",alamda);
			printf("%10s %10s %10s %10s %10s %10s %10s %10s %10s\n", 
				"a[1]","a[2]","a[3]","a[4]","a[5]","a[6]","a[7]","a[8]","a[9]");
			for (m=1;m<=MA;m++) printf("%11.3e",apar[m]);
			printf("\n");
		    }
		    /* sort eigenvalues */
	 	    if (apar[1] > apar[3] && apar[1] > apar[2]) { 
			temp=apar[3]; apar[3]=apar[1]; apar[1]=temp; }
		    else if (apar[2] > apar[3] && apar[2] > apar[1]) { 
			temp=apar[3]; apar[3]=apar[2]; apar[2]=temp; }
		    if (apar[1] > apar[2]) { temp=apar[2]; apar[2]=apar[1]; apar[1]=temp; }

		    niter++;
		    olprob=lprob;
		    if ( mrqprior(dindx,dmeas,ndir,apar,lista,aavg,MA,covar,alpha,&lprob,&alamda,bval, bvec) ) niter=MAXITER;
		/* more flexible option */
		    if (lprob - olprob < -EPS) itst=0;
		    else if (lprob-olprob < DELCHI && lprob-olprob > -EPS) itst++;
		/* original style 
		    if (lprob < olprob) itst=0;
		    else if (fabs(lprob-olprob) < DELCHI) itst++;
		*/
		    if (alamda < 1.0e-10 || niter >= MAXITER) break;
		}
		if (niter >= MAXITER) continue;

		alamda=0.0;
		if ( mrqprior(dindx,dmeas,ndir,apar,lista,aavg,MA,covar,alpha,&lprob,&alamda,bval,bvec) ) continue;

		/* angle transformation section */
		if (range_angle(&apar[4], &apar[5], &apar[6])) {
			/*
			printf("Warning(1): angle out of range: %f %f %f\n",apar[4],apar[5],apar[6]);
			*/
			continue;
		}

		if (prtflg) {
			printf("%2d lprob %10.3f %4d ",l, lprob, k);
			for (m=1;m<=MA;m++) printf("%11.3e ",apar[m]);
			printf("\n");
			printf("                      ");
			for (m=1;m<=MA;m++) {
				if (covar[m][m] >= 0.0) printf("(%10.3e)",sqrt(covar[m][m]));
				else printf("(%10.3e)",0.0);
			}
			printf("\n");
		}

		/* store best one so far */
		lprob_array[l]=lprob;
		if (lprob > lprobopt) { 
			lprobopt = lprob; 
			iopt = l; 
			for (k=1; k<=MA; k++) {
				parval[k] = apar[k];
				if (covar[k][k] > 0.0) parstd[k] = sqrt(covar[k][k]);
				else parstd[k] = 0.0;
			}
		}
	} /* end population loop */
	
	for (nlprobopt=l=0;l<NPOP; l++) if (fabs(lprobopt - lprob_array[l]) <= 20.0) nlprobopt++;
	parval[0] = (float)nlprobopt; parstd[0] = 0;

	/* deallocate LM  memory */
	free_matrix(alpha,1,MA,1,MA);
	free_matrix(covar,1,MA,1,MA);
	free_vector(dmeas,1,ndir);
	free_vector(dindx,1,ndir);
	free_ivector(lista,1,MA); 

	if (lprobopt == LRGNEG) return 1;
	else return 0;
}


/* full tensor single fiber model */
void dti_full(float x, float a[], float *y, float dyda[], float *bval, float **bvec)
{
	int i,j,k,l,m,n,indx;
	double sa,ca,sb,cb,sg,cg,ex,ey,ez;
	double drot,drotda,drotdb,drotdg,drotdl1,drotdl2,drotdl3;
	double rot[3][3],dmatd[3][3],dmat[3][3];
	double s2a,c2a,s2b,c2b,s2g,c2g;
	
	indx = (int) x;

	/* angel transformation section 
	while (a[4] < 0.0) { a[4] += 2.0*PI; }
	while (a[4] > 2.0*PI) { a[4] -= 2.0*PI; }
	while (a[5] < 0.0) { a[5] += 2.0*PI; }
	while (a[5] > 2.0*PI) { a[5] -= 2.0*PI; }
	while (a[6] < 0.0) { a[6] += 2.0*PI; }
	while (a[6] > 2.0*PI) { a[6] -= 2.0*PI; }
	if (a[5] > PI) a[5] += -PI;
	if (a[4] >= PI) {
		a[4] += -PI;
		a[5] = PI - a[5];
	}
	if (a[6] > PI) a[6] += -PI;
	*/

	sa = sin((double) a[4]); ca = cos((double) a[4]);
	sb = sin((double) a[5]); cb = cos((double) a[5]);
	sg = sin((double) a[6]); cg = cos((double) a[6]);
	s2a = sin((double)2.0*a[4]); c2a = cos((double)2.0*a[4]);
	s2b = sin((double)2.0*a[5]); c2b = cos((double)2.0*a[5]);
	s2g = sin((double)2.0*a[6]); c2g = cos((double)2.0*a[6]);
	
	/* form rotation matrix Arfken p.180 */
	rot[0][0] = cg*cb*ca - sg*sa;
	rot[0][1] = cg*cb*sa + sg*ca;
	rot[0][2] = -cg*sb;
	rot[1][0] = -sg*cb*ca - cg*sa;
	rot[1][1] = -sg*cb*sa + cg*ca;
	rot[1][2] = sg*sb;
	rot[2][0] = sb*ca;
	rot[2][1] = sb*sa;
	rot[2][2] = cb;
	
	/* init the diagonal matrix */
	for (i=0; i<3; i++) 
		for (j=0; j<3; j++) dmatd[i][j]=0.0;
	if (a[1]<=a[2] && a[2]<=a[3]) {
		dmatd[0][0] = a[1]; dmatd[1][1] = a[2]; dmatd[2][2] = a[3];
	}
	else if (a[1]<=a[3] && a[3]<=a[2]) {
		dmatd[0][0] = a[1]; dmatd[1][1] = a[3]; dmatd[2][2] = a[2];
	}
	else if (a[2]<=a[1] && a[1]<=a[3]) {
		dmatd[0][0] = a[2]; dmatd[1][1] = a[1]; dmatd[2][2] = a[3];
	}
	else if (a[2]<=a[3] && a[3]<=a[1]) {
		dmatd[0][0] = a[2]; dmatd[1][1] = a[3]; dmatd[2][2] = a[1];
	}
	else if (a[3]<=a[1] && a[1]<=a[2]) {
		dmatd[0][0] = a[3]; dmatd[1][1] = a[1]; dmatd[2][2] = a[2];
	}
	else if (a[3]<=a[2] && a[2]<=a[1]) {
		dmatd[0][0] = a[3]; dmatd[1][1] = a[2]; dmatd[2][2] = a[1];
	}
	
	/* rotate the tensor into the lab frame */
	for (m=0; m<3; m++)
	for (n=0; n<3; n++) {
		dmat[m][n]=0.0;
		for (k=0; k<3; k++)
		for (l=0; l<3; l++) {
			dmat[m][n] += rot[k][m]*dmatd[k][l]*rot[l][n];
		}
	}
	
/*********************************************************************************/
/* This next section explicitely reproduces the above dmat for testing purposes */
/*
	for (m=0; m<3; m++) for (n=0; n<3; n++)  printf("%d %d test %f\n",m,n,dmat[m][n]);
	for (m=0; m<3; m++) for (n=0; n<3; n++)  dmat[m][n]=0.0;
	dmat[0][0] += a[1]*( cg*cg*cb*cb*ca*ca - 2.0*sg*cg*cb*sa*ca + sg*sg*sa*sa);
	dmat[0][1] += a[1]*( cg*cg*cb*cb*sa*ca - sg*cg*cb*sa*sa + sg*cg*cb*ca*ca - sg*sg*sa*ca);
	dmat[0][2] += a[1]*( -cg*cg*sb*cb*ca + sg*cg*sb*sa);
	dmat[1][0] += a[1]*( cg*cg*cb*cb*sa*ca + sg*cg*cb*ca*ca - sg*cg*cb*sa*sa - sg*sg*sa*ca);
	dmat[1][1] += a[1]*( cg*cg*cb*cb*sa*sa + 2.0*sg*cg*cb*sa*ca + sg*sg*ca*ca);
	dmat[1][2] += a[1]*( -cg*cg*sb*cb*sa - sg*cg*sb*ca);
	dmat[2][0] += a[1]*( -cg*cg*sb*cb*ca + sg*cg*sb*sa);
	dmat[2][1] += a[1]*( -cg*cg*sb*cb*sa - sg*cg*sb*ca);
	dmat[2][2] += a[1]*( cg*cg*sb*sb);

	dmat[0][0] += a[2]*( sg*sg*cb*cb*ca*ca + 2.0*sg*cg*cb*sa*ca + cg*cg*sa*sa);
	dmat[0][1] += a[2]*( sg*sg*cb*cb*sa*ca + sg*cg*cb*sa*sa - sg*cg*cb*ca*ca - cg*cg*sa*ca);
	dmat[0][2] += a[2]*( -sg*sg*sb*cb*ca - sg*cg*sb*sa);
	dmat[1][0] += a[2]*( sg*sg*cb*cb*sa*ca - sg*cg*cb*ca*ca + sg*cg*cb*sa*sa - cg*cg*sa*ca);
	dmat[1][1] += a[2]*( sg*sg*cb*cb*sa*sa - 2.0*sg*cg*cb*sa*ca + cg*cg*ca*ca);
	dmat[1][2] += a[2]*( -sg*sg*sb*cb*sa + sg*cg*sb*ca);
	dmat[2][0] += a[2]*( -sg*sg*sb*cb*ca - sg*cg*sb*sa);
	dmat[2][1] += a[2]*( -sg*sg*sb*cb*sa + sg*cg*sb*ca);
	dmat[2][2] += a[2]*( sg*sg*sb*sb);

	dmat[0][0] += a[3]*( sb*sb*ca*ca);
	dmat[0][1] += a[3]*( sb*sb*sa*ca);
	dmat[0][2] += a[3]*( sb*cb*ca);
	dmat[1][0] += a[3]*( sb*sb*sa*ca);
	dmat[1][1] += a[3]*( sb*sb*sa*sa);
	dmat[1][2] += a[3]*( sb*cb*sa);
	dmat[2][0] += a[3]*( sb*cb*ca);
	dmat[2][1] += a[3]*( sb*cb*sa);
	dmat[2][2] += a[3]*( cb*cb);
	for (m=0; m<3; m++) for (n=0; n<3; n++)  printf("%d %d test %f\n",m,n,dmat[m][n]);
*/
/*************************************************************************************/

	/* calc the basic drot element for later calculations */
	ex = bvec[indx][1]*dmat[0][0] + bvec[indx][2]*dmat[0][1] + bvec[indx][3]*dmat[0][2];
	ey = bvec[indx][1]*dmat[1][0] + bvec[indx][2]*dmat[1][1] + bvec[indx][3]*dmat[1][2];
	ez = bvec[indx][1]*dmat[2][0] + bvec[indx][2]*dmat[2][1] + bvec[indx][3]*dmat[2][2];
	drot = bvec[indx][1]*ex + bvec[indx][2]*ey + bvec[indx][3]*ez;

	/* derivative of dmat with respect to alpha */
	for (m=0; m<3; m++) for (n=0; n<3; n++)  dmat[m][n]=0.0;
	dmat[0][0] += a[1]*( cg*cg*cb*cb*(-s2a) - 2.0*sg*cg*cb*c2a + sg*sg*s2a);
	dmat[0][1] += a[1]*( cg*cg*cb*cb*c2a - sg*cg*cb*s2a + sg*cg*cb*(-s2a) - sg*sg*c2a);
	dmat[0][2] += a[1]*( -cg*cg*sb*cb*(-sa) + sg*cg*sb*ca);
	dmat[1][0] += a[1]*( cg*cg*cb*cb*c2a + sg*cg*cb*(-s2a) - sg*cg*cb*s2a - sg*sg*c2a);
	dmat[1][1] += a[1]*( cg*cg*cb*cb*s2a + 2.0*sg*cg*cb*c2a + sg*sg*(-s2a));
	dmat[1][2] += a[1]*( -cg*cg*sb*cb*ca - sg*cg*sb*(-sa));
	dmat[2][0] += a[1]*( -cg*cg*sb*cb*(-sa) + sg*cg*sb*ca);
	dmat[2][1] += a[1]*( -cg*cg*sb*cb*ca - sg*cg*sb*(-sa));
	dmat[2][2] += a[1]*(  0.0 );

	dmat[0][0] += a[2]*( sg*sg*cb*cb*(-s2a) + 2.0*sg*cg*cb*c2a + cg*cg*s2a);
	dmat[0][1] += a[2]*( sg*sg*cb*cb*c2a + sg*cg*cb*s2a - sg*cg*cb*(-s2a) - cg*cg*c2a);
	dmat[0][2] += a[2]*( -sg*sg*sb*cb*(-sa) - sg*cg*sb*ca);
	dmat[1][0] += a[2]*( sg*sg*cb*cb*c2a - sg*cg*cb*(-s2a) + sg*cg*cb*s2a - cg*cg*c2a);
	dmat[1][1] += a[2]*( sg*sg*cb*cb*s2a - 2.0*sg*cg*cb*c2a + cg*cg*(-s2a));
	dmat[1][2] += a[2]*( -sg*sg*sb*cb*ca + sg*cg*sb*(-sa));
	dmat[2][0] += a[2]*( -sg*sg*sb*cb*(-sa) - sg*cg*sb*ca);
	dmat[2][1] += a[2]*( -sg*sg*sb*cb*ca + sg*cg*sb*(-sa));
	dmat[2][2] += a[2]*( 0.0 );

	dmat[0][0] += a[3]*( sb*sb*(-s2a));
	dmat[0][1] += a[3]*( sb*sb*c2a);
	dmat[0][2] += a[3]*( sb*cb*(-sa));
	dmat[1][0] += a[3]*( sb*sb*c2a);
	dmat[1][1] += a[3]*( sb*sb*s2a);
	dmat[1][2] += a[3]*( sb*cb*ca);
	dmat[2][0] += a[3]*( sb*cb*(-sa));
	dmat[2][1] += a[3]*( sb*cb*ca);
	dmat[2][2] += a[3]*( 0.0 );

	ex = bvec[indx][1]*dmat[0][0] + bvec[indx][2]*dmat[0][1] + bvec[indx][3]*dmat[0][2];
	ey = bvec[indx][1]*dmat[1][0] + bvec[indx][2]*dmat[1][1] + bvec[indx][3]*dmat[1][2];
	ez = bvec[indx][1]*dmat[2][0] + bvec[indx][2]*dmat[2][1] + bvec[indx][3]*dmat[2][2];
	drotda = bvec[indx][1]*ex + bvec[indx][2]*ey + bvec[indx][3]*ez;

	/* derivative of dmat with respect to beta */
	for (m=0; m<3; m++) for (n=0; n<3; n++)  dmat[m][n]=0.0;
	dmat[0][0] += a[1]*( cg*cg*(-s2b)*ca*ca - 2.0*sg*cg*(-sb)*sa*ca );
	dmat[0][1] += a[1]*( cg*cg*(-s2b)*sa*ca - sg*cg*(-sb)*sa*sa + sg*cg*(-sb)*ca*ca );
	dmat[0][2] += a[1]*( -cg*cg*c2b*ca + sg*cg*cb*sa);
	dmat[1][0] += a[1]*( cg*cg*(-s2b)*sa*ca + sg*cg*(-sb)*ca*ca - sg*cg*(-sb)*sa*sa );
	dmat[1][1] += a[1]*( cg*cg*(-s2b)*sa*sa + 2.0*sg*cg*(-sb)*sa*ca );
	dmat[1][2] += a[1]*( -cg*cg*c2b*sa - sg*cg*cb*ca);
	dmat[2][0] += a[1]*( -cg*cg*c2b*ca + sg*cg*cb*sa);
	dmat[2][1] += a[1]*( -cg*cg*c2b*sa - sg*cg*cb*ca);
	dmat[2][2] += a[1]*( cg*cg*s2b);

	dmat[0][0] += a[2]*( sg*sg*(-s2b)*ca*ca + 2.0*sg*cg*(-sb)*sa*ca );
	dmat[0][1] += a[2]*( sg*sg*(-s2b)*sa*ca + sg*cg*(-sb)*sa*sa - sg*cg*(-sb)*ca*ca );
	dmat[0][2] += a[2]*( -sg*sg*c2b*ca - sg*cg*cb*sa);
	dmat[1][0] += a[2]*( sg*sg*(-s2b)*sa*ca - sg*cg*(-sb)*ca*ca + sg*cg*(-sb)*sa*sa );
	dmat[1][1] += a[2]*( sg*sg*(-s2b)*sa*sa - 2.0*sg*cg*(-sb)*sa*ca );
	dmat[1][2] += a[2]*( -sg*sg*c2b*sa + sg*cg*cb*ca);
	dmat[2][0] += a[2]*( -sg*sg*c2b*ca - sg*cg*cb*sa);
	dmat[2][1] += a[2]*( -sg*sg*c2b*sa + sg*cg*cb*ca);
	dmat[2][2] += a[2]*( sg*sg*s2b);

	dmat[0][0] += a[3]*( s2b*ca*ca);
	dmat[0][1] += a[3]*( s2b*sa*ca);
	dmat[0][2] += a[3]*( c2b*ca);
	dmat[1][0] += a[3]*( s2b*sa*ca);
	dmat[1][1] += a[3]*( s2b*sa*sa);
	dmat[1][2] += a[3]*( c2b*sa);
	dmat[2][0] += a[3]*( c2b*ca);
	dmat[2][1] += a[3]*( c2b*sa);
	dmat[2][2] += a[3]*( (-s2b));

	ex = bvec[indx][1]*dmat[0][0] + bvec[indx][2]*dmat[0][1] + bvec[indx][3]*dmat[0][2];
	ey = bvec[indx][1]*dmat[1][0] + bvec[indx][2]*dmat[1][1] + bvec[indx][3]*dmat[1][2];
	ez = bvec[indx][1]*dmat[2][0] + bvec[indx][2]*dmat[2][1] + bvec[indx][3]*dmat[2][2];
	drotdb = bvec[indx][1]*ex + bvec[indx][2]*ey + bvec[indx][3]*ez;

	/* derivative of dmat with respect to gamma */
	for (m=0; m<3; m++) for (n=0; n<3; n++)  dmat[m][n]=0.0;
	dmat[0][0] += a[1]*( (-s2g)*cb*cb*ca*ca - 2.0*c2g*cb*sa*ca + s2g*sa*sa);
	dmat[0][1] += a[1]*( (-s2g)*cb*cb*sa*ca - c2g*cb*sa*sa + c2g*cb*ca*ca - s2g*sa*ca);
	dmat[0][2] += a[1]*( -(-s2g)*sb*cb*ca + c2g*sb*sa);
	dmat[1][0] += a[1]*( (-s2g)*cb*cb*sa*ca + c2g*cb*ca*ca - c2g*cb*sa*sa - s2g*sa*ca);
	dmat[1][1] += a[1]*( (-s2g)*cb*cb*sa*sa + 2.0*c2g*cb*sa*ca + s2g*ca*ca);
	dmat[1][2] += a[1]*( -(-s2g)*sb*cb*sa - c2g*sb*ca);
	dmat[2][0] += a[1]*( -(-s2g)*sb*cb*ca + c2g*sb*sa);
	dmat[2][1] += a[1]*( -(-s2g)*sb*cb*sa - c2g*sb*ca);
	dmat[2][2] += a[1]*( (-s2g)*sb*sb);

	dmat[0][0] += a[2]*( s2g*cb*cb*ca*ca + 2.0*c2g*cb*sa*ca + (-s2g)*sa*sa);
	dmat[0][1] += a[2]*( s2g*cb*cb*sa*ca + c2g*cb*sa*sa - c2g*cb*ca*ca - (-s2g)*sa*ca);
	dmat[0][2] += a[2]*( -s2g*sb*cb*ca - c2g*sb*sa);
	dmat[1][0] += a[2]*( s2g*cb*cb*sa*ca - c2g*cb*ca*ca + c2g*cb*sa*sa - (-s2g)*sa*ca);
	dmat[1][1] += a[2]*( s2g*cb*cb*sa*sa - 2.0*c2g*cb*sa*ca + (-s2g)*ca*ca);
	dmat[1][2] += a[2]*( -s2g*sb*cb*sa + c2g*sb*ca);
	dmat[2][0] += a[2]*( -s2g*sb*cb*ca - c2g*sb*sa);
	dmat[2][1] += a[2]*( -s2g*sb*cb*sa + c2g*sb*ca);
	dmat[2][2] += a[2]*( s2g*sb*sb);

	dmat[0][0] += a[3]*( 0.0 );
	dmat[0][1] += a[3]*( 0.0 );
	dmat[0][2] += a[3]*( 0.0 );
	dmat[1][0] += a[3]*( 0.0 );
	dmat[1][1] += a[3]*( 0.0 );
	dmat[1][2] += a[3]*( 0.0 );
	dmat[2][0] += a[3]*( 0.0 );
	dmat[2][1] += a[3]*( 0.0 );
	dmat[2][2] += a[3]*( 0.0 );

	ex = bvec[indx][1]*dmat[0][0] + bvec[indx][2]*dmat[0][1] + bvec[indx][3]*dmat[0][2];
	ey = bvec[indx][1]*dmat[1][0] + bvec[indx][2]*dmat[1][1] + bvec[indx][3]*dmat[1][2];
	ez = bvec[indx][1]*dmat[2][0] + bvec[indx][2]*dmat[2][1] + bvec[indx][3]*dmat[2][2];
	drotdg = bvec[indx][1]*ex + bvec[indx][2]*ey + bvec[indx][3]*ez;

	/* derivative of dmat with respect to lambda1 */
	dmat[0][0] = cg*cg*cb*cb*ca*ca - 2.0*sg*cg*cb*sa*ca + sg*sg*sa*sa;
	dmat[0][1] = cg*cg*cb*cb*sa*ca - sg*cg*cb*sa*sa + sg*cg*cb*ca*ca - sg*sg*sa*ca;
	dmat[0][2] = -cg*cg*sb*cb*ca + sg*cg*sb*sa;
	dmat[1][0] = cg*cg*cb*cb*sa*ca + sg*cg*cb*ca*ca - sg*cg*cb*sa*sa - sg*sg*sa*ca;
	dmat[1][1] = cg*cg*cb*cb*sa*sa + 2.0*sg*cg*cb*sa*ca + sg*sg*ca*ca;
	dmat[1][2] = -cg*cg*sb*cb*sa - sg*cg*sb*ca;
	dmat[2][0] = -cg*cg*sb*cb*ca + sg*cg*sb*sa;
	dmat[2][1] = -cg*cg*sb*cb*sa - sg*cg*sb*ca;
	dmat[2][2] = cg*cg*sb*sb;

	ex = bvec[indx][1]*dmat[0][0] + bvec[indx][2]*dmat[0][1] + bvec[indx][3]*dmat[0][2];
	ey = bvec[indx][1]*dmat[1][0] + bvec[indx][2]*dmat[1][1] + bvec[indx][3]*dmat[1][2];
	ez = bvec[indx][1]*dmat[2][0] + bvec[indx][2]*dmat[2][1] + bvec[indx][3]*dmat[2][2];
	drotdl1 = bvec[indx][1]*ex + bvec[indx][2]*ey + bvec[indx][3]*ez;

	/* derivative of dmat with respect to lambda2 */
	dmat[0][0] = sg*sg*cb*cb*ca*ca + 2.0*sg*cg*cb*sa*ca + cg*cg*sa*sa;
	dmat[0][1] = sg*sg*cb*cb*sa*ca + sg*cg*cb*sa*sa - sg*cg*cb*ca*ca - cg*cg*sa*ca;
	dmat[0][2] = -sg*sg*sb*cb*ca - sg*cg*sb*sa;
	dmat[1][0] = sg*sg*cb*cb*sa*ca - sg*cg*cb*ca*ca + sg*cg*cb*sa*sa - cg*cg*sa*ca;
	dmat[1][1] = sg*sg*cb*cb*sa*sa - 2.0*sg*cg*cb*sa*ca + cg*cg*ca*ca;
	dmat[1][2] = -sg*sg*sb*cb*sa + sg*cg*sb*ca;
	dmat[2][0] = -sg*sg*sb*cb*ca - sg*cg*sb*sa;
	dmat[2][1] = -sg*sg*sb*cb*sa + sg*cg*sb*ca;
	dmat[2][2] = sg*sg*sb*sb;

	ex = bvec[indx][1]*dmat[0][0] + bvec[indx][2]*dmat[0][1] + bvec[indx][3]*dmat[0][2];
	ey = bvec[indx][1]*dmat[1][0] + bvec[indx][2]*dmat[1][1] + bvec[indx][3]*dmat[1][2];
	ez = bvec[indx][1]*dmat[2][0] + bvec[indx][2]*dmat[2][1] + bvec[indx][3]*dmat[2][2];
	drotdl2 = bvec[indx][1]*ex + bvec[indx][2]*ey + bvec[indx][3]*ez;

	/* derivative of dmat with respect to lambda3 */
	dmat[0][0] = sb*sb*ca*ca;
	dmat[0][1] = sb*sb*sa*ca;
	dmat[0][2] = sb*cb*ca;
	dmat[1][0] = sb*sb*sa*ca;
	dmat[1][1] = sb*sb*sa*sa;
	dmat[1][2] = sb*cb*sa;
	dmat[2][0] = sb*cb*ca;
	dmat[2][1] = sb*cb*sa;
	dmat[2][2] = cb*cb;

	ex = bvec[indx][1]*dmat[0][0] + bvec[indx][2]*dmat[0][1] + bvec[indx][3]*dmat[0][2];
	ey = bvec[indx][1]*dmat[1][0] + bvec[indx][2]*dmat[1][1] + bvec[indx][3]*dmat[1][2];
	ez = bvec[indx][1]*dmat[2][0] + bvec[indx][2]*dmat[2][1] + bvec[indx][3]*dmat[2][2];
	drotdl3 = bvec[indx][1]*ex + bvec[indx][2]*ey + bvec[indx][3]*ez;

	/* put it all together */
	*y = a[7]*exp(-bval[indx]*drot) + a[8]*exp(-bval[indx]*a[9]);
	dyda[1] = a[7]*exp(-bval[indx]*drot)*(-bval[indx]*drotdl1);
	dyda[2] = a[7]*exp(-bval[indx]*drot)*(-bval[indx]*drotdl2);
	dyda[3] = a[7]*exp(-bval[indx]*drot)*(-bval[indx]*drotdl3);
	dyda[4] = a[7]*exp(-bval[indx]*drot)*(-bval[indx]*drotda);
	dyda[5] = a[7]*exp(-bval[indx]*drot)*(-bval[indx]*drotdb);
	dyda[6] = a[7]*exp(-bval[indx]*drot)*(-bval[indx]*drotdg);
	dyda[7] = exp(-bval[indx]*drot);
	dyda[8] = exp(-bval[indx]*a[9]);
	dyda[9] = -bval[indx]*a[8]*exp(-bval[indx]*a[9]);
}


/*
** mrqprior - improved Levenberg/Marquardt from mrqmin in NR
** 	added relaxation parameter for slower convergence 
**		changes involve the relax parameter
**	added ability to use priors for the parameter estimation
*/
int mrqprior(float x[], float y[],int ndata,float a[],int ia[],
	float aavg[],int ma,float **covar,float **alpha,float *logprob,
	float *alamda, float *bval, float **bvec)
{
	int j,k,l,m;
	static int mfit;
	static float *da,*atry,**oneda,*beta,ologprob,relax;

	/* initialization */
	if (*alamda < 0.0) {
		atry = vector(1,ma);
		beta = vector(1,ma);
		da = vector(1,ma);
		for (mfit=0, j=1; j<=ma; j++) if (ia[j]) mfit++;
		oneda = matrix(1,mfit,1,1);
		*alamda = 0.001;
		mrqpcof(x,y,ndata,a,ia,aavg,ma,alpha,beta,logprob,bval,bvec);
		ologprob = (*logprob);
		for (j=1; j<=ma; j++) atry[j]=a[j];
		relax = 1.0e-6;
	}

	/* Alter fitting matrix by augmenting diagonal elements */
	for (j=0,l=1; l<=ma; l++) {
	    if (ia[l]) {
		for (j++,k=0,m=1; m<=ma; m++) {
			if (ia[m]) {
				k++; covar[j][k] = alpha[j][k];
			}
		}
		covar[j][j] = alpha[j][j]*(1.0+(*alamda));
		oneda[j][1] = beta[j];
	    }
	}

	/* Solve the matrix */

	if (gaussjordan(covar,mfit,oneda,1)) {
		free_matrix(oneda,1,mfit,1,1);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		free_vector(atry,1,ma);
		return(1); 
	}
/*
	gaussj(covar,mfit,oneda,1);
*/
	for (j=1; j<=mfit; j++) da[j] = oneda[j][1];

	/* finished, evaluate the covariance matrix */
	if (*alamda == 0.0) {
		covsrt(covar, ma,ia,mfit);
		free_matrix(oneda,1,mfit,1,1);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		free_vector(atry,1,ma);
		return(0);
	}

	/* Did the trial succed? */
	for (j=0, l=1; l<=ma; l++) {
		if (ia[l]) {
			atry[l] = a[l] + da[++j]*relax;
		}
	}
	mrqpcof(x,y,ndata,atry,ia,aavg,ma,covar,da,logprob,bval,bvec);

	/* Accept new solution */
	if (*logprob > ologprob) {
		relax = sqrt(relax);
		*alamda *= 0.1;
		ologprob = (*logprob);
		for (j=0, l=1; l<=ma; l++) {
		    if (ia[l]) {
			for (j++,k=0,m=1; m<=ma; m++) {
				if (ia[m]) {
					k++; alpha[j][k] = covar[j][k];
				}
			}
			beta[j] = da[j];
			a[l] = atry[l];
		    }
		}
	}
	/* failure, increase lambda and return */
	else {
		*alamda *= 10.0;
		*logprob = ologprob;
	}
	return(0); /* normal return */
}

/* 
** mrqpcof - designed to go with mrqprior
**	does not use sig[] but estimates it from residuals
**	ia[j]==2 add prior based on the following possible functions:
**		log(w/(w0^2 + w^2)) min is 0, peak at w0
**		1.0/(1.0 + exp(-h*(w - w_min))) min is w_min, no peak
**	ia[j]==1 do not add prior
**	ia[j]==0 do not search on this parameter
**	alpha and beta are NOT angels (see NR derivation)
**	factor of -1 is included in the beta calculation
**	There is a factor of 2 cancellation in the alpha, beta calcs
*/
void mrqpcof(float x[],float y[],int ndata,float a[],int ia[],
	float aavg[],int ma, float **alpha, float beta[],float *logprob, float *bval, float **bvec)
{
	int i,j,k,l,m,mfit=0;
	float ymod,sig2i,sig2,dy,den,*dyda;
	float h=2.0;  /* h is transition const for prior function 2 */

	dyda = vector(1,ma);

	/* initialization */
	for (j=1; j<=ma; j++) if (ia[j]) mfit++;
	for (j=1; j<=mfit; j++) {
		for (k=1; k<=j; k++) alpha[j][k] = 0.0;
		beta[j] = 0.0;
	}

	/* loop over data, calculate the likelihood function */
	sig2 = 0.0;
	for (i=1; i<=ndata; i++) {
		dti_full(x[i],a,&ymod,dyda,bval,bvec);
		dy = y[i] - ymod;
		sig2 += dy*dy;

		for (j=0, l=1; l<=ma; l++) {
		    if (ia[l]) {
			for (j++, k=0,m=1;m<=l;m++)
				if (ia[m]) alpha[j][++k] += dyda[l]*dyda[m];
			beta[j] += dy*dyda[l];
		    }
		}
	}
	(*logprob) = -((float)ndata/2.0)*log(sig2/2.0);

	/* normalize with 1/sig2 and add prior */
	sig2i = (float)ndata/sig2;
	for (j=0, l=1; l<=ma; l++) {
	    if (ia[l]) {
		j++;
		den = aavg[j]*aavg[j] + a[j]*a[j];
		for (k=0,m=1;m<=l;m++) {
			if (ia[m]) {
				++k;
				alpha[j][k] *= sig2i;
				/* function 1 2nd derivative */
				if (j==k && ia[j]==2) alpha[j][k] += (4.0*a[j]*a[j]/(den*den) - 2.0/den - 1.0/(a[j]*a[j]));
			}
		}
		beta[j] *= sig2i;

		/* function 1 1st derivative */
		if (ia[j]==2) beta[j] += -(1.0/a[j] - 2.0*a[j]/den);
		
		/* function 1 */
		if (ia[j]==2) (*logprob) += log(a[j]/den);
		
	    }
	}

	/* copy to symmetric matrix */
	for (j=2;j<=mfit; j++) 
		for (k=1; k<j; k++) alpha[k][j] = alpha[j][k];

	free_vector(dyda,1,ma);
}


/* invert a matrix using gauss jordan method NR p. 36 */
int gaussjordan(float **a, int n, float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i, icol,irow,j,k,l,ll;
	float big,dum,pivinv,temp;
	float swap;
	
	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for(i=1;i<=n;i++) {
		big=0.0;
		for(j=1;j<=n;j++) {
			if(ipiv[j] !=1) {
				for(k=1;k<=n;k++) {
					if(ipiv[k] == 0) {
						if(fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} 
					else if (ipiv[k] > 1) { /* singular matrix */
						free_ivector(ipiv,1,n);
						free_ivector(indxr,1,n);
						free_ivector(indxc,1,n);
						return(1);
					}
				}
			}
		}
		++(ipiv[icol]);
		if(irow != icol) {
			for(l=1;l<=n;l++) {swap=a[irow][l];a[irow][l]=a[icol][l];a[icol][l]=swap;}
			for(l=1;l<=m;l++) {swap=b[irow][l];b[irow][l]=b[icol][l];b[icol][l]=swap;}
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) { /* singular matrix */
			free_ivector(ipiv,1,n);
			free_ivector(indxr,1,n);
			free_ivector(indxc,1,n);
			return(1);
		}
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *=pivinv;
		for (l=1;l<=m;l++) b[icol][l] *=pivinv;
		for(ll=1;ll<=n;ll++) {
			if(ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for(l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for(l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
		}
	}
	for (l=n;l>=1;l--) {
		if(indxr[l] != indxc[l]) {
			for(k=1;k<=n;k++) {
				swap=a[k][indxr[l]];a[k][indxr[l]]=a[k][indxc[l]];a[k][indxc[l]]=swap;
			}
		}
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
	
	return(0); /* normal return */
}
