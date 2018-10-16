/*$Header: /home/usr/shimonyj/diff4dfp/RCS/dtensor.c,v 1.12 2008/04/14 01:32:53 adrian Exp $*/
/*$Log: dtensor.c,v $
 * Revision 1.12  2008/04/14  01:32:53  adrian
 * compute FA and store in ang[6] and in D_PARAMS_2
 *
 * Revision 1.11  2008/04/14  00:05:59  adrian
 * added inv_svd_thresh and removed inv_svd from eigen_calc
 *
 * Revision 1.10  2007/08/30  04:56:43  avi
 * JSSutil.h compliant
 *
 * Revision 1.9  2004/11/23  21:33:56  avi
 * fix memory leaks in inv_svd() and calc_xyz()
 *
 * Revision 1.8  2004/07/14  19:36:46  adrian
 * added dparam_calc_noswap that doesn't sort the eigenvectors
 *
 * Revision 1.7  2003/02/05  00:07:36  shimonyj
 * dparam_calc()
 *
 * Revision 1.6  2001/07/03  22:27:33  shimonyj
 * calc_xyz()
 *
 * Revision 1.5  2001/06/23  02:14:02  avi
 * suppress "set w_element 3 to zero" message
 *
 * Revision 1.4  2001/04/03  20:31:00  shimonyj
 * eliminate near singular error message
 *
 * Revision 1.3  2000/12/19  03:53:24  avi
 * #defines
 *
 * Revision 1.2  2000/12/14  22:57:50  shimonyj
 * copyright
 **/
/*******************************************************************************************/
/* Copyright 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008	 	   */
/* Washington University, Mallinckrodt Institute of Radiology.				   */
/* All Rights Reserved.                                                                    */
/* This software may not be reproduced, copied, or distributed without written             */
/* permission of Washington University. For further information contact J.S. Shimony.      */
/*******************************************************************************************/

extern double atan2_custom (double y, double x);

/*
functions to calculate the diffusion tensor
given the experimental measurements
list of updates:
10-5-00 eigen_calc modified to return dbar and asig as part of the ang array
10-6-00 inv_svd changed to return 1 if matrix was near singular, 0 otherwise
12-14-00 nih version, header file eliminated
NOTE: The following numerical recipes routines where used and need replacement
svdcmp, jacobi, vector, free_vector, matrix, free_matrix 
*/

/* 
** function: set_svd
** inputs: the experimental inputs
**	n - number of measurements, including I0
**	b[n] - array of b values for each measurement
** 	q[n][3] - array of vectors for each measurement
** output: a[n][7] matrix needed for diffusion calculation
** action: setup matrix needed to solve 
** the least square/max likelihood problem 
** notes: 
** 1. only needs to be calculated once per experiment
** 2. assumes the data arrays start at 1
** 3. assumes the b values are in units of sec/mm2 divided by 1000.0
** 4. assumes the q vectors are normalized
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <JSSutil.h>

void
set_svd (n, b, q, aa)
int	n; 
float	*b, **q, **aa;
{
	int i;

	for (i=1; i<=n; ++i) {
		aa[i][1] = -b[i] * q[i][1] * q[i][1]; 
		aa[i][2] = -b[i] * q[i][2] * q[i][2]; 
		aa[i][3] = -b[i] * q[i][3] * q[i][3]; 
		aa[i][4] = -2.0 * b[i] * q[i][1] * q[i][2]; 
		aa[i][5] = -2.0 * b[i] * q[i][1] * q[i][3]; 
		aa[i][6] = -2.0 * b[i] * q[i][2] * q[i][3]; 
		aa[i][7] = 1.0;
	}
}

/* 
** function: calc_svd
** inputs: the experimental inputs
**	n - number of measurements, including I0
**	in[n] - array of measurement values
** 	a[n][7] - array previously prepared in set_svd
** output: 
**	dvec[7] - array of the diffusion tensor results
** 1. dxx 2. dyy 3. dzz 4. dxy 5. dxz 6. dyz 7. ln(I0)
** action: set up and invert the svd problem, solve for tensor
** notes:  
** 1. assumes the data arrays start at 1
** 2. uses an approx for the st. dev.
*/

void 
calc_svd (n, in, aa, dvec)
int	n; 
float	*in, **aa, *dvec;
{
	int	i, j, k, imax;
	float   *ll, **a, **ainv;

	ll = vector (1, n);
	a = matrix (1, n, 1, 7);
	ainv = matrix (1, 7, 1, n);

/* create the ll vector with the results */
	for (i = 1; i <= n; ++i) {
		if (in[i] > 0.0) ll[i] = log(in[i]);
		else ll[i] = 1.0e-8;
/* 1st approx for standard deviation */
		ll[i] *= in[i];
	}

/* create the matrix to be inverted, including the 1st approx std */
	for (i=1; i<=n; ++i)
		for (j=1; j<=7; ++j)
			a[i][j] = aa[i][j] * in[i];

/* perform the svd inversion */
	inv_svd(a, n, 7, ainv);

/* calculate the diffusion values */
	for (i=1; i<=7; i++) {
		dvec[i] = 0.0;
		for (j=1; j<=n; j++) {
			dvec[i] += ainv[i][j] * ll[j];
		}
	}

	free_matrix (ainv, 1, 7, 1, n);
	free_matrix (a, 1, n, 1, 7);
	free_vector (ll, 1, n);
}

/* 
** function: eigen_calc
** inputs: dvec[7] - the diffusion vector as produced by calc_svd
** 1. dxx 2. dyy 3. dzz 4. dxy 5. dxz 6. dyz 7. ln(I0)
** output: 
** dinv[3][3] - the inverse diffusion matrix
** eigenval[3] - the eigenvalues, sorted by size
** eigenvec[3][3] - the eigenvectors, sorted by size, not transposed
** ang[3] - angular directions 
** ang[4] - returns dbar
** ang[5] - returns asig
** ang[6] - returns FA
** action: 
** notes:  
** 1. 
** 2.
*/

void
eigen_calc(dvec, dinv, eigenval, eigenvec, ang)
float *dvec, **dinv, *eigenval, **eigenvec, *ang;
{
	int	inv_svd();
	int	i,j, nj, imax;
	float 	**dten, max;

	dten = matrix(1,3,1,3);

/* storage of tensor for calculations */
	dten[1][1] = dvec[1]; dten[1][2] = dvec[4]; dten[1][3] = dvec[5];
	dten[2][1] = dvec[4]; dten[2][2] = dvec[2]; dten[2][3] = dvec[6];
	dten[3][1] = dvec[5]; dten[3][2] = dvec[6]; dten[3][3] = dvec[3];

/* Call to jacobi to calculate the eigenvalues */
	jacobi (dten, 3, eigenval, eigenvec, &nj);

/* sort eigenvalues by size, largest in 3 , permute the rest */
	/* find the largest, store in imax */
	max = -1.0e9;
	for (i=1; i<=3; ++i) 
	    if (eigenval[i]>=max) {
		max = eigenval[i];
		imax = i;
	    }
	/* swap into position 3 if needed */
	if (imax != 3) {
		max = eigenval[3];
		eigenval[3] = eigenval[imax];
		eigenval[imax] = max;
		for (i=1; i<=3; ++i) {
			max = eigenvec[i][3];
			eigenvec[i][3] = eigenvec[i][imax];
			eigenvec[i][imax] = max;
		}
	}
	/* swap 1 and 2 if needed */
	if (eigenval[2] < eigenval[1]) {
		max = eigenval[2];
		eigenval[2] = eigenval[1];
		eigenval[1] = max;
		for (i=1; i<=3; ++i) {
			max = eigenvec[i][2];
			eigenvec[i][2] = eigenvec[i][1];
			eigenvec[i][1] = max;
		}
	}

/* fill the ang array assume the eigenvec array is NOT transposed */
/* The angle order is alfa,beta,gamma from Arfken */
	ang[1] = atan2_custom(eigenvec[2][3],eigenvec[1][3]);
	ang[2] = acos(eigenvec[3][3]);
	ang[3] = atan2_custom(-1.0*eigenvec[3][2],eigenvec[3][1]);

/* in the case of beta=0.0 only alfa+gamma can be solved */
/* arbitrarily set alfa=0.0 */
	if (ang[2]==0.0) {
		ang[1] = 0.0;
		ang[3] = atan2_custom(eigenvec[2][1],eigenvec[1][1]);
	}

/* return other important diffusion values */
	ang[4] = (eigenval[1]+eigenval[2]+eigenval[3])/3.0;
	ang[5] = sqrt(((eigenval[1]-ang[4])*(eigenval[1]-ang[4])+
		       (eigenval[2]-ang[4])*(eigenval[2]-ang[4])+
		       (eigenval[3]-ang[4])*(eigenval[3]-ang[4]))/6.0)/ang[4];
	ang[6] = sqrt(((eigenval[1]-ang[4])*(eigenval[1]-ang[4])+
		       (eigenval[2]-ang[4])*(eigenval[2]-ang[4])+
		       (eigenval[3]-ang[4])*(eigenval[3]-ang[4]))/(2./3.)) 
		       / sqrt (eigenval[1]*eigenval[1] + eigenval[2]*eigenval[2] + eigenval[3]*eigenval[3]);
	free_matrix (dten, 1, 3, 1, 3);
}

typedef struct {
        float	ang[4];
	int 	imaj;
        float	Dbar;
        float	Asig;
	float 	FA;
        float	Anu;
        float	prolat;
        float	eta;
        float	eps;
        float	Amaj;
        float	Amin;
} D_PARAMS_2;

/* 
** function: dparam_calc
** derived from eigen_calc with the ang return modified to 
** return a structire with more complete diffusion data
** inputs: dvec[7] - the diffusion vector as produced by calc_svd
** 1. dxx 2. dyy 3. dzz 4. dxy 5. dxz 6. dyz 7. ln(I0)
** output: 
** dinv[3][3] - the inverse diffusion matrix
** eigenval[3] - the eigenvalues, sorted by size
** eigenvec[3][3] - the eigenvectors, sorted by size, not transposed
** dparam - D_PARAMS structure defined in tgrid11 with the follwing:
**	ang[3],Dbar,Asig,Anu,prolat,eta,eps,Amaj,Amin
*/
void
dparam_calc(float *dvec, float **dinv, float *eigenval, float **eigenvec, D_PARAMS *dparam)
{
	int	inv_svd();
	int i,j, nj, imax;
	float **dten, max, eta, eps, dbar;

	dten = matrix(1,3,1,3);

/* storage of tensor for calculations */
	dten[1][1] = dvec[1]; dten[1][2] = dvec[4]; dten[1][3] = dvec[5];
	dten[2][1] = dvec[4]; dten[2][2] = dvec[2]; dten[2][3] = dvec[6];
	dten[3][1] = dvec[5]; dten[3][2] = dvec[6]; dten[3][3] = dvec[3];

/* invert */
	inv_svd(dten, 3, 3, dinv); 

/* Call to jacobi to calculate the eigenvalues */
	jacobi (dten, 3, eigenval, eigenvec, &nj);

/* sort eigenvalues by size, largest in 3 , permute the rest */
	/* find the largest, store in imax */
	max = -1.0e9;
	for (i=1; i<=3; ++i) 
	    if (eigenval[i]>=max) {
		max = eigenval[i];
		imax = i;
	    }
	/* swap into position 3 if needed */
	if (imax != 3) {
		max = eigenval[3];
		eigenval[3] = eigenval[imax];
		eigenval[imax] = max;
		for (i=1; i<=3; ++i) {
			max = eigenvec[i][3];
			eigenvec[i][3] = eigenvec[i][imax];
			eigenvec[i][imax] = max;
		}
	}
	/* swap 1 and 2 if needed */
	if (eigenval[2] < eigenval[1]) {
		max = eigenval[2];
		eigenval[2] = eigenval[1];
		eigenval[1] = max;
		for (i=1; i<=3; ++i) {
			max = eigenvec[i][2];
			eigenvec[i][2] = eigenvec[i][1];
			eigenvec[i][1] = max;
		}
	}

/* fill the ang array assume the eigenvec array is NOT transposed */
/* The angle order is alfa,beta,gamma from Arfken */
	dparam->ang[1] = atan2_custom(eigenvec[2][3],eigenvec[1][3]);
	dparam->ang[2] = acos(eigenvec[3][3]);
	dparam->ang[3] = atan2_custom(-1.0*eigenvec[3][2],eigenvec[3][1]);

/* in the case of beta=0.0 only alfa+gamma can be solved */
/* arbitrarily set alfa=0.0 */
	if (dparam->ang[2]==0.0) {
		dparam->ang[1] = 0.0;
		dparam->ang[3] = atan2_custom(eigenvec[2][1],eigenvec[1][1]);
	}

/* return other important diffusion values */
	dparam->Dbar = dbar = (eigenval[1]+eigenval[2]+eigenval[3])/3.0;
	dparam->Asig = sqrt(((eigenval[1]-dbar)*(eigenval[1]-dbar)+
			(eigenval[2]-dbar)*(eigenval[2]-dbar)+
			(eigenval[3]-dbar)*(eigenval[3]-dbar))/6.0)/dbar;
	dparam->prolat = (2.0*eigenval[2]-eigenval[1]-eigenval[3])/(eigenval[3]-eigenval[1]);
	dparam->eps = eps = (eigenval[2] - eigenval[1])/2.0;
	dparam->Amin = dparam->eps / dparam->Dbar;
	dparam->eta = eta = (eigenval[3] - (eigenval[2]+eigenval[1])/2.0)/3.0;
	dparam->Amaj = dparam->eta / dparam->Dbar;
	dparam->Anu = pow(eta*(eta*eta - eps*eps), 1.0/3.0)/ dparam->Dbar;

	free_matrix (dten, 1, 3, 1, 3);
}

/* 
** function: dparam_calc_noswap
** derived from eigen_calc with the ang return modified to 
** return a structire with more complete diffusion data
** and no eigenvec sort 
** inputs: dvec[7] - the diffusion vector as produced by calc_svd
** 1. dxx 2. dyy 3. dzz 4. dxy 5. dxz 6. dyz 7. ln(I0)
** output: 
** dinv[3][3] - the inverse diffusion matrix
** eigenval[3] - the eigenvalues
** eigenvec[3][3] - the eigenvectors not transposed
** dparam - D_PARAMS_2 structure defined in tgrid11 with the follwing:
**	ang[3],imaj,Dbar,Asig,FA,Anu,prolat,eta,eps,Amaj,Amin
*/
void
dparam_calc_noswap(float *dvec, float **dinv, float *eigenval, float **eigenvec, D_PARAMS_2 *dparam)
{
	int	inv_svd();
	int i,j, nj, imax;
	float **dten, max, eta, eps, dbar;

	dten = matrix(1,3,1,3);

/* storage of tensor for calculations */
	dten[1][1] = dvec[1]; dten[1][2] = dvec[4]; dten[1][3] = dvec[5];
	dten[2][1] = dvec[4]; dten[2][2] = dvec[2]; dten[2][3] = dvec[6];
	dten[3][1] = dvec[5]; dten[3][2] = dvec[6]; dten[3][3] = dvec[3];

/* invert */
	inv_svd(dten, 3, 3, dinv); 

/* Call to jacobi to calculate the eigenvalues */
	jacobi (dten, 3, eigenval, eigenvec, &nj);

/* sort eigenvalues by size, largest in 3 , permute the rest */
	/* find the largest, store in imax */
	max = -1.0e9;
	for (i=1; i<=3; ++i) 
	    if (eigenval[i]>=max) {
		max = eigenval[i];
		imax = i;
	    }
	dparam->imaj = imax;

/* Angles not calculated becuase no sort */
	dparam->ang[1] = dparam->ang[2] = dparam->ang[3] = 0.;

/* return other important diffusion values */
	dparam->Dbar = dbar = (eigenval[1]+eigenval[2]+eigenval[3])/3.0;
	dparam->Asig = sqrt(((eigenval[1]-dbar)*(eigenval[1]-dbar)+
			(eigenval[2]-dbar)*(eigenval[2]-dbar)+
			(eigenval[3]-dbar)*(eigenval[3]-dbar))/6.0)/dbar;
	dparam->FA = sqrt(((eigenval[1]-dbar)*(eigenval[1]-dbar)+
			(eigenval[2]-dbar)*(eigenval[2]-dbar)+
			(eigenval[3]-dbar)*(eigenval[3]-dbar))/(2./3.)) 
			/ sqrt (eigenval[1]*eigenval[1] + eigenval[2]*eigenval[2] + eigenval[3]*eigenval[3]);
	dparam->prolat = (2.0*eigenval[2]-eigenval[1]-eigenval[3])/(eigenval[3]-eigenval[1]);
	dparam->eps = eps = (eigenval[2] - eigenval[1])/2.0;
	dparam->Amin = dparam->eps / dparam->Dbar;
	dparam->eta = eta = (eigenval[3] - (eigenval[2]+eigenval[1])/2.0)/3.0;
	dparam->Amaj = dparam->eta / dparam->Dbar;
	dparam->Anu = pow(eta*(eta*eta - eps*eps), 1.0/3.0)/ dparam->Dbar;

	free_matrix (dten, 1, 3, 1, 3);
}

/* 
** inv_svd
** invert matrix a using svd algorithm
** based on numerical recipes book
*/

int
inv_svd (a, rows, cols, a_inv)
int rows, cols;
float	**a, **a_inv; 
{
	int i,j,k, flg=0;
	float	**u, *w, **v;	
	float	wmax, wmin;	

	u = matrix (1, rows, 1, cols);
	v = matrix (1, cols, 1, cols);
	w = vector (1, cols);

/* copy matrix to save it */
	for (i=1; i<=rows; i++)
		for (j=1; j<=cols; j++)
			u[i][j] = a[i][j];

/* do singular value decomposition */
	svdcmp (u, rows, cols, w, v); 

/* get max sigular value */
	wmax = 0.0;
	for (i=1; i<=cols; i++)
		if (w[i] > wmax) wmax = w[i];
	
/* set typical threshold and correct for smaller elements */
	wmin = wmax * 1.0e-6;
	for (i=1; i<=cols; i++) {
		if (w[i] > wmin) w[i] = 1.0/w[i];
		else { 
			w[i] = 0.0; 
			if (0) printf("inv_svd: set w-element %d to zero\n",i); 
			flg=1;
		}
	}

/* calculate the inverted matrix */
	for (i=1; i<=rows; i++) {
		for (j=1; j<=cols; j++) {
			a_inv[j][i] = 0.0;
			for (k=1; k<=cols; k++) {
				a_inv[j][i] += v[j][k] * w[k] * u[i][k];	
			}
		}
	}

	free_vector(w,1,cols);
	free_matrix(v,1,cols,1,cols);
	free_matrix(u,1,rows,1,cols);
	return(flg);
}

/* 
** inv_svd
** invert matrix a using svd algorithm
** based on numerical recipes book
*/

int
inv_svd_thresh (float **a, int rows, int cols, float **a_inv, float thresh)
{
	int i,j,k, flg=0;
	float	**u, *w, **v;	
	float	wmax, wmin;	

	u = matrix (1, rows, 1, cols);
	v = matrix (1, cols, 1, cols);
	w = vector (1, cols);

/* copy matrix to save it */
	for (i=1; i<=rows; i++)
		for (j=1; j<=cols; j++)
			u[i][j] = a[i][j];

/* do singular value decomposition */
	svdcmp (u, rows, cols, w, v); 

/* get max sigular value */
	wmax = 0.0;
	for (i=1; i<=cols; i++)
		if (w[i] > wmax) wmax = w[i];
	
/* set typical threshold and correct for smaller elements */
	wmin = wmax * thresh;
	for (i=1; i<=cols; i++) {
		if (w[i] > wmin) w[i] = 1.0/w[i];
		else { 
			w[i] = 0.0; 
			if (0) printf("inv_svd: set w-element %d to zero\n",i); 
			flg=1;
		}
	}

/* calculate the inverted matrix */
	for (i=1; i<=rows; i++) {
		for (j=1; j<=cols; j++) {
			a_inv[j][i] = 0.0;
			for (k=1; k<=cols; k++) {
				a_inv[j][i] += v[j][k] * w[k] * u[i][k];	
			}
		}
	}

	free_vector(w,1,cols);
	free_matrix(v,1,cols,1,cols);
	free_matrix(u,1,rows,1,cols);
	return(flg);
}

/* 
** atan2_custom
** atan2 that preserves quadrant 
*/
double 
atan2_custom (double y, double x)
{ 
	double res;

	if (x!=0.0 && y!=0.0) res = atan2 (y, x);
	else if (x == 0.0 && y != 0.0) {
    		if (y > 0.0) res =  atan2 (1.0, 0.0);
		if (y < 0.0) res = -atan2 (1.0, 0.0);
	}
	else if (x != 0.0 && y == 0.0) {
		if (x > 0.0) res = 0.0;
		if (x < 0.0) res = 2.0 * atan2 (1.0, 0.0);
	}
	else if (x == 0.0 && y == 0.0) res = 0.0;

	return (res);
}

/* 
** function: calc_xyz
** inputs: the experimental inputs
**	n - number of measurements, including I0
**	b[n] - array of b values for each measurement
** 	q[n][3] - array of vectors for each measurement
**	in[n] - array of measurement values
** output: 
**	dvec[7] - array of the diffusion tensor results
** 1. dxx 2. dyy 3. dzz 4. dbar 5. AW (anisotropy weighting)
** action: solve for the diffusion values in the case of 
** 	measurements made only in the orthogonal directions
** output: calculated diffusion values in dvec
** notes: 
** 1. assumes the data arrays start at 1
** 2. assumes the b values are in units of sec/mm2 divided by 1000.0
** 3. assumes the q vectors are normalized, and point only in x,y, and z
*/
void
calc_xyz (n, b, q, in, dvec)
int	n; 
float	*b, **q, *in, *dvec;
{
	int	i, j, k, ind;
	float   *ll, sx[4], sy[4], sxy[4], sx2[4], ssig2[4];

	ll = vector (1, n);

/* create the ll vector with the logarithm of the results */
	for (i = 1; i <= n; ++i) {
		if (in[i] > 0.0) ll[i] = log(in[i]);
		else ll[i] = 1.0e-8;
	}

/* initialize all the summation arrays */
	for (i = 1; i <= 3; ++i) {
		sx[i] = sy[i] = sxy[i] = sx2[i] = ssig2[i] = 0.0;
	}

/* loop through the data and increment the summation values */
/* set index for x,y,z and the b=0 contributes to all three */
	for (i = 1; i <= n; ++i) {
		for (j = 1; j <= 3; ++j) {
			if (b[i] != 0.0 && q[i][j] == 0.0) continue;

/* standard deviation = invesre of in[i] */
			ssig2[j] += in[i]*in[i];
			sx[j] += -b[i]*in[i]*in[i];
			sy[j] += ll[i]*in[i]*in[i];
			sxy[j] += -b[i]*ll[i]*in[i]*in[i];
			sx2[j] += b[i]*b[i]*in[i]*in[i];
		}
	}

/* solve the least square fit */
/* calc dx, dy, dz */
	for (i = 1; i <= 3; ++i) {
		if (ssig2[i]*sx2[i] - sx[i]*sx[i]==0.0) {
			dvec[i] = 0.0;
			printf("Error(caclc_xyz): division by zero in D%d\n",i);
		}
		else
		dvec[i] = (ssig2[i]*sxy[i] - sx[i]*sy[i])/(ssig2[i]*sx2[i] - sx[i]*sx[i]);
	}
/* calc dbar and AW orthogonal */
	dvec[4] = (dvec[1] + dvec[2] + dvec[3])/3.0;
	dvec[5] = sqrt( (dvec[1]-dvec[4])*(dvec[1]-dvec[4]) +
			(dvec[2]-dvec[4])*(dvec[2]-dvec[4]) +
			(dvec[3]-dvec[4])*(dvec[3]-dvec[4]) )/(dvec[4]*sqrt(6.0));

		free_vector (ll, 1, n);
}
