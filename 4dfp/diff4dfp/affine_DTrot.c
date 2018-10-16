/*$Header: /home/usr/shimonyj/diff4dfp/RCS/affine_DTrot.c,v 1.2 2008/09/16 02:59:29 avi Exp $*/
/*$Log: affine_DTrot.c,v $
 * Revision 1.2  2008/09/16  02:59:29  avi
 * #include <JSSutil.h>
 *
 * Revision 1.1  2005/12/08  03:16:37  avi
 * Initial revision
 **/
static char rcsid[] = "$Id: affine_DTrot.c,v 1.2 2008/09/16 02:59:29 avi Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <JSSutil.h>

/*
Vnorm
return normalized vector a[1:n]
*/
void Vnorm(float *a, int n)
{
	int i;
	float sum=0.0;
	
	for (i=1; i<=n; i++) sum += a[i]*a[i];
	sum = sqrt(sum);
	for (i=1; i<=n; i++) a[i] /= sum;
}

/*
Vdot
return dot product of 2 vectors a[1:n], b[1:n]
*/
float Vdot(float *a, float *b, int n)
{
	int i;
	float sum=0.0;
	
	for (i=1; i<=n; i++) sum += a[i]*b[i];
	return(sum);
}

/*
Vcross
return cross product of 2 vectors a[1:3], b[1:3] in c[1:3]
*/
void Vcross(float *a, float *b, float *c)
{
	c[1] = a[2]*b[3] - a[3]*b[2];
	c[2] = a[3]*b[1] - a[1]*b[3];
	c[3] = a[1]*b[2] - a[2]*b[1];
}

/*
MVmul
Matrix multiply by vector a*b=c with low index = 1  
a[nra][ncarb], b[ncarb], c[nra], as Mmul with ncb=1 
*/
void MVmul(float **a, float *b, int nra, int ncarb, float *c)
{
	int i,k;
	
	for (i=1; i<=nra; i++) {
		c[i] = 0.0;
		for (k=1; k<=ncarb; k++)
			c[i] += a[i][k]*b[k];
	}
}

/*
rotate_tensor
rotate tensor d using matrix rot 
input: 3x3 matrix rot, 3x3 matrix d
output: 3x3 tensor rdrt
all indecis from 1 to 3
*/
void rotate_tensor(float **rot, float **d, float **rdrt)
{
	int k,l,m,n;

	for (m=1; m<=3; m++)
	for (n=1; n<=3; n++) {
		rdrt[m][n]=0.0;
		for (k=1; k<=3; k++)
		for (l=1; l<=3; l++) {
			rdrt[m][n] += rot[k][m]*d[k][l]*rot[l][n];
		}
	}
}

/*
rotate_matrix
rotate matrix d using matrix rot 
input: 3x3 matrix rot, 3x3 matrix d
output: 3x3 matrix rotd
all indecis from 1 to 3
*/
void rotate_matrix(float **rot, float **d, float **rotd)
{
	int k,l,m,n;

	for (m=1; m<=3; m++)
	for (n=1; n<=3; n++) {
		rotd[m][n]=0.0;
		for (k=1; k<=3; k++) {
			rotd[m][n] += rot[m][k]*d[k][n];
		}
	}
}

/* 
make_rotmat
create rotation matrix around axis r by given angle 
input: normed axis of rotation r (index 1-3)
	cosine of rotation angle cang
output: rotation matrix rot (index 1-3)
*/
void make_rotmat(float *r, float cang, float **rot)
{
	int k,l,m,n;
	float sa,ca,sb,cb;
	float **Rz, **Rr;
	
	Rz = matrix(1,3,1,3);
	Rr = matrix(1,3,1,3);
	
	sa = r[2]/sqrt(r[1]*r[1]+r[2]*r[2]); 
	ca = r[1]/sqrt(r[1]*r[1]+r[2]*r[2]);
	sb = sqrt(r[1]*r[1]+r[2]*r[2]); 
	cb = r[3];
	
	/* form rotation matrix Arfken p.180 */
	Rr[1][1] = cb*ca;
	Rr[1][2] = cb*sa;
	Rr[1][3] = -sb;
	Rr[2][1] = -sa;
	Rr[2][2] = ca;
	Rr[2][3] = 0.0;
	Rr[3][1] = sb*ca;
	Rr[3][2] = sb*sa;
	Rr[3][3] = cb;
	
	if (cang > 1.0 || cang < -1.0) cang = 1.0;
	/* init the inverse diagonal matrix */
	Rz[1][1] = cang;
	Rz[1][2] = -sqrt(1.0 - cang*cang);
	Rz[1][3] = 0.0;
	Rz[2][1] = sqrt(1.0 - cang*cang);
	Rz[2][2] = cang;
	Rz[2][3] = 0.0;
	Rz[3][1] = 0.0;
	Rz[3][2] = 0.0;
	Rz[3][3] = 1.0;

	
	/* rotate the tensor into the lab frame */
	for (m=1; m<=3; m++)
	for (n=1; n<=3; n++) {
		rot[m][n]=0.0;
		for (k=1; k<=3; k++)
		for (l=1; l<=3; l++) {
			rot[m][n] += Rr[k][m]*Rz[k][l]*Rr[l][n];
		}
	}
	
	free_matrix(Rz,1,3,1,3);
	free_matrix(Rr,1,3,1,3);
}

/* 
affine_DTrot
correct the eigenvalue matrix for affine rotation
algorithm from Alexander IEEE TMI 20:1131-1139, 2001
input: affine matrix atm (index 1-3)
	eignevector matrix evec (index 1-3)
output: corrected eigenvector matrix cevec (index 1-3)
*/
void affine_DTrot(float **atm, float **evec, float **cevec)
{
	int i;
	float	*e1,*e2,*n1,*n2,*raxis,*Pn2,*N1,*N2,costh;
	float	**R1, **R2, **RT;

	e1 = vector(1,3);
	e2 = vector(1,3);
	n1 = vector(1,3);
	n2 = vector(1,3);
	raxis = vector(1,3);
	Pn2 = vector(1,3);
	N1 = vector(1,3);
	N2 = vector(1,3);	

	R1 = matrix(1,3,1,3);
	R2 = matrix(1,3,1,3);
	RT = matrix(1,3,1,3);	

	/* get major and mid evector, store in e1 and e2 */
	for (i=1; i<=3; i++) {
		e1[i] = evec[i][3];
		e2[i] = evec[i][2];
	}
	
	MVmul(atm, e1, 3, 3, n1);
	MVmul(atm, e2, 3, 3, n2);
	Vnorm(n1, 3);
	Vnorm(n2, 3);
	
	costh = Vdot(e1, n1, 3);
	Vcross(e1, n1, raxis);
	Vnorm(raxis, 3);
	
	make_rotmat(raxis, costh, R1);
	
	costh = Vdot(n1, n2, 3);
	for (i=1; i<=3; i++) Pn2[i] = n2[i] - costh*n1[i];
	Vnorm(Pn2, 3);
	
	/* testing: N1 = n1 */
	MVmul(R1, e1, 3, 3, N1);
	MVmul(R1, e2, 3, 3, N2);
	Vnorm(N1, 3);
	Vnorm(N2, 3);
	
	costh = Vdot(N2, Pn2, 3);
	make_rotmat(N1, costh, R2);
	
	rotate_matrix(R1, evec, RT);
	rotate_matrix(R2, RT, cevec);
	
	/* deallocate memory */
	free_vector(e1,1,3);
	free_vector(e2,1,3);
	free_vector(n1,1,3);
	free_vector(n2,1,3);
	free_vector(raxis,1,3);
	free_vector(Pn2,1,3);
	free_vector(N1,1,3);
	free_vector(N2,1,3);	

	free_matrix(R1,1,3,1,3);
	free_matrix(R2,1,3,1,3);
	free_matrix(RT,1,3,1,3);
	
	return;
}

