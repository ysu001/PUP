/*$Header: /home/usr/shimonyj/supercomp/RCS/JSSnrutil.c,v 1.2 2007/08/28 03:51:33 avi Exp $*/
/*$Log: JSSnrutil.c,v $
 * Revision 1.2  2007/08/28  03:51:33  avi
 * remove #include <omp.h> and inconsistently defined err[rwm]
 *
 * Revision 1.1  2007/08/28  03:33:22  avi
 * Initial revision
 **/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <JSSutil.h>

#define NR_END 1
#define FREE_ARG char*

/* byte swap routines */
void byte_swap(char *in)
{
	char temp;

	temp=in[0]; in[0]=in[3]; in[3]=temp;
	temp=in[1]; in[1]=in[2]; in[2]=temp;
}

/* error handlers */
void nrerror(char error_text[]) 
{
	fprintf(stderr,"%s\n", error_text);
	exit(1);
}

/* allocate and free an int vector range nl..nh */
int *ivector(long nl, long nh)
{
	int *v;
	
	v = (int *)malloc( ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("Error(ivector): allocation failure\n");
	return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

/* allocate and free an unsigned long vector range nl..nh */
unsigned long *lvector(long nl, long nh)
{
	unsigned long *v;
	
	v = (unsigned long *)malloc( ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("Error(lvector): allocation failure\n");
	return v-nl+NR_END;
}

void free_lvector(unsigned long *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

/* allocate and free a float vector range nl..nh */
float *vector(long nl, long nh)
{
	float *v;
	
	v = (float *)malloc( ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("Error(vector): allocation failure\n");
	return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

/* allocate and free a double vector range nl..nh */
double *dvector(long nl, long nh)
{
	double *v;
	
	v = (double *)malloc( ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("Error(dvector): allocation failure\n");
	return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

/* allocate and free a complex vector range nl..nh */
fcomplex *Cvector(long nl, long nh)
{
	fcomplex *v;
	
	v = (fcomplex *)malloc( ((nh-nl+1+NR_END)*sizeof(fcomplex)));
	if (!v) nrerror("Error(Cvector): allocation failure\n");
	return v-nl+NR_END;
}

void free_Cvector(fcomplex *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

/* allocate and free a int matrix with range m[nrl..nrh][ncl..nch] */
int **imatrix(long nrl, long nrh, long ncl ,long nch)
{
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	int **m;
	
	/* allocate pointers to rows */
	m = (int **) malloc( ((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("Error(imatrix): allocation failure\n");
	m += NR_END;
	m -= nrl;
	
	/* allocate rows and pointers to them */
	m[nrl] = (int *) malloc( ((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("Error(imatrix): allocation failure\n");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl ,long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


/* allocate and free a float matrix with range m[nrl..nrh][ncl..nch] */
float **matrix(long nrl, long nrh, long ncl ,long nch)
{
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	float **m;
	
	/* allocate pointers to rows */
	m = (float **) malloc( ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("Error(matrix): allocation failure\n");
	m += NR_END;
	m -= nrl;
	
	/* allocate rows and pointers to them */
	m[nrl] = (float *) malloc( ((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("Error(matrix): allocation failure\n");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

void free_matrix(float **m, long nrl, long nrh, long ncl ,long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* allocate and free a double matrix with range m[nrl..nrh][ncl..nch] */
double **dmatrix(long nrl, long nrh, long ncl ,long nch)
{
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	double **m;
	
	/* allocate pointers to rows */
	m = (double **) malloc( ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("Error(dmatrix): allocation failure\n");
	m += NR_END;
	m -= nrl;
	
	/* allocate rows and pointers to them */
	m[nrl] = (double *) malloc( ((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("Error(dmatrix): allocation failure\n");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl ,long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* allocate and free a fcomplex matrix with range m[nrl..nrh][ncl..nch] */
fcomplex **Cmatrix(long nrl, long nrh, long ncl ,long nch)
{
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	fcomplex **m;
	
	/* allocate pointers to rows */
	m = (fcomplex **) malloc( ((nrow+NR_END)*sizeof(fcomplex*)));
	if (!m) nrerror("Error(Cmatrix): allocation failure\n");
	m += NR_END;
	m -= nrl;
	
	/* allocate rows and pointers to them */
	m[nrl] = (fcomplex *) malloc( ((nrow*ncol+NR_END)*sizeof(fcomplex)));
	if (!m[nrl]) nrerror("Error(Cmatrix): allocation failure\n");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

void free_Cmatrix(fcomplex **m, long nrl, long nrh, long ncl ,long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

