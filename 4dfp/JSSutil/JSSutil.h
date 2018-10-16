/*$Header: /home/usr/shimonyj/JSSutil/RCS/JSSutil.h,v 1.3 2008/11/10 23:41:52 shimonyj Exp $*/
/*$Log: JSSutil.h,v $
 * Revision 1.3  2008/11/10  23:41:52  shimonyj
 * added gaussj_flg
 *
 * Revision 1.2  2007/08/30  04:57:52  avi
 * int inv_svd_thresh()
 *
 * Revision 1.1  2007/08/28  03:33:30  avi
 * Initial revision
 **/

/* definitions */
#define PI 	3.14159265358979323846
#define TWOPI	6.28318530717958647692

static float sqarg;
#define SQR(a) ((sqarg=(a)) == 0.0 ? 0.0 : sqarg*sqarg)

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1)>(maxarg2) ? (maxarg1) : (maxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1)<(iminarg2) ? (iminarg1) : (iminarg2))

/* structures */
typedef struct FCOMPLEX {float r,i;} fcomplex;
/* structure diffusion data */
typedef struct {
        float	ang[4];
        float	Dbar;
        float	Asig;
        float	Anu;
        float	prolat;
        float	eta;
        float	eps;
        float	Amaj;
        float	Amin;
} D_PARAMS;

/* nrutil */
void byte_swap(char *in);


int *ivector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
int **imatrix (long nrl, long nrh, long ncl ,long nch);
float *vector(long nl, long nh);
float **matrix (long nrl, long nrh, long ncl ,long nch);
double *dvector(long nl, long nh);
double **dmatrix (long nrl, long nrh, long ncl ,long nch);
fcomplex *Cvector(long nl, long nh);
fcomplex **Cmatrix (long nrl, long nrh, long ncl ,long nch);

void free_ivector(int *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl ,long nch);
void free_vector(float *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl ,long nch);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl ,long nch);
void free_Cvector(fcomplex *v, long nl, long nh);
void free_Cmatrix(fcomplex **m, long nrl, long nrh, long ncl ,long nch);

/* srandom */
float ran1(long *idum);
float expdev(long *idum);
float gasdev(long *idum);

/* lin_algebra */
void gaussj(float **a, int n, float **b, int m);
int  gaussj_flg(float **a, int n, float **b, int m);
void covsrt(float **covar, int ma, int ia[], int mfit);
void jacobi(float **a,int n,float d[],float **v,int *nrot);
void eigsrt(float d[], float **v, int n);
void svdcmp(float **a, int m, int n, float w[], float **v);

/* dtensor */
void set_svd (int n, float *b,float **q,float **aa);
void calc_svd (int n,float *in,float **aa,float *dvec);
void eigen_calc(float *dvec,float **dinv,float *eigenval,float **eigenvec,float *ang);
void eigen_calc_noswap(float *dvec,float **dinv,float *eigenval,float **eigenvec,float *ang);
void dparam_calc(float *dvec, float **dinv, float *eigenval, float **eigenvec, D_PARAMS *dparam);
int inv_svd		(float **a,int rows,int cols,float **a_inv);
int inv_svd_thresh	(float **a,int rows,int cols,float **a_inv, float thresh);
int inv_svd_filter	(float **a,int rows,int cols,float **a_inv, float thresh, int nkeep);
double atan2_custom (double y, double x);
void calc_xyz (int n,float *b,float **q,float *in,float *dvec);

