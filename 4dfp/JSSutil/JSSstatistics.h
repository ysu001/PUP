/*$Header: /home/usr/shimonyj/JSSutil/RCS/JSSstatistics.h,v 1.1 2008/12/30 03:37:52 avi Exp $*/
/*$Log: JSSstatistics.h,v $
 * Revision 1.1  2008/12/30  03:37:52  avi
 * Initial revision
 **/

void fit(float x[], float y[], int ndata, float sig[], int mwt, float *a,
	float *b, float *siga, float *sigb, float *chi2, float *q, float *lcc);

/* f-test using NR p.619 */
void ftest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *f, float *prob);

/* return ave and var of distribution NR p. 617 */
void avevar(float data[], unsigned long n, float *ave, float *var);

/* returns log of the gamma function NR p. 214 */
float gammln(float xx);

/* returns the incomplete beta function NR p. 227 */
float betai(float a, float b, float x);

/* evaluates continued fraction for betai NR p. 227 */
float betacf(float a, float b, float x);

/* paired t-test NR p. 618 small prob indicates significant difference in means */
void tptest(float data1[],float data2[],unsigned long n,float *t,float *prob);

/* chi square statistic NR p. 621 */
void chsone(float bins[],float ebins[],int nbins,int knstrn,float *df,float *chsq,float *prob);

void chstwo(float bins1[],float bins2[],int nbins,int knstrn,float *df,float *chsq,float *prob);

/* return the incomplete gamma function P and Q NR p. 218 */
float gammp(float a, float x);

float gammq(float a, float x);

void gser(float *gamser, float a, float x, float *gln);

void gcf(float *gammcf, float a, float x, float *gln);

/* error function routines, NR pg. 220 */
float erfunc(float x);

float erffc(float x);

float erfcc(float x);
