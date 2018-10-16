/*$Header: /data/petsun4/data1/src_solaris/librms/RCS/librms.h,v 1.3 2009/03/08 01:20:15 avi Exp $*/
/*$Log: librms.h,v $
 * Revision 1.3  2009/03/08  01:20:15  avi
 * dfft_(), drealt_(), dgeigen()
 *
 * Revision 1.2  2009/01/20  00:51:11  avi
 * dmatopr.f
 *
 * Revision 1.1  2007/09/18  23:10:40  avi
 * Initial revision
 **/

/************/
/* imgpad.f */
/************/
int	npad_	(int *n, int *m);
void	imgpad_ (float *imag, int *nx, int *ny, int *nz, float *imagp, int *nxp, int *nyp, int *nzp);
void	imgdap_ (float *imag, int *nx, int *ny, int *nz, float *imagp, int *nxp, int *nyp, int *nzp);

/*************/
/* gauss3d.f */
/*************/
void	gauss3d_ (float *imag, int *nx, int *ny, int *nz, float *cmppix, float *fhalf);

/**************/
/* cgauss3d.c */
/**************/
void	gauss3d  (float *imag, int *nx, int *ny, int *nz, float *cmppix, float *fhalf);

/***************/
/* img2lmask.c */
/***************/
void	img2lmask (int *pnx, int *pny, int *pnz, float *imag,   int *mask, float *mmppix, float *pfhalf, float *pcrit);

/***************/
/* imag2mask.c */
/***************/
void	imag2mask (int *pnx, int *pny, int *pnz, float *imag, short *mask, float *mmppix, float *pfhalf, float *pcrit);

/*************/
/* frmsmri.f */
/*************/
void	param2warp_ (int *mode, float *param, float *a);
void	warp2param_ (int *mode, float *a, float *param);

/***************/
/* param6opr.f */
/***************/
void	img2vrt_ (float *mmppix, float *center, float *vox2ras);
void	vrt2img_ (float *mmppix, float *center, float *ras2vox);
void	ang2rot_ (float *ang, float *rot);
void	rot2ang_ (float *rot, float *ang);

/************************/
/* fftsol.f or fftsun.f */
/************************/
void fft_   (float *a, float *b, int *nseg, int *n, int *nspn, int *isn);
void realt_ (float *a, float *b, int *nseg, int *n, int *nspn, int *isn);

/*************/
/* dfftsol.f */
/*************/
void dfft_   (double *a, double *b, int *nseg, int *n, int *nspn, int *isn);
void drealt_ (double *a, double *b, int *nseg, int *n, int *nspn, int *isn);

/************/
/* matopr.f */
/************/
void	matmul_	(float *a, float *b, float *c, int *n);
void	matinv_	(float *a, int *n, float *det);
void	matcop_ (float *a, float *b, int *n);
void	transpos_ (float *a, float *b, int *n);
void	geninv_	(float *a, float *g, float *b, int *n, float *det);
void	matmulg_ (float *a, int *l, int *m, float *b, int *n, float *c);

/*************/
/* dmatopr.f */
/*************/
void	dmatmul_   (double *a, double *b, double *c, int *n);
void	dmatinv_   (double *a, int *n, double *det);
void	dmatcop_   (double *a, double *b, int *n);
void	dtranspos_ (double *a, double *b, int *n);

/*************/
/* dmatinv.f */
/*************/
void	dmatinv_ (double *a, int *n, double *det);

/***********/
/* eigen.f */
/***********/
void	eigen_ (float *e, float *w, int *n);

/************/
/* deigen.f */
/************/
void	deigen_ (double *e, double *w, int *n);

/*************/
/* dgeigen.c */
/*************/
void	dgeigen (double *sig1, double *sig2, double *lambda, double *E, int *pn);

/**************/
/* determ12.f */
/**************/
float	determ12_ (float *a, float *det);

/************/
/* deconv.f */
/************/
void	conv_   (float *f, float *g, float *h, int *n);
void	deconv_ (float *h, float *g, float *f, int *n);
void	gconv_  (float *f, float *g, float *h, int *n, int *ldirection);
