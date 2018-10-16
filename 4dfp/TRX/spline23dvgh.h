/*$Header: /data/petsun4/data1/src_solaris/TRX/RCS/spline23dvgh.h,v 1.2 2007/04/22 03:34:28 avi Exp $*/
/*$Log: spline23dvgh.h,v $
 * Revision 1.2  2007/04/22  03:34:28  avi
 * eliminate splint2dvgh_testf_ ()
 *
 * Revision 1.1  2007/04/01  00:04:01  avi
 * Initial revision
 **/

/*****************/
/* spline2dvgh.f */
/*****************/
extern void	spline2dvgh_rcs_ (void);
extern void	splint2dv_test_	(void);
extern void	splint2dvgh_testh_ (void);
extern void	splint2dvgh_testg_ (void);
extern void	splint2dv_	(float *imag, int *nx, int *ny, float *d2xi, float *d2yi, float *x, float *y,
					float *v);
extern void	splint2dvg_	(float *imag, int *nx, int *ny, float *d2xi, float *d2yi, float *x, float *y,
					float *v, float *grad);
extern void	splint2dvgh_	(float *imag, int *nx, int *ny, float *d2xi, float *d2yi, float *x, float *y,
					float *v, float *grad, float *hess);

/*****************/
/* spline3dvgh.f */
/*****************/
extern void spline3dvgh_rcs_ (void);
extern void splint3dvl_test_ (void);
extern void splint3dv_test_ (void);
extern void splint3dvgh_test_h (void);
extern void splint3dvgh_test_g (void);
extern void splinex_	(float *imag, int *nx, int *ny, int *nz, float *d2xi);
extern void spliney_	(float *imag, int *nx, int *ny, int *nz, float *d2yi);
extern void splinezf_	(float *imag, int *nx, int *ny, int *nz, float *d2zi);
extern void splineza_	(float *imag, int *nx, int *ny, int *nz, float *d2zi);
extern void splint3dv_	(float *imag, int *nx, int *ny, int *nz, float *d2xi, float *d2yi, float *d2zi,
				float *x, float *y, float *z, float *v);
extern void splint3dvg_ (float *imag, int *nx, int *ny, int *nz, float *d2xi, float *d2yi, float *d2zi,
				float *x, float *y, float *z, float *v, float *grad);
extern void splint3dvl_ (float *imag, int *nx, int *ny, int *nz, float *d2xi, float *d2yi, float *d2zi,
				float *x, float *y, float *z, float *v, float *del2);
extern void splint3dvgh_(float *imag, int *nx, int *ny, int *nz, float *d2xi, float *d2yi, float *d2zi,
				float *x, float *y, float *z, float *v, float *grad, float *hess);
extern void splint3dvgl_(float *imag, int *nx, int *ny, int *nz, float *d2xi, float *d2yi, float *d2zi,
				float *x, float *y, float *z, float *v, float *grad, float *del2);
