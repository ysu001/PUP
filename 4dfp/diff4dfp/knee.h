/* constants */
#define NGRAY   4096      /* number of gray levels */

#define DY	256       /* y displacement = no. cols */
#define DY2     512       /* y displacement = no. cols */
#define DZ      65536     /* z displacement = no. 2D pixels */
#define M       1.5       /* power coefficient */
#define SPAN    15        /* size of low pass filter operation */
#define FZ_SCALE 10000    /* integer scaling */

/* constants for image overflow */
#define OV_ABS  0
#define OV_CUT  1

/* constants for fuzzy array use */
#define FF_NULL 0
#define FF_DISP 1
#define FF_AVG  2

/* constants for blob search */
#define AREA_MIN 64
#define HOLE_MIN 1
#define MAXBLOB  1024
#define MAXOBJ   32

/* constants for perim_trace */
#define CCW_FLAG 0
#define CW_FLAG  1
#define MAXPTS   8192	
#define MAXLINE  64
#define MAXDEFF  16

/* constants for contour map */
#define CON_NMAX 1024
#define CON_STEP 16
#define CON_AREA 32

#define MOD8(a) ((a)%8)
#define ABSV(a) (((a)>0)?((a)):(-(a)))
#define MAXV(a,b) (((a)>(b))?((a)):((b)))
#define MINV(a,b) (((a)<(b))?((a)):((b)))
#define RTOD(a) ((a)*(180.0/M_PI))
#define DTOR(a) ((a)*(M_PI/180.0))
#define IND(x,y,z) ((x)+(y)*DY+(z)*DZ)

/* structures */

/* convex defficiency struct */
typedef struct cdeff_data {
    int area, len, circ;
    int n1,x1,y1,n2,x2,y2;
    double ang;
} CDEFF;

typedef struct line_data {
    int indx,len;
    double ang,var;
    int x1,y1,x2,y2;
} FITLINE;
    
typedef struct cont_data {
    int xbeg, ybeg, hvbeg;
    int z, area, perim;
    int parent, nson, depth;
} CONT;

#define POINT struct point_data
typedef POINT {
    int x, y, d, i;
    POINT *next;
} *PPOINT;

#define BLOB struct blob_data
typedef BLOB {
    int    xbeg, ybeg, zbeg;
    int    nblob, parent;
    int    m000, m100, m010, m001;
    double m110, m101, m011, m200, m020, m002;
    int    nholes, nbays;
    BLOB   *hole, *bay, *next;
    PPOINT perim, con_hull;
} *PBLOB;

typedef struct moment_dat {
    int    xbeg, ybeg, zbeg, fgflg;
    int    m000, m100, m010, m001;
    double m110, m101, m011, m200, m020, m002;
} MOMENT;

/* global arrays */
PBLOB   objects[MAXOBJ];
CONT    cont[CON_NMAX];
CDEFF   cdeff[MAXDEFF];
FITLINE fitline[MAXLINE];

/* knee.c */
void	latknee_pros(short int *orig,short int *grad,
	int, int, int ncol,int npix, int *res);
void	edge_clr(short int *orig,int orig_grad,int ncol,int nrow);
void	latknee_femur(short int *grad,int seg_grad,int ncol,int npix,
		int x1,int y1,int x2,int y2,double ang);
void	latknee_tibfib(short int *grad,int seg_grad,int ncol,int npix,
		int x1,int y1,int x2,int y2,double ang);
int	num_femur_bay(int);

/* morph.c - utilities */
void    get_image(unsigned char *image,char *filename,int size);
void    put_image(unsigned char *image,char *filename,int size);
int	dicom_hdr(short int *src,int nhdr,int *ncol,int *nrow,int *npix);
void	inv_chk_img(short int *src,int npix);
void    compress4(short int *src,int *ncol,int *npix);
void    vax_to_uchar(), int_to_char();
void    prt_hex_2D();

/* contour.c */
int     con_map();
void    tree_prt();
void    con_trace(), con_fill();
void    con_parents(), con_depth();
int     trace_flag(), parent_flag();

/* fuzzy.c */
void	seg_fuzz(short int *src,int npart,int ipart,double eps,
	double *seed,int npix);
void    defuzz();
void    avgpix(), avgfuzz();

void	sector_mask(short int *src, short int *dst,int x1,int y1,
	double dir,double wid,int ncol,int npix);
void	sector_fill(short int *dst, double rho, double ang,
	int id, int ncol, int npix);
double  ang_norm(int x1,int y1,int x2,int y2);
double  ang_slope(int x1,int y1,int x2,int y2);
double  rho_norm(int x1,int y1,double *ang);
void	hough_line(short int *src, int ncol, int nrow);
void	hough_circ(short int *mask,int *xp,int *yp,int npts,int ncol,int npix);
double  norm_ang(double ang);
double  norm_angp(double ang);

/* morph.c */
void    conv3x3(short int *src,short int *dst,int *conv,int norm,
		int ncol, int npix);
void    conv3x3x3(short int *src,short int *dst,int *conv,int norm,
		int ncol, int npix, int nvol);
void    sobelx(short int *src,short int *dst, int ncol,int npix);
void    sobely(short int *src,short int *dst, int ncol,int npix);
void    mag_img(short int *srcx,short int *srcy,short int *dst,int npix);
void    linear_img(short int *src0,short int *src1,short int *dst,
	int c0, int c1, int c2, int ov_flag, int npix);
void    scale_img(short int *src, short int *dst,
	int c0, int c1, int ov_flag, int npix);
void    grad_img(short int *src,short int *dst,int scale,int ncol,int npix);
void    segment_img(short int *src,short int *dst,int lo,int hi,int npix);
void    segment_img_bin(short int *src,short int *dst,int n,int npix);
void    and_img(short int *src0,short int *src1,short int *dst,int npix);
void	or_img(), xor_img();
void	not_img(short int *src, short int *dst, int npix);
void    dilate_erode(short int *src,short int *dst,int n,int ncol,int npix);
void    dilate_erode_bin(short int *src,short int *dst,int n,int ncol,int npix);
void    dilate_erode_3D();
void    open_img(short int *src,short int *dst,int n,int ncol,int npix);
void    close_img(short int *src,short int *dst,int n,int ncol,int npix);
void    profile_img();
void 	hist_img(short int *img, int npix), draw_hist(int *bin,int max);

/* feat.c */
int 	find_blob(short int *src, short int *dst, int ncol, int npix);
int 	find_blob_rv(short int *src, short int *dst, int ncol, int npix);
int 	find_blob_2d(short int *src, short int *dst, int ncol, int npix);
int 	find_blob_3d(short int *src, short int *dst, int ncol, int npix);
int     calc_blobs(int nblob, int *lub);
int     calc_blobs_fill(int nblob, int *lub);
int	init_moments(short int *bimg,int i,int nblob,
	int nprnt,int fgflg,int ncol);
void    add_moments(short int *bimg, int i, int nblob, int ncol), connect_reg(int n, int m);
PBLOB   alloc_blob();

int	perim_trace(short int *src, short int *dst, int xbeg, int ybeg,
	int *xp,int *yp,int id, int cw_flg, int ncol,int npix);
int	convex_hull(short int *dst,int *xp,int *yp,int npts,int ncol,int npix);
void	trace_line(int x1,int y1,int x2,int y2,int *xo,int *yo);
void	draw_line(short int *dst,int x1,int y1,int x2,int y2,int ncol);
void	draw_circle(short int *dst,int xcen,int ycen,double rad,
	int ncol,int npix);
double  fit_line(int *xp,int *yp,int npts);
int	find_line(short int *dst,int *xp,int *yp,int npts,
	    int minlen,double maxvar,int ncol,int npix);

void    amoeba();
double  amotry();
