/*$Header: /data/petsun4/data1/src_solaris/nonlin/RCS/morphc_4dfp.c,v 1.15 2007/07/04 03:24:10 avi Exp $*/
/*$Log: morphc_4dfp.c,v $
 * Revision 1.15  2007/07/04  03:24:10  avi
 * default LAMBDA 1.0 -> 2.0
 * comment out f_init() and f_exit() for gcc
 *
 * Revision 1.14  2007/07/03  05:10:32  avi
 * endian compliant i/o
 *
 * Revision 1.13  2004/11/05  22:14:07  rsachs
 * Fixed the margin computation code.
 *
 * Revision 1.12  2004/10/31  06:45:40  avi
 * call to get_date_log repalce with call to get_time_usr (eliminate -lmri dependence)
 * numerous errors corrected
 *
 * Revision 1.11  2004/10/29  21:31:19  rsachs
 * volume labels in main output rec file
 *
 * Revision 1.10  2004/10/29  19:42:58  rsachs
 * Added voxdim, print_iterdat.
 *
 * Revision 1.8  2004/10/21  00:33:15  rsachs
 * Some improvements. Cosmetics.
 *
 * Revision 1.6  2004/10/19  21:20:18  rsachs
 * Once more.
 *
 * Revision 1.5  2004/10/19  20:50:54  rsachs
 * Replaced 'Get4dfpDimN' with 'get_4dfp_dimo'.
 *
 * Revision 1.4  2004/10/19  20:36:48  rsachs
 * Eliminated some volumes, to save memory.
 *
 * Revision 1.3  2004/10/19  19:15:30  rsachs
 * Removed 'writeifh'. Added 'setprog'. Minor cosmetics.
 *
 * Revision 1.2  2004/10/14  23:30:54  rsachs
 * Initial Revision.
 **/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <Getifh.h>
#include <rec.h>
#include <endianio.h>

#define NVOL	15		/* number of output volumes in main 4dfp output */
#define MAXL	256
#define NIMG	3
#define LAMBDA	2.0
#define MU	2.0
#define MITR	64		/* default number of iterations */
#define DELINT  32		/* default movie interval */
#define T       0.06
#define FWHMP   5.0
#define FWHMD   1.0

/**********************************/
/* type and structure definitions */
/**********************************/
typedef float typ;

typedef struct {
	typ factor, crit, mu, lambda, dt, fwhmd, fwhmp;
	int ii, id, ip;
} QDAT;

typedef struct {
	typ bf[3], vel[3], sol[3], disp[3];
} ITERDAT;

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: morphc_4dfp.c,v 1.15 2007/07/04 03:24:10 avi Exp $";
static char	program[MAXL];

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void read_file_float (char *filename, float *stack, int dimension, int isbig) {
	FILE *fp;
 
	if (!(fp = fopen (filename, "rb"))
	|| gread ((char *) stack, sizeof (float), dimension, fp, isbig)
	|| fclose (fp)) errr (program, filename);
}

void getrange (char *string, float *minval, float *maxval) {
              char    *str;
    
              str = strstr (string, "to");
              if (str) {
                      *str = '\0';
                      *minval = atof (string);
                      *maxval = atof (str + 2);
              } else {
                      *minval = 0.0;
                      *maxval = atof (string);
              }
}

void print_date_log (char *str) {
	char dat[64];

	get_time_usr (dat);
	printf ("%s %s\n", str, dat);
}

/************/
/* imgpad.f */
/************/
extern int	npad_	(long *n, long *m);
extern void	imgpad_ (float *imag, long *nx, long *ny, long *nz, float *imagp, long *nxp, long *nyp, long *nzp);
extern void	imgdap_ (float *imag,  int *nx,  int *ny,  int *nz, float *imagp,  int *nxp,  int *nyp,  int *nzp);
extern void	delu_ (typ *v, typ *u, typ *du, typ *vn, long *nx, 
			long *ny, long *nz, float *voxdim, float *dt); 
extern void	ntarg_ (typ *u, typ *img1, typ *img2, float *voxdim, long *nx,
			long *ny, long *nz);
extern void	copyimg_ (typ *img0, typ *img1, long *nx, long *ny, long *nz); 
extern void	mskpac_ (short *mask, long *nx, long *ny, long *nz, short *maskp, long *nxp, long *nyp, long *nzp); 
extern void	parzen_ (typ *img, long *nx, long *ny, long *nz, float *voxdim, typ *fwhm); 
extern void	norms_ (typ *v, typ *vn, int *n); 
extern void	b2v4_ (float *b, long *nx, long *ny, long *nz, 
			float *voxdim, float *mu, float *lambda, float *v);
extern void	chkvel_ (typ *v, typ *b, typ *err, typ *vn, float *voxdim, long *nx,
			long *ny, long *nz, typ *lambda, typ *mu);
extern void	norms (typ *img, typ *nr, short *mask, int size, char *str);
extern void	inorms (short *mask, long size, char *str);
extern void	piic (typ *img0, typ *img1, typ *img2, typ *imgf, typ *imgh, typ *imgd,
			typ *bb, short *mask0, short *mask1, typ *voxdim, int nx, int ny,
			int nz, int step, typ *factor, typ *crit, typ fwhm);  
extern void	mbodyforce_ (float *img1, short *mask1, long *nx, long *ny,
			long *nz, float *img2, short *mask2, float *voxdim, float *b);
extern typ	imagerr (typ *img0, typ *img1, short *mask0, short *mask1, int nx, int ny, int nz, char *str);
extern void	dvrgnc_ (typ *v, typ *div, typ *voxdim, long *nx, long *ny, long *nz);

void outimg (typ *img, int nx, int ny, int nz, typ *imgp, int nxp, int nyp, int nzp, FILE *fp, int dim, char control) {
/* note inconsistent casting of nx, ny, nz, nxp, nyp, nzp */
	imgdap_ (img, &nx, &ny, &nz, imgp, &nxp, &nyp, &nzp);
	gwrite ((char *) img, sizeof(typ), dim, fp, control);
}  

void print_iterdat (char *program, char *imgroot, ITERDAT *iterdat, int nr) {
	ITERDAT	 *D;
	FILE     *fp;
	char      datfile[MAXL];
	int       j;

	strcpy (datfile, imgroot);
	strcat (datfile,".dat"); 
	if (!(fp = fopen(datfile, "w"))) errw (program, datfile); 
	fprintf(fp, "#    Bodyforce                     Velocity                    Solution                      Displacement\n");
	for (j = 0; j < nr; j++) {
		D = iterdat + j;
		fprintf(fp, "%3d  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e   %8.2e   %8.2e\n", j,
			D->bf[0], D->bf[1], D->bf[2], D->vel[0], D->vel[1], D->vel[2],
			D->sol[0],D->sol[1],D->sol[2],D->disp[0],D->disp[1],D->disp[2]);
	}
}

/**********************************************/
/* qstep performs one linear integration step */
/**********************************************/
void qstep (typ **img, typ **arrptr, short **maskp, float *voxdim,   
            long *imgpdim, QDAT *qdat, ITERDAT *iterdat) {
	int	ii, id, ip, j, k, step = 1, m;
	long	int nxp, nyp, nzp, size, jj;
	typ	*img1, *img2, *img3, *img4, *img5, *imgd, *imgf, *imgh;
	typ	*u, *bb, *b, *v, *er, *du, fwhmd, fwhmp, factor, crit;
	typ	nr[3], vn[3], mu, lambda, dt;

	fwhmd = qdat-> fwhmd; fwhmp = qdat-> fwhmp;
	img1 = img[0]; img2 = img[1]; img3 = img[2]; img4 = img[3]; 
	img5 = img[4]; imgd = img[5]; imgf = img[6]; imgh = img[7]; 

	u = arrptr[0]; bb = arrptr[1]; b = arrptr[2]; 
	v = arrptr[5]; er = arrptr[6]; du = arrptr[7];
	factor = qdat->factor; crit = qdat->crit; mu = qdat->mu;
	lambda = qdat->lambda; dt   = qdat->dt;   ii = qdat->ii;
	id     = qdat->id;     ip   = qdat->ip; 

	nxp = imgpdim[0]; nyp = imgpdim[1]; nzp = imgpdim[2]; size = nxp * nyp * nzp; m = 3*size;

	printf ("\nITER #  %d \n", ii);
	print_date_log ("\nBeginning the iteration... ");
	ntarg_ (u, img5, img4, voxdim, &nxp, &nyp, &nzp); 
	print_date_log ("\nAfter  ntarg()... ");
	if (!(ii%ip)) {
		piic (img1, img4, img3, imgf, imgh, imgd, bb, maskp[0], maskp[1],
                      voxdim, nxp, nyp, nzp, step, &factor, &crit, fwhmp);
	} else {
		for (j = 0; j < size; j++) img3[j] = img4[j] * imgh[j] * imgf[j];
	}
	imagerr (img3,img4,maskp[0],maskp[1],nxp,nyp,nzp,"img3 & img4 ");
	norms (img4, nr, maskp[0], size, "Following ntarg_ Norm(img4): ");
	norms (img3, nr, maskp[0], size, "Following ntarg_ Norm(img3): ");
	mbodyforce_ (img1,maskp[0],&nxp,&nyp,&nzp,img3,maskp[1],voxdim,b);
	print_date_log ("\nFollowing mbodyforce...");
	printf ("qstep: B-vector:  "); norms_ (b,  iterdat[ii].bf, &m);
	print_date_log ("\nBefore b2v4_ ...");
        b2v4_ (b, &nxp, &nyp, &nzp, voxdim, &mu, &lambda, v);
 	print_date_log ("\nBefore chkvel_ ");	
       	chkvel_ (v,b,er,iterdat[ii].sol,voxdim,&nxp,&nyp,&nzp,&lambda,&mu);
       	printf ("qstep: Velocity vector: "); 
	norms_ (v, iterdat[ii].vel, &m);
	norms (du, nr, maskp[0], size, "du (Before delu_): ");
        print_date_log ("\nBefore delu()...");
	delu_ (v, u, du, iterdat[ii].disp, &nxp, &nyp, &nzp, voxdim, &dt);
/************************************************/
/* smooth displacement field every id iterations*/
/************************************************/
        if (!(ii%id)) for (k = 0; k < 3; k++) parzen_ (u + k*size, &nxp, &nyp, &nzp, voxdim, &fwhmd); 
}

void qstep2 (typ **img, typ **arrptr, short **maskp, float *voxdim, 
             long *imgpdim, QDAT *qdat, ITERDAT *iterdat) {
	static typ	*arrptr2[8], *imga[8];
	static typ	*ba, *va, *ua, *dua, vn[3];
	static typ	*img3, *img4;
	static typ	*b, *v, *u, *du, dt, th;
	static long	nxp, nyp, nzp, m, j;

	nxp = imgpdim[0]; nyp = imgpdim[1]; nzp = imgpdim[2]; 
	m = 3 * nxp * nyp * nzp;

	for (j = 0; j < 8; j++) {
 		arrptr2[j] = arrptr[j];
		imga[j]    = img[j];
 	}
	
	if (!(ba   = (typ *) calloc (m, sizeof(typ)))
       	||  !(dua  = (typ *) calloc (m, sizeof(typ)))
       	||  !(ua   = (typ *) calloc (m, sizeof(typ)))
       	||  !(va   = (typ *) calloc (m, sizeof(typ)))) errm ("qstep2");

	img3 = img[2]; img4 = img[3];
	u  = arrptr[0]; b = arrptr[2]; 
	v = arrptr[5]; du = arrptr[7]; 
 
	arrptr2[0] = ua;  arrptr2[2] = ba; 
	arrptr2[5] = va; arrptr2[7] = dua;
 
	for (j = 0; j < m; j++) {
		ua[j] = u[j]; va[j] = v[j]; dua[j] = du[j]; ba[j] = b[j];
 	}

	dt = qdat->dt; th = dt/2.; qdat->dt = th;
	qstep (imga, arrptr2, maskp, voxdim, imgpdim, qdat, iterdat);
	delu_ (va, ua, dua, vn, &nxp, &nyp, &nzp, voxdim, &dt);
	for (j = 0; j < m; j++) {
        	u[j] = ua[j]; v[j] = va[j]; du[j] = dua[j]; b[j] = ba[j];
 	}
 	qdat->dt = dt;
	free (ba); free (dua); free (ua); free (va);
}
   
int main (int argc, char **argv) {
/************/
/* imag I/O */
/************/
	FILE		*fp_img, *fp_out, *fp, *fp_diff, *fp_div;
	IFH		ifh[NIMG], ifhm;
	char		imgroot[NIMG+2][MAXL], mskroot[2][MAXL]; 
        char            imgfile[NIMG+2][MAXL], ifhfile[NIMG+2][MAXL];

/**************/
/* processing */
/**************/
        typ             *img[8], *arrptr[8];
        float           voxdim[NIMG][3], vn[3], fwhmp = FWHMP, fwhmd = FWHMD;
        float           bb[4], factor = 1.7, crit = 1.e-03, *er, margin;
        float           *imag, *img1, *img2, *img3, *b, *v; 
        float           *img4, *img5, *imgd, *imgf, *imgh, *du, *u, *div, dt = T;
        float           minval[2], maxval[2], nr[3], totdiv;
        short           *mask[2], *maskp[2];
        int             imgdim[NIMG][4], mskdim[NIMG][4], orientation, diffdim[4];
	int		isbig, isbigm;
	char		control = '\0';
        int             ii, iter, kj, step = 1, mitr = MITR, nterms = 1, orientm;
        int             m, divdim[4], ip = 1, id = 512;
        long		nx, ny, nz, dimension, size, imgpdim[4];
	long		ix, iy, iz;
	long		nxp, nyp, nzp, jj, marg[3]; 
	QDAT		qdat;
	ITERDAT		*iterdat;
	char		*lables[NVOL] = {"body force x",     "body force y",     "body force z",
					 "velocity x",       "velocity y",       "velocity z",
					 "displacement x",   "displacement y",   "displacement z",
					 "checkvel error x", "checkvel error y", "checkvel error z",
					 "gain field corrected source image",
					 "target vs. doubly corrected source difference",
					 "divergence of displacement field"};

/***************/
/* fluid model */
/***************/
	float		mu = MU, lambda = LAMBDA;

/***********/
/* utility */
/***********/
	char		*ptr, command[MAXL];
	int		i, j, k, c;
	double		p, q;
 
/*********/
/* flags */
/*********/
	int		status = 0, thrsh_flg[2] = {0, 0};
        short           reg = 1, movie_interval = DELINT;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/*	f_init ();	open FORTRAN I/O */

/************************/
/* process command line */
/************************/
        for (k = 0, i = 1; i < argc; i++) {
                if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
                                case 'r': reg = 0;			break;
                                case 'd': dt = atof (ptr); 		*ptr = '\0'; break;
                                case 'e': fwhmd = atof (ptr);		*ptr = '\0'; break; 
                                case 'f': fwhmp = atof (ptr);		*ptr = '\0'; break; 
                                case 'g': movie_interval = atoi (ptr);	*ptr = '\0'; break;
				case 'j': ip = atoi (ptr);		*ptr = '\0'; break;
				case 'k': id = atoi (ptr);		*ptr = '\0'; break;
				case 'l': lambda = atof (ptr);		*ptr = '\0'; break;
				case 'm': mu = atof (ptr);		*ptr = '\0'; break;
                                case 'o': mitr = atoi (ptr);		*ptr = '\0'; break;
                                case 's': getrange (ptr, minval + 0, maxval + 0);	/* img0 counted values */
                               		  thrsh_flg[0] = 1;		*ptr = '\0'; break;
                                case 't': getrange (ptr, minval + 1, maxval + 1);	/* img1 counted values */
                                	  thrsh_flg[1] = 1;		*ptr = '\0'; break;
			}
                } else switch (k) {
			case 0: getroot (argv[i], imgroot[0]); k++; break;
                	case 1: getroot (argv[i], mskroot[0]); k++; break;
	      		case 2: getroot (argv[i], imgroot[1]); k++; break;
               		case 3: getroot (argv[i], mskroot[1]); k++; break;
	        	case 4: getroot (argv[i], imgroot[2]); k++; break;
                }
	}
	if (k < 5) {
		printf ("Usage:\t%s target_imag target_mask source_imag source_mask output_imag\n", program);
                printf ("or:   \t%s target_imag none        source_imag none        output_imag\n", program);
		printf ("e.g., \t%s test_4_0    none        test_4_1    none        test_4_out\n",  program);
		printf ("\toptions\n");
		printf ("\t-d<flt>\tspecify time step (default = %6.4f)\n", T);
		printf ("\t-e<flt>\tspecify displacement field FWHM in mm (default = %6.4f)\n", FWHMD);
		printf ("\t-f<flt>\tspecify gain field FWHM in mm (default = %6.4f)\n", FWHMP);
		printf ("\t-g<int>\twrite div and diff volumes every so many iterations (default = %d)\n", movie_interval); 
		printf ("\t-j<int>\tgain field correct every so many iterations (default = %d)\n",ip);
		printf ("\t-k<int>\tsmooth displacement field every so many iterations (default = %d)\n",id);
		printf ("\t-l<flt>\tspecify lambda (default = %.4f)\n", LAMBDA);
		printf ("\t-m<flt>\tspecify mu (default = %.4f)\n", MU);
		printf ("\t-o<int>\tspecify number of iterations (default = %d)\n", MITR);
		printf ("\t-r\tdisable spatial morph\n");
                printf ("\t-s<flt>[to<flt>] specify range of values for target_imag\n"); 
                printf ("\t-t<flt>[to<flt>] specify range of values for source_imag\n"); 
		exit (1);
	}

/***********************/
/* get 4dfp dimensions */
/***********************/
	for (k = 0; k < NIMG; k++) {
		sprintf (imgfile[k], "%s.4dfp.img", imgroot[k]);
		sprintf (ifhfile[k], "%s.4dfp.ifh", imgroot[k]);
		if (k < 2) {
			if (get_4dfp_dimoe (imgfile[k], imgdim[k], voxdim[k], &orientation, &isbig) < 0
		        || Getifh (ifhfile[k], ifh + k)) errr (program, imgroot[k]);
			nx = imgdim[0][0];
			ny = imgdim[0][1];
			nz = imgdim[0][2];
			dimension = nx * ny * nz;
                	Getifh (imgfile[k], ifh + k);
			orientation = ifh[k].orientation;
			printf ("Reading: %s\n", imgfile[k]);
			printf ("dimensions:%9d%10d%10d\n", nx, ny, nz);
			printf ("mmppix:   %10.4f%10.4f%10.4f\n", ifh[k].mmppix[0], ifh[k].mmppix[1], ifh[k].mmppix[2]);
			printf ("center:   %10.4f%10.4f%10.4f\n", ifh[k].center[0], ifh[k].center[1], ifh[k].center[2]);
			if (!k) imag = (float *) calloc (dimension, sizeof (float));
			mask[k]	= (short *) calloc (dimension, sizeof (short));
			if (!imag || !mask[k]) errm (program);
                	if (strcmp (mskroot[k], "none")) {
                   		sprintf (imgfile[k], "%s.4dfp.img", mskroot[k]);
		   		get_4dfp_dimoe (imgfile[k], mskdim[k], voxdim[k], &orientm, &isbigm);
                   		Getifh (imgfile[k], &ifhm);
                   		kj = orientation - ifhm.orientation;
                   		for (i = 0; i < 3; i++) kj |= (mskdim[k][i] != imgdim[k][i]);
                   		for (i = 0; i < 3; i++) kj |= (fabs(ifhm.mmppix[i] - ifh[k].mmppix[i]) > 0.0001);
                   		if (kj) {
                        		fprintf (stderr, "%s: %s %s mismatch\n", program, imgroot[k], mskroot[k]);
                      			exit (-1);
                   		}
                   		printf ("Reading: %s", imgfile[k]);
                   		read_file_float (imgfile[k], imag, dimension, isbigm);
		   		for (i = 0; i < dimension; i++) mask[k][i] = (short) imag[i];
                	} else {
                   		for (i = 0; i < dimension; i++) mask[k][i] = 1; 
                	}
                	sprintf (imgfile[k], "%s.4dfp.img", imgroot[k]); 
                	read_file_float (imgfile[k], imag, dimension, isbig);
                	if (thrsh_flg[k]) {
                      		for (i = 0; i < dimension; i++) {            
                        		if (imag[i] < minval[k] || imag[i] > maxval[k]) mask[k][i] = 0;
                      		}
                   	}
                   	for (kj = i = 0; i < dimension; i++) if (mask[k][i]) kj++;
                   	if ((float) kj / (float) dimension < 0.1) {
                     		printf ("%s fraction of within-mask voxels (%.4f) must be at least 0.1\n", imgfile[k]);
                     		exit (-1); 
                   	}
               	        printf (" (%d voxels)\n", kj);
                }                      /*    if( k < 2 )        */
		switch (k) {
		case 0: break;	
		case 1:	status = ifh[k].orientation - ifh[0].orientation;
			status != strcmp (ifh[k].imagedata_byte_order, ifh[0].imagedata_byte_order);
			for (i = 0; i < 3; i++) {
				status |= imgdim[k][i] - imgdim[0][i];
				status |= (fabs (voxdim[k][i] - voxdim[0][i]) > 0.0001);
			}
			if (status) {
				fprintf (stderr, "%s: %s %s dimension/endian mismatch\n", program, imgroot[k], imgroot[0]);
				exit (-1);
			}
			break;
		case 2:	for (i = 0; i < 3; i++) {
				imgdim[k][i] = imgdim[0][i];
				voxdim[k][i] = voxdim[0][i];
			}
			if (!control) control = (isbig) ? 'b' : 'l';
			break;
                default: break;
		}
	}

/********************/
/* allocate buffers */
/********************/
/****************************************************************************/
/* Explanation: The constant .7223517245.                                   */
/* To implement the Parzen window the following function is used:           */
/* f(u) = 1 - 6 * u^2 + 6 * u^3,   where:  |u| <= 0.5                       */
/* The value u0 = 0.3611758622 is the value for which f(u) assumes the      */
/* value 0.5. That is, f(u0) = 0.5.                                         */
/* We chose the value 2*u0 as a reasonable value for the margin calculation */
/* 2*u0 = .7223517245.                                                      */
/****************************************************************************/     
	margin = 2. * fwhmp / .7223517245 + 1.; 
	marg[0] = margin / voxdim[0][0];   nxp = npad_ (&nx, marg + 0);
	marg[1] = margin / voxdim[0][1];   nyp = npad_ (&ny, marg + 1);
	marg[2] = margin / voxdim[0][2];   nzp = npad_ (&nz, marg + 2);
	nxp = nx + 16; nyp = ny + 16; nzp = nz + 16;  /* fixed padding */
	printf ("image dimensions %d %d %d padded to %d %d %d \n", nx, ny, nz, nxp, nyp, nzp);
	if (!(img1     =    (typ *) calloc (nxp*nyp*nzp*1, sizeof (typ)))
	||  !(img2     =    (typ *) calloc (nxp*nyp*nzp*1, sizeof (typ)))
	||  !(img3     =    (typ *) calloc (nxp*nyp*nzp*1, sizeof (typ)))
	||  !(img4     =    (typ *) calloc (nxp*nyp*nzp*1, sizeof (typ)))
	||  !(img5     =    (typ *) calloc (nxp*nyp*nzp*1, sizeof (typ)))
	||  !(imgd     =    (typ *) calloc (nxp*nyp*nzp*1, sizeof (typ)))
	||  !(imgf     =    (typ *) calloc (nxp*nyp*nzp*1, sizeof (typ)))
	||  !(imgh     =    (typ *) calloc (nxp*nyp*nzp*1, sizeof (typ)))
	||  !(div      =    (typ *) calloc (nxp*nyp*nzp*1, sizeof (typ)))
	||  !(b        =    (typ *) calloc (nxp*nyp*nzp*3, sizeof (typ)))
	||  !(v        =    (typ *) calloc (nxp*nyp*nzp*3, sizeof (typ)))
	||  !(er       =    (typ *) calloc (nxp*nyp*nzp*3, sizeof (typ)))
	||  !(du       =    (typ *) calloc (nxp*nyp*nzp*3, sizeof (typ)))
	||  !(u        =    (typ *) calloc (nxp*nyp*nzp*3, sizeof (typ)))
	||  !(maskp[0] =    (short *) calloc (nxp*nyp*nzp*1, sizeof (short)))
	||  !(maskp[1] =    (short *) calloc (nxp*nyp*nzp*1, sizeof (short)))
	||  !(iterdat  =    (ITERDAT *) calloc (mitr, sizeof (ITERDAT)))) errm (program);
	mskpac_ (mask[0], &nx, &ny, &nz, maskp[0], &nxp, &nyp, &nzp); 
	mskpac_ (mask[1], &nx, &ny, &nz, maskp[1], &nxp, &nyp, &nzp); 
	size = nxp * nyp * nzp;
	inorms (maskp[0], size, "from the main pgm, Norms of maskp0: ");
	inorms (maskp[1], size, "from the main pgm, Norms of maskp1: ");

/***************/
/* read images */
/***************/
	printf ("Reading: %s\n", imgfile[0]);
	if (!(fp_img = fopen (imgfile[0], "rb"))
	|| gread ((char *) imag, sizeof (float), dimension, fp_img, isbig)
	|| fclose (fp_img)) errr (program, imgfile[0]);
	imgpad_	(imag, &nx, &ny, &nz, img1, &nxp, &nyp, &nzp);
	printf ("Reading: %s\n", imgfile[1]);
	if (!(fp_img = fopen (imgfile[1], "rb"))
	|| gread ((char *) imag, sizeof (float), dimension, fp_img, isbig)
	|| fclose (fp_img)) errr (program, imgfile[1]);
	imgpad_	(imag, &nx, &ny, &nz, img2, &nxp, &nyp, &nzp);
        norms (img1, nr, maskp[0], size, "main pgm(img1): ");
        norms (img2, nr, maskp[0], size, "main pgm(img2): ");
        norms (img3, nr, maskp[0], size, "main pgm(img3): ");
        norms (img4, nr, maskp[0], size, "main pgm(img4): ");
        norms (imgf, nr, maskp[0], size, "main pgm(imgf): ");
        norms (imgh, nr, maskp[0], size, "main pgm(imgh): ");
        norms (imgd, nr, maskp[0], size, "main pgm(imgd): ");

/*******************************/
/* prepare diff and div movies */
/*******************************/
        sprintf (imgfile[NIMG], "%s_diff.4dfp.img", imgroot[NIMG-1]);
        sprintf (ifhfile[NIMG], "%s_diff.4dfp.ifh", imgroot[NIMG-1]);
	if (!(fp_diff = fopen (imgfile[NIMG], "wb"))) errw (program, imgfile[NIMG]);
        sprintf (imgfile[NIMG+1], "%s_div.4dfp.img", imgroot[NIMG-1]);
        sprintf (ifhfile[NIMG+1], "%s_div.4dfp.ifh", imgroot[NIMG-1]);
	if (!(fp_div = fopen (imgfile[NIMG+1], "wb"))) errw (program, imgfile[NIMG+1]);
        for (j = 0; j < 3; j++) diffdim[j] = divdim[j] = imgdim[0][j];
	diffdim[3] = divdim[3] = mitr / movie_interval;    

/***********/
/* compute */
/***********/
	print_date_log ("\nBefore piic...  \n"); 
	piic (img1, img2, img5, imgf, imgh, imgd, bb, maskp[0], maskp[1], 
              voxdim[0], nxp, nyp, nzp, step, &factor, &crit, fwhmp);
	print_date_log ("\nAfter piic... \n");
	copyimg_ (img3,img2,&nxp,&nyp,&nzp);
	imagerr (img1,img3,maskp[0],maskp[1],nxp,nyp,nzp,"img1 & img3 ");
	imagerr (img1,img2,maskp[0],maskp[1],nxp,nyp,nzp,"img1 & img2 ");
	imagerr (img2,img3,maskp[0],maskp[1],nxp,nyp,nzp,"img2 & img3 ");
	imgpdim[0] = nxp; imgpdim[1] = nyp; imgpdim[2] = nzp;
	img[0] = img1; img[1] = img2; img[2] = img3; img[3] = img4; img[4] = img5;
	img[5] = imgd; img[6] = imgf; img[7] = imgh;
	arrptr[0] = u; arrptr[1] = bb; arrptr[2] = b; 
	arrptr[5] = v; arrptr[6] = er; arrptr[7] = du;
	qdat.factor = factor; qdat.crit  = crit;  qdat.mu = mu;
	qdat.lambda = lambda; qdat.fwhmd = fwhmd; qdat.dt = dt;
	qdat.fwhmp  = fwhmp;  qdat.id    = id;    qdat.ip = ip;
	
	m = nxp * nyp * nzp * 3; printf ("\nNr. of unknowns= %ld\n", m);
	for (ii = 0; ii < mitr; ii++) {
		qdat.ii = ii;
		qstep2 (img, arrptr, maskp, voxdim[0], imgpdim, &qdat, iterdat);
		if ( !((ii+1) % movie_interval) ) {
  			for (j = 0; j < size; j++ ) imgd[j] = img1[j] - img3[j];
			outimg (imag, nx, ny, nz, imgd, nxp, nyp, nzp, fp_diff, dimension, control);
       			norms (img3, nr, maskp[0], size, "Norms of the diff image(img):");
		} 
        	dvrgnc_ (u,div,voxdim[0],&nxp,&nyp,&nzp);
		if ( !( (ii + 1) % movie_interval ) ) /*write div movie.*/ 
			outimg (imag, nx, ny, nz, div, nxp, nyp, nzp, fp_div, dimension, control);
  	}
	print_iterdat (program, imgroot[NIMG-1], iterdat, mitr);
        for (totdiv = jj = 0; jj < nxp*nyp*nzp; jj++) totdiv += div[jj];
        printf ("Integral(Div) = %10.4f\n", totdiv);
       	norms (div, nr, maskp[0], size, "the divergence: ");

/****************/
/* write images */
/****************/
	printf ("Writing: %s\n", imgfile[NIMG-1]);
	if (!(fp_out = fopen (imgfile[NIMG-1], "wb"))) errw (program, imgfile[NIMG-1]);
	for (k = 0; k < 3; k++) {
 		outimg (imag,nx,ny,nz,b+k*size,nxp,nyp,nzp,fp_out,dimension, control);
	}
	for (k = 0; k < 3; k++) {
        	outimg (imag,nx,ny,nz,v+k*size,nxp,nyp,nzp,fp_out,dimension, control);
	}
   	for (k = 0; k < 3; k++) {
       		outimg (imag,nx,ny,nz,u+k*size,nxp,nyp,nzp,fp_out,dimension, control);
    	}
	for (k = 0; k < 3; k++) {
		outimg (imag,nx,ny,nz,er+k*size,nxp,nyp,nzp,fp_out,dimension, control);
	}
	outimg (imag,nx,ny,nz,img4,nxp,nyp,nzp,fp_out,dimension, control);
	outimg (imag,nx,ny,nz,imgd,nxp,nyp,nzp,fp_out,dimension, control);
	outimg (imag,nx,ny,nz,div, nxp,nyp,nzp,fp_out,dimension, control);
	if (fclose (fp_out)) errw (program, imgfile[NIMG-1]);

/*******/
/* ifh */
/*******/
        imgdim[NIMG-1][3] = NVOL;
	writeifhe (program, imgroot[NIMG-1], imgdim[NIMG-1], voxdim[NIMG-1], ifh[0].orientation, control);
	writeifhe (program, imgfile[NIMG],   diffdim,        voxdim[NIMG-1], ifh[0].orientation, control);
	writeifhe (program, imgfile[NIMG+1], divdim,         voxdim[NIMG-1], ifh[0].orientation, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr             %s", imgroot[NIMG-1]); status = system (command);
	sprintf (command, "ifh2hdr -r-100to100 %s", imgfile[NIMG]);   status = system (command);
	sprintf (command, "ifh2hdr -r-1to1     %s", imgfile[NIMG+1]); status = system (command);

/*******/
/* rec */
/*******/
	startrecle (imgfile[NIMG-1], argc, argv, rcsid, control);
	catrec     (imgfile[0]);
	catrec     (imgfile[1]);
	for (k = 0; k < NVOL; k++) {
		sprintf (command, "Volume %2d %s\n", k + 1, lables[k]); printrec (command);
	}
        endrec ();

/***************/
/* free memory */
/***************/
	free (imag); free (mask[0]); free (mask[1]);
	free (img1); free (img2); free (img3); free (img4); free (img5); free (imgd); free (imgf); free (imgh);
	free (div); free (b); free (v); free (er); free (du); free (u);
        free (maskp[0]); free (maskp[1]); free (iterdat);
/*	f_exit ();	close FORTRAN I/O */
	exit (status);
}
