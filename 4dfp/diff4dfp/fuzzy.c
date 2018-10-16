#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <knee.h>

/* fuzzy.c - all the fuzzy functions */
/*
** seg_fuzz - perform the fuzzy segmentation
** inputs:
** npix - the number of pixels in the image
** src - the source image 
** npart - the number of partitions 
** eps - min error for convergence
** seed - array with npart elements with initial avg cluster values
** ipart - indicate which of the partitions to monitor for convergence
*/

void
seg_fuzz(short int *src,int npart,int ipart,double eps,double *seed,int npix) 
{
void calc_avg(short int *img,double *g_avg,double *fuzz,int npart,int npix);
void calc_dist(short int *img,double *g_avg,double *g_dst,int npart,int npix);
void calc_fuzz(double *g_dst,double *fuzz,int npart,int npix);
    int i, j;
    double avg_old,g_avg[32], *g_dst, *fuzz;

/* allocate arrays and seed, do not use integers */
    g_dst = (double *) malloc(npix*npart*sizeof(double));
    fuzz = (double *) malloc(npix*npart*sizeof(double));
    if (g_dst==NULL || fuzz==NULL) {
	printf("Error - allocation memory for seg_fuzz\n");
	return;
    }
    for (i=0; i<npart; ++i){
	g_avg[i] = seed[i];
    }
    i=0;
    avg_old = g_avg[ipart];
/* simple loop for no. of iterations or until convergance */
    do {
	++i;
	avg_old = g_avg[ipart];
        calc_dist(src,g_avg,g_dst,npart,npix);
        calc_fuzz(g_dst,fuzz,npart,npix);  
        calc_avg(src,g_avg,fuzz,npart,npix); 
/*
        printf("iteration %d avg %d = %6.3f\n",i,ipart,g_avg[ipart]);
        for (j=0; j<npart; ++j) printf("%d->%6.3f ", j,g_avg[j]);
	printf("\n");
*/
    }  while (fabs(g_avg[ipart]-avg_old) > eps);

    for (j=0; j<npart; ++j) {
	seed[j]=g_avg[j];
        printf("%d->%6.3f\n", j,g_avg[j]);
    }

/* different uses for the fuzz array 
    if (flag == FF_DISP) defuzz(npix, ipart, fuzz, dst);
    else if (flag == FF_AVG) avgfuzz(npix,ipart,oimg,fuzz,dst,
        SPAN,DY-SPAN,SPAN, SPAN,DY-SPAN,SPAN, 0,1,0);
*/

    free(g_avg); free(g_dst); free(fuzz);
}


/*
** calc_avg - calculate the cluster (partition) centers 
** note: integers may cause division by 0 error
*/
void
calc_avg(short int *img,double *g_avg,double *fuzz,int npart,int npix)
{
    int i,j;
    double gnum, den, tmp;

    for (i=0; i<npart; ++i) {
        gnum = den = 0.0;
        for (j=0; j<npix; ++j) {
/*
            tmp = powrr(*(fuzz+i*npix+j),M);
*/
            tmp = (*(fuzz+j*npart+i))*(*(fuzz+j*npart+i));
            gnum += tmp * (*(img+j));
            den += tmp;
        }
if (den==0.0) {
	printf("div by 0 in avg %d %d \n",i,j);
	return;
}
        *(g_avg+i) = gnum / den;
    }
}

/*
** calc_dist - calculate the distance from cluster average to each pixel 
** alternate method: use the average over 3x3 neighbors
*/
void
calc_dist(short int *img,double *g_avg,double *g_dst,int npart,int npix)
{
    int i,j;
    double gsum, tmp;

    for (j=0; j<npix; ++j) {
        gsum = 0.0;
        for (i=0; i<npart; ++i) {
            *(g_dst+j*npart+i) = tmp = ((double)*(img+j)-(*(g_avg+i)))*
		((double)*(img+j)-(*(g_avg+i)));
	    gsum += tmp;
        }
if (gsum==0.0) {
	printf("div by 0 in dist %d \n",j);
        return;
}
        for (i=0; i<npart; ++i) {
            *(g_dst+j*npart+i) /= gsum;
        }
    }
}

/*
** calc_fuzz - fill the fuzzy array for next iteration
*/
void
calc_fuzz(double *g_dst,double *fuzz,int npart,int npix)
{
    int i,j,k;
    double gsum,tmp;

    for (k=0; k<npix; ++k)  {
        for (i=0; i<npart; ++i) {
            gsum = 0.0;
	    tmp = *(g_dst+k*npart+i);
            for (j=0; j<npart; ++j) { 
/*
                gsum += powrr(g_dst[i*npix+k]/g_dst[j*npix+k],M-1);
*/
                gsum += tmp/(*(g_dst+k*npart+j));
            }
            *(fuzz+k*npart+i) = 1.0/gsum;
        }
    }
}


/*
** name  - defuzz
** use   - convert a fuzzy partition to an image
** inputs:
** npixel - number of pixels
** ipart - the partition to use in the defuzz
** fuzz - access to the fuzzy partition array
** outputs:
** dst - an image created based on the fuzzy partition
*/
void 
defuzz(int npixel, int ipart, double *fuzz, unsigned char *dst)
{
    int i,index;

    index = ipart*npixel;
    for (i=0; i<npixel; ++i, ++dst) {
        *dst = *(fuzz+index+i)*NGRAY;
    }
}


/*
** avgfuzz - get fuzzy average of pixels in image under a mask
*/
void
avgfuzz(int npixel,int ipart, unsigned char *src,double *fuzz, unsigned char *dst,
    int xmin, int xmax, int xspan, int ymin, int ymax,int yspan,int zmin,int zmax,int zspan)
{
    int i,j,k, x,y,z;
    int xsmin, xsmax, ysmin, ysmax, zsmin, zsmax;
    double sum, num;

    fuzz += ipart*npixel;

    for (z=zmin; z<zmax; ++z) 
      for (y=ymin; y<ymax; ++y) 
        for (x=xmin; x<xmax; ++x) {
          zsmin=z-zspan; zsmax=z+zspan;
          ysmin=y-yspan; ysmax=y+yspan;
          xsmin=x-xspan; xsmax=x+xspan;
          sum = 0.0; num = 0.0;
          for (k=zsmin; k<=zsmax; ++k) 
            for (j=ysmin; j<=ysmax; ++j) 
              for (i=xsmin; i<=xsmax; ++i) {
                sum += *(src+k*DZ+j*DY+i)* *(fuzz+k*DZ+j*DY+i);
                num += *(fuzz+k*DZ+j*DY+i);
          }
          if (num==0.0) *(dst+z*DZ+y*DY+x) = 0x00;
          else *(dst+z*DZ+y*DY+x) = (unsigned char) (sum/num);
    }
}

/*
** name - hough_line
** use  - find lines in the image
** 1024*2048 = 2097152
** 1024*4096 = 4194304
*/
void 
hough_line(short int *src, int ncol, int nrow)
{
    int i,j,ix,iy,irho,iang;
    int size = 2097152;
    short int *tmp, *hough;
    double cs[1024],sn[1024],deg[1024];
    double r1,t1,r2,t2,r3,t3;

/*
    hough = (short int *)malloc(size*sizeof(short int));
    if (hough == (short int *)NULL) {
	printf("Error - allocating hough\n");
	return;
    }
*/
/* init the hough, the cos and sin arrays */
/*
    tmp=(short int *)hough; 
    for (i=0; i<size; ++i, ++tmp) *tmp=0;
    for (i=0; i<1024; ++i) { 
	cs[i] = cos((double)i*M_PI/512.0 - M_PI);
	sn[i] = sin((double)i*M_PI/512.0 - M_PI);
	deg[i] = ((double)i*M_PI/512.0 - M_PI)*180/M_PI;
    }
*/
/* calc lines femur */
    t1 = ang_norm(244,326,408,439);
    r1 = rho_norm(244,326,&t1);

    t2 = ang_slope(244,326,408,439)+15.0*M_PI/180.0;
    t3 = ang_slope(244,326,408,439)-30.0*M_PI/180.0;
    r2 = rho_norm(244,326,&t2);
    r3 = rho_norm(408,439,&t3);

/* tibia 
    t1 = ang_norm(486,651,619,443);
    r1 = rho_norm(486,651,&t1);

    t2 = ang_slope(486,651,619,443)-15.0*M_PI/180.0;
    t3 = ang_slope(486,651,619,443)+30.0*M_PI/180.0;
    r2 = rho_norm(486,651,&t2);
    r3 = rho_norm(619,443,&t3);
*/
    printf("ang 1 %f %f 2 %f %f 3 %f %f\n",r1,RTOD(t1),
	r2,RTOD(t2),r3,RTOD(t3));

/* create a mask 
    sector_flag(src,r1,t1,r2,t2,r3,t3,6,ncol,nrow);
    draw_line(src,244,326,408,439,ncol);
*/
return;


/* cont into the hough array 
    for (ix=0; ix<ncol; ++ix) {
    for (iy=0; iy<nrow; ++iy) {
*/
/* for femur 
	if (!sector_flag(ix,iy,r1,t1,1,r2,t2,1,r3,t3,0)) continue;
*/
/* for the tibia 
	if (!sector_flag(ix,iy,r1,t1,0,r2,t2,0,r3,t3,0)) continue;

	if (*(src+iy*ncol+ix)==0) continue;
	for (j=0; j<1024; ++j) {
*/
/* cut on angle femur 
	    if (j>607 || j<507) continue;
*/
/* tibia 
	    if (j>903 || j<803) continue;

	    irho = (int)((double)ix*cs[j] + (double)iy*sn[j]);
	    if (irho>=0 && irho<2048)
	    (*(hough+2048*j+irho))++;
	}
    }}
*/
/*
put_image((unsigned char *)hough, "hough.img",2*size);
*/
/*
    tmp=(short int *)hough;  
    for (i=0; i<size; ++i, ++tmp) {
femur line cut 
	if (*tmp>=300) { 
*/
/* tibia 
	if (*tmp>=300) { 
	    irho=i%2048; iang=i/2048;
	    printf("%d rho %d ang %d %f\n",
		*tmp,irho,iang,deg[iang]);
/* vary y and find x 
	    if (iang<128 || iang>896 || (iang>384 && iang<640)) { 
		for (iy=0; iy<nrow; ++iy) {
		    ix=(int)(((double)irho - (double)iy*sn[iang])/cs[iang]);
		    if (ix>=0 && ix<ncol)
			*(src+iy*ncol+ix) = NGRAY-1;
		}
	    }
/* vary x and find y 
	    else {  
		for (ix=0; ix<ncol; ++ix) {
		    iy=(int)(((double)irho - (double)ix*cs[iang])/sn[iang]);
		    if (iy>=0 && iy<nrow)
			*(src+iy*ncol+ix) = NGRAY-1;
		}
	    }
	}
    }
*/
}


/* 
** name - ang_norm, rho_norm
** use - returns the angle and rho of normal representation of line 
** 	the angle can be from -pi to +pi, and rho>0
*/
double 
ang_norm(int x1, int y1, int x2, int y2)
{
    double dx,dy,ang;

    if (x1==x2 && y1==y2) return(0.0);
    dx = (double)(x2-x1);
    dy = (double)(y2-y1);
    ang = atan2(-1.0*dx, dy);
    return(ang);
}

/* 
** name - ang_slope
** use - returns the slope angle of line determined by 
** 	two points between -pi/2 and +pi/2
*/
double 
ang_slope(int x1, int y1, int x2, int y2)
{
    double dx,dy,ang;

    if (x1==x2 && y1==y2) return(0.0);
    dx = (double)(x2-x1);
    dy = (double)(y2-y1);
    ang = atan2(dy, dx);
    if (ang>M_PI/2.0) ang -= M_PI;
    else if (ang<M_PI/-2.0) ang += M_PI;
    return(ang);
}

double 
rho_norm(int x1, int y1, double *ang)
{
    double rho;

    rho = (double)x1*cos(*ang)+(double)y1*sin(*ang);
    if (rho<0.0) { 
	if (*ang<=0.0) *ang += M_PI;
	else *ang -= M_PI;
        rho *= -1.0;
    }
    return(rho);
}


/*
** name-sector_mask
** use- The dst is the src anded with mask. 
**	Mask is defined by a point, direction, and width
** 	The mask is created by 3 lines in parametric form
**	The sector is selected by id = [0,7].
*/
void
sector_mask(short int *src, short int *dst,int x1,int y1,
	double dir,double wid,int ncol,int npix)
{
    int i, id;
    short int *tmp;
    double r1,r2,r3,t1,t2,t3;
    int x2,x3,x4,y2,y3,y4;

    tmp = (short int *)malloc(npix*sizeof(short int));
    if (tmp == (short int *)NULL) {
	printf("Error - allocating mem in sector_mask\n");
	return;
    }
    for (i=0; i<npix; ++i) *(tmp+i)=0;

/* calculate the params for the major line */
    t1 = norm_ang(dir);
    r1 = rho_norm(x1,y1,&t1);
/* calc points for the other lines */
    x2 =(int)((double)x1 - sin(dir)*wid);
    x3 =(int)((double)x1 + sin(dir)*wid);
    y2 =(int)((double)y1 + cos(dir)*wid);
    y3 =(int)((double)y1 - cos(dir)*wid);
    t2 = norm_ang(dir+0.5*M_PI);
    t3 = norm_ang(dir+0.5*M_PI);
    r2 = rho_norm(x2,y2,&t2);
    r3 = rho_norm(x3,y3,&t3);

/* fill the mask array */
    sector_fill(tmp,r1,t1,4,ncol,npix);
    sector_fill(tmp,r2,t2,2,ncol,npix);
    sector_fill(tmp,r3,t3,1,ncol,npix);

/* find the id of the sector of choice */
    x4 = x1 + 10.0*cos(dir);
    y4 = y1 + 10.0*sin(dir);
    id = *(tmp+ (int)y4*ncol + (int)x4);

    for (i=0; i<npix; ++i) {
	if (*(tmp+i)==id) *(tmp+i)=NGRAY-1;
	else *(tmp+i)=0;
    }

put_image((unsigned char *)tmp, "temp.img",2*npix);

    and_img(src,tmp,dst,npix);
}

/*
** name-sector_fill
** use- Add id to each element of dst that is further from the origin
**	than the line defined by rho and ang. The size of dst is defined
**	by ncol and npix.
*/
void 
sector_fill(short int *dst,double rho,double ang,int id,int ncol,int npix)
{
    int ix,iy,nrow;
    double sn,cs,x,y;

    nrow = npix/ncol;
    sn = sin(ang);
    cs = cos(ang);

    if (ang>=0.25*M_PI && ang<0.75*M_PI) {
	for (ix=0; ix<ncol; ++ix) {
            x = (double)ix;
	    y = (rho - x*cs)/sn;
	    if ((int)y>=nrow) continue;
            else if (y<0.0) y=0.0;
	    for (iy=(int)y; iy<nrow; ++iy) 
		*(dst+iy*ncol+ix) += id;
	}
    }
    else if (ang<-0.25*M_PI && ang>= -0.75*M_PI) {
	for (ix=0; ix<ncol; ++ix) {
            x = (double)ix;
	    y = (rho - x*cs)/sn;
	    if ((int)y>nrow) y=nrow;
	    else if (y<0.0) continue;
	    for (iy=0; iy<(int)y; ++iy) 
		*(dst+iy*ncol+ix) += id;
	}
    }
    else if (ang>= -0.25*M_PI && ang<0.25*M_PI) {
	for (iy=0; iy<nrow; ++iy) {
            y = (double)iy;
	    x = (rho - y*sn)/cs;
	    if ((int)x>=ncol) continue;
	    else if (x<0.0) x=0.0;
	    for (ix=(int)x; ix<ncol; ++ix) 
		*(dst+iy*ncol+ix) += id;
	}
    }
    else {
	for (iy=0; iy<nrow; ++iy) {
            y = (double)iy;
	    x = (rho - y*sn)/cs;
	    if ((int)x>ncol) x=ncol;
	    else if (x<0.0) continue;
	    for (ix=0; ix<(int)x; ++ix) 
		*(dst+iy*ncol+ix) += id;
	}
    }
}


/*
** name - hough_circ
** use  - find circles in the image
** 1024*2048 = 2097152
** 1024*4096 = 4194304
*/
void 
hough_circ(short int *mask,int *xp,int *yp,int npts,int ncol,int npix)
{
    int i,j,ix,iy,irho,iang,ind;
    int size = 2097152;
    short int *tmp, *hough;
    double r1,t1,r2,t2,r3,t3;

    tmp = (short int *)malloc(npix*sizeof(short int));
    hough = (short int *)malloc(size*sizeof(short int));
    if (hough == (short int *)NULL) {
	printf("Error - allocating hough\n");
	return;
    }

/* mask out part of the mask */
    for (i=0; i<npix; ++i) {
	if (i%ncol>320) *(mask+i)=0;
    }
    for (j=0; j<npts; ++j) {
	*(tmp+yp[j]*ncol+xp[j])=NGRAY-1;
    }

/* clear the hough array */
    for (i=0; i<size; ++i) *(hough+i)=0;
/* loop on all pixels */
    ind=0;
    for (i=0; i<npix; ++i) {
/* if not in mask skip */
	if (*(mask+i)==0) continue;
	++ind;
	if (ind>=4*2048) continue;
	ix=i%ncol; iy=i/ncol;
	for (j=0; j<npts; ++j) {
	    if (xp[j]>320) continue;
	    irho = (ix-xp[j])*(ix-xp[j])+(iy-yp[j])*(iy-yp[j]);
	    r1 = sqrt((double)irho);
	    irho = (int)rint(r1);
	    if (irho<100) 
	        (*(hough+ind*256+irho))++;
	}
    }

/* stats */
    hist_img(hough, size);
/* print out */
    ind=0;
    for (i=0; i<npix; ++i) {
	if (*(mask+i)==0) continue;
	++ind;
	if (ind>=4*2048) continue;
	ix=i%ncol; iy=i/ncol;
        for (j=70; j<128; ++j) {
	    if (*(hough+ind*256+j)>45) {
		printf("circ %d %d rad %d val %d\n",
		ix,iy,j,*(hough+ind*256+j));
		draw_circle(tmp,ix,iy,(double)j,ncol,npix);
		*(tmp+i)=NGRAY-1;
	    }
	}
    }
    put_image((unsigned char *)tmp, "tmp1.img",2*npix);
}

/* 
** name - norm_ang, norm_angp
** normalize the angle value
** norm_ang: to range of +/-PI while preserving direction
** 	assume input between +/-2.0*PI
** norm_angp: to range of 0-PI losing direction 
** 	assume input between +/-PI
*/
double 
norm_ang(double ang)
{
    if (ang>M_PI) return(ang - 2.0*M_PI);
    if (ang<-M_PI) return(ang + 2.0*M_PI);
    else return(ang);
}
double 
norm_angp(double ang)
{
    if (ang<0.0) return(ang+M_PI);
    if (ang>=M_PI) return(ang-M_PI);
    else return(ang);
}
