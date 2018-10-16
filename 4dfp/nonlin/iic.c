/* 2003.02.14 - module to do Intensity Inhomogeneity Correction (IIC).*/
/*$Header: /data/petsun4/data1/src_solaris/nonlin/RCS/iic.c,v 2.3 2004/10/28 21:55:51 rsachs Exp $*/
/*$Log: iic.c,v $
 * Revision 2.3  2004/10/28  21:55:51  rsachs
 * Undefined #Debug.
 *
 * Revision 2.2  2004/09/28  23:31:15  rsachs
 * Reinstalled computing the diff-image. Aesthetics.
 *
 * Revision 2.1  2003/11/26  23:21:39  rsachs
 * Modified computation of the diff image in function imgerr.
 *
 * Revision 2.0  2003/04/22  20:48:34  rsachs
 * New final version. Abandoning the attempt to take advantage of the gcomp matrix symmetry.
 *
 * Revision 1.8  2003/04/21  23:45:53  rsachs
 * Final Version for now.
 *
 * Revision 1.7  2003/04/17  19:03:58  rsachs
 * Unsuccessful optimisation effort.
 *
 * Revision 1.6  2003/04/15  15:56:17  rsachs
 * About to perform radical surgery, to improve speed.
 *
 * Revision 1.5  2003/04/10  23:57:35  rsachs
 * Recovered after deletion.
 *
 * Revision 1.4  2003/04/08  19:06:56  rsachs
 * Works well. Needs optimisation.
 *
 * Revision 1.3  2003/03/18  19:27:41  rsachs
 * Final version b4 adding masks..
 *
 * Revision 1.2  2003/03/10  21:58:10  rsachs
 * Module to do Intensity Inhomogeneity Correction (iic) on paired images.
 **/
 
 #include <stdlib.h>
 #include <stdio.h>
 #include <math.h>
 typedef float typ; 
 #define TWOPI 6.283185307179586
/* #define Debug */

 extern void	ludcmp (typ **a, int n, int *mdx, typ *d);
 extern void 	lubksb (typ **a, int n, int *mdx, typ *b); 
 extern void 	hcomp (typ *img0 , typ *img1, typ *imgg,typ *imgh,
  		short *mask0,short *mask1,float *voxdim,typ *bb,int nx,int ny, int nz,
  		int step);
 extern void 	gcomp (typ *img0, typ *img1, typ *imgg, typ *imgh, 
        	typ *imgd, short *mask0, short *mask1, typ *bb, typ *bbb, 
        	float *voxdim, int nx, int ny, int nz, int nterm, int step); 
 extern void 	imgerr (typ *img0, typ *img1, typ *imgg, typ *imgh, 
        	typ *imgd, short *mask0, short *mask1, int nx, int ny, int nz, 
        	typ *err, char *errstr); 
 typ 		gu (typ u, typ Lu, int i); 
 void 		tables (typ **xtab, typ **ytab, typ **ztab, float *voxdim, typ Lx,
             	typ Ly, typ Lz, int nx, int ny, int nz, int nterm1);
 void 		alloctab (typ ***ptab, int m, int n, char *s);
 extern void 	norms (typ *img, typ *nr, short *mask, int size, char *str);
 extern void 	inorms (short *mask, long size, char *str);
/*----------------------------------------------------------------------
 *  void iic() - function to do iic.
 *  Definition of Arguments:
 * img0,img1   - the input images. img1 = the image that need b corrected.
 * img2        - the output image. Note: img0,img1 r not modified.
 * imgg,imgh   - the g & h volumes.
 * imgd        - the difference image. 
 *--------------------------------------------------------------------*/
  void iic (typ *img0, typ *img1, typ *img2, typ *imgg, typ *imgh, typ 
            *imgd, typ *bb, typ *bbb, short *mask0, short *mask1, float *voxdim, 
            int nx, int ny, int nz, int nterm, int step, int iter)
{
	char errstr[80];
	int j, size, ndx, M;
	typ err[3], Lx, Ly, Lz, nr[3];

	size = nx * ny * nz;             /* initialising the h & g volumes */
	for (ndx = 0; ndx < size; ndx++)         
      		imgg[ndx] = imgh[ndx] = 1.0;

	imgerr (img0,img1,imgg,imgh,imgd,mask0,mask1,nx,ny,nz,err,"0 & 1");
	norms (imgd,nr,mask0,size," imgdi "); 

	for (j = 0; j < iter; j++) {
		sprintf (errstr, "\nAfter %d iterations", j);
   		hcomp (img0,img1,imgg,imgh,mask0,mask1,voxdim,bb,nx,ny,nz,step);
   		gcomp (img0,img1,imgg,imgh,imgd,mask0,mask1,bb,bbb,voxdim,nx,ny,nz,nterm,step); 
   		for (ndx = 0; ndx < size; ndx++) {
         		if (!( mask0[ndx] && mask1[ndx])) continue;
         		img2[ndx] = imgg[ndx] * imgh[ndx] * img1[ndx];
   		} 
		norms (img2,nr,mask0,size," img2 "); 
   		imgerr (img0,img1,imgg,imgh,imgd,mask0,mask1,nx,ny,nz,err,"0&1 final");
		norms (imgd,nr,mask0,size," imgdf "); 
 	}
}                                          /*  End iic()   */  
 
/*----------------------------------------------------------------------
// void gcomp - the g-function computer.
// Definition of the variables:
// Input: 
//       img0, img - input images. img is the 1 that need b corrected.
//       bb  - the linear gradient vector. h(x,y,z)= b0+b1*x+b2*y+b3*z.
//       bbb - vector of Fourier series coefficients. The iic output.
//       nx, ny, nz - nr. of image voxels in the respective direction.
//       Lx, Ly, Lz - image's dimensions in the respective directions.
//       nterm - number of terms taken in the Fourier expansion.
//       step - density control parameter.
// Output:
//       imgg - the g-volume. 
//--------------------------------------------------------------------*/

void gcomp (typ *img0, typ *img1, typ *imgg, typ *imgh, typ *imgd,  
            short *mask0, short *mask1, typ *bb, typ *bbb, float *voxdim, 
            int nx, int ny, int nz, int nterm, int step)
{
	static typ 	**a, d, **xtab, **ytab, **ztab, gnorms[3];
	static typ 	hhh, *ztabk, *ytabj, *xtabi, gx, gy, gz, gix, giy, giz;
	static typ 	g, gi, Lx, Ly, Lz, x, y, z, val, valsq, val0, vm;
	static typ 	hvv, hvs, hvsg, giyz;
 	static int 	*mdx, i, j, k, i1, j1, k1, i2, j2 ,k2, nterm1, ii, jj, M;
	static int 	jnt, knt, nterm2;
 	int 		ndx = 0;

 	printf ("\nEntering gcomp...\n");
	system ("date"); 
 	nterm1 = 2 * nterm + 1; nterm2 = nterm1 * nterm1;
 	M  = nterm1 * nterm1 * nterm1;
 	a = (typ **) calloc (M, sizeof(typ *));
 	if( !a ) errm ("iic");  
	for (j = 0; j < M; j++) {
     		a[j] = (typ *) calloc (M, sizeof(typ));
     		if (!a[j] ) errm("iic");
 	}

	mdx = (int *) calloc (M, sizeof(int)); 
	if (!mdx) errm ("iic"); 

	for (i = 0; i < M; i++) {
    		bbb[i] = 0.0;
    		for (j = 0; j < M; j++)
       			a[i][j] = 0.0;
 	}

	alloctab (&xtab, nx, nterm1, "The x-table"); 
	alloctab (&ytab, ny, nterm1, "The y-table"); 
	alloctab (&ztab, nz, nterm1, "The z-table"); 
	Lx = nx * voxdim[0]; Ly = ny * voxdim[1]; Lz = nz * voxdim[2];
	tables (xtab, ytab, ztab, voxdim, Lx, Ly, Lz, nx, ny, nz, nterm1);

	for (k = 0; k < nz; k += step) {
     		z = k * voxdim[2]; ztabk = ztab[k]; 
     		for (j = 0; j < ny; j += step) {
         		y = j * voxdim[1]; ytabj = ytab[j];
         			for (i = 0; i < nx; i += step) {
             				ndx = i + nx * ( j + ny * k );
             				x = i * voxdim[0]; xtabi = xtab[i];
             				if (!( mask0[ndx] && mask1[ndx])) continue; 
             				val  = img1[ndx];   valsq = val  * val;
             				val0 = img0[ndx];   /* hhh = hh(x, y, z); */
             				hhh  = imgh[ndx]; hvv = hhh * val * val0;
             				hvs  = hhh * hhh * valsq;
    
             				for (k1 = 0; k1 < nterm1; k1++) {
                 				gz = ztabk[k1]; 
                 				for (j1 = 0; j1 < nterm1; j1++) {
                    					gy = ytabj[j1]; 
                     					for (i1 = 0; i1 < nterm1; i1++) {
                         					gx = xtabi[i1]; 
                         					g = gx * gy * gz;     hvsg = hvs * g;
                         					ii = nterm1 * (nterm1 * k1 + j1) + i1;
                         					bbb[ii] += hvv * g; /*bbb[ii] += hhh * val * val0 * g;*/ 
                             					for (k2 = 0, knt = 0; k2 < nterm1; k2++, knt += nterm2) {
                             						giz = ztab[k][k2];  
                             						for (j2 = jnt = 0; j2 < nterm1; j2++, jnt+= nterm1) {
                                 						giy = ytab[j][j2]; giyz = giy * giz * hvsg;
                                 						for (i2 = 0; i2 < nterm1; i2++) {
                                     							jj = knt + jnt + i2;
                                     							if (jj >= ii) a[ii][jj] += giyz * xtab[i][i2];
                                     							else a[ii][jj] = a[jj][ii];
                                 						}
                             						}
                         					}
                     					}
                 				}
             				}
         			}
     			}
 		}


/* Solve now the M x M linear equation system. */
	printf ("Tackling the 4x4 linear system...\n");
	system ("date"); 
	ludcmp (a, M, mdx, &d);
	lubksb (a, M, mdx, bbb);

	printf ("\niic solution (the bbb-vector)\n");
	for (i = 0; i < M; i++) {
     		printf ("%10.4f  ", bbb[i]);
 	}
/*   Compute the g-volume.     */
	for (k = 0; k < nz; k++) {
   		z = k * voxdim[2];
   		for (j = 0; j < ny; j++) {
     			y = j * voxdim[1];
     			for (i = 0; i < nx; i++) {
       				x = i * voxdim[0];
       				for (k1 = 0; k1 < nterm1; k1++) {
         				gz = ztab[k][k1];
         				for (j1 = 0; j1 < nterm1; j1++) {
           					gy = xtab[j][j1];
           					for (i1 = 0; i1 < nterm1; i1++) {
             						ndx = i1 + nx * (j1 + ny * k1);
             						gx = xtab[i][i1];
             						ii = nterm1 * (nterm1 * k1 + j1) + i1;
             						imgg[ndx] = bbb[ii] * gx * gy * gz;  
           					}
         				}
       				}
     			}
   		}
 	}
/* freeing pointers. */ 
	for (j = 0; j < M; j++) free(a[j]);
	free (a); free (mdx);
	free (xtab); free (ytab); free (ztab);
	norms (imgg, gnorms, mask0, nx*ny*nz, " imgg (at end of gcomp)"); 
	printf ("Exiting gcomp\n");
	system ("date");
}                                     /*  End gcomp  */

void tables (typ **xtab, typ **ytab, typ **ztab, float *voxdim,
             typ Lx, typ Ly, typ Lz, int nx, int ny, int nz, int nterm1)
{
	typ 	dim, x, y, z, Lu, u, **ptab, *ptabj;
	int 	i, j, k, jj, kj, nu;

 	for (i = 0; i < 3; i++) {
    		if (i == 0)  { nu = nx; Lu = Lx; ptab = xtab; }
    		if (i == 1) { nu = ny; Lu = Ly; ptab = ytab; }
    		if (i == 2) { nu = nz; Lu = Lz; ptab = ztab; }
    		dim = voxdim[i];
    		for (jj = 0, u = 0.0; jj < nu; jj++, u += dim) {
         		ptabj = ptab[jj]; 
         		for (kj = 0; kj < nterm1; kj++) { 
              			if (kj != 0) {
               				j = kj % 2; 
                  			k = (kj + 1) / 2;
                  			if (j == 0) 
                      				ptabj[kj] = sin (TWOPI * u * k / Lu);
                  			else 
                      				ptabj[kj] = cos(TWOPI * u * k / Lu);
               			} else {
                      			ptabj[kj] = 1.0;
               			} 
         		}
    		}
 	}
}                      /*  End tables()  */ 

typ gu (typ u, typ Lu, int i)
{
	int 	j, k;

	if (i != 0) { 
   		j = i % 2;           
   		k = (i + 1) / 2;     
   		if (j == 0) {
       			return (sin(TWOPI * u * (float)k / Lu));
   		} else {
       			return (cos(TWOPI * u * (float)k / Lu));
   		}      /* Endelse */
 	} else {
   		return 1.0;
 	}        
}         /* End gu() */ 

void alloctab (typ ***ptab, int m, int n, char *s)
{
	int 	j;
	typ 	**ptabl;

	ptabl = (typ **) calloc (m, sizeof(typ *));
	if (!ptabl) errm ("alloctab");

	for (j = 0; j < m; j++) {
    		ptabl[j] = (typ *) calloc (n, sizeof(typ));
    		if (!ptabl[j]) errm ("alloctab");
 	}
 	*ptab = ptabl;
}                           /* End alloctab */ 

/*---------------------------------------------------------------------
// void hcomp - function to compute the linear gradient that may be
// present in an image. Returns, in vector b, the coefficients of the 
// form b0 + b1 * x + b2 * y + b3 * z that minimise the error.
// Uses: Values of the g-function (volume imgg).
// Output:
//        imgh - the h-volume.
//-------------------------------------------------------------------*/
void hcomp (typ *img0, typ *img1, typ *imgg, typ *imgh, short *mask0, 
	    short *mask1, float *voxdim, typ *bb, int nx, int ny, int nz, int step) 
{
	static 	typ a[4][4], *pa[4], gv, gvs, d, x, y, z, val, valsq, val0, vm;
	static 	typ hnorms[3];
	int 	mdx[4], i, j, k, ndx = 0;
	static 	int iter = 0, size;

	for (i = 0; i < 4; i++) {
   		bb[i] = 0.0;   
   		for (j = 0; j < 4; j++)
     			a[i][j] = 0.0;
 	}

	size = nx * ny * nz;
	inorms (mask0,(long)size,"\nnorms of mask0(b4 hcomp): ");
	inorms (mask1,(long)size,"\nnorms of mask1(b4 hcomp): ");
  
	norms (img0,hnorms,mask0,size,"img0(b4 hcomp): ");
	norms (img1,hnorms,mask0,size,"img1(b4 hcomp): ");
	norms (imgg,hnorms,mask0,size,"imgf(b4 hcomp): ");
	norms (imgh,hnorms,mask0,size,"imgh(b4 hcomp): ");

	for (k = 0; k < nz; k += step) {
     		z = k * voxdim[2]; 
     		for (j = 0; j < ny; j += step) {
         		y = j * voxdim[1];
         		for (i = 0; i < nx; i += step) {
             			ndx = i + nx * ( j + ny * k );
             			x = i * voxdim[0];
             			gv   =  imgg[ndx];  gvs = gv * gv; 
             			val  =  img1[ndx];  valsq  = val  * val; 
             			val0 =  img0[ndx];  /* ndx++; */
             			if (!( mask0[ndx] && mask1[ndx])) continue; 

             			a[0][0] += valsq;      
             			a[0][1] += gvs * x * valsq; 
             			a[0][2] += gvs * y * valsq; 
             			a[0][3] += gvs * z * valsq;  
             			bb[0]   += gv  * val * val0;
       
             			a[1][0] += gvs * x * valsq;  
             			a[1][1] += gvs * x * x * valsq; 
             			a[1][2] += gvs * x * y * valsq;
             			a[1][3] += gvs * x * z * valsq; 
             			bb[1]   += gv  * x * val * val0; 
       
             			a[2][0] += gvs * y * valsq; 
             			a[2][1] += gvs * x * y * valsq; 
             			a[2][2] += gvs * y * y * valsq;
             			a[2][3] += gvs * y * z * valsq;
             			bb[2]   += gv  * y * val * val0;

             			a[3][0] += gvs * z * valsq; 
             			a[3][1] += gvs * x * z * valsq;
             			a[3][2] += gvs * y * z * valsq;
             			a[3][3] += gvs * z * z * valsq;
             			bb[3]   += gv  * z * val * val0;
     			}
   		}
 	}

/* solve now the 4 x 4 linear system. */
	printf ("\nvoxdim=%10.4f %10.4f %10.4f\n",voxdim[0],voxdim[1],voxdim[2]);

 	for (j = 0; j < 4; j++)  pa[j] = a[j];
 		#ifdef Debug
   			printf ("hcomp: the coefficients matrix\n");
   			for (i = 0; i < 4; i++) {
      				putchar ('\n');
      				for (j = 0; j < 4; j++)
           				printf ("%10.4e   ", a[i][j]);
				      	printf ("         %10.4e", bb[i]);
   			} 
 		#endif
	ludcmp (pa, 4, mdx, &d);
	lubksb (pa, 4, mdx, bb);

	iter++; 
	printf ("\nIter # %d The linear gradient solution (bb-vector)\n",iter);
	for (i = 0; i < 4; i++) {
     		printf ("%10.4f  ", bb[i]);
 	}
/* Compute the  h-volume.  */
 	for (k = 0; k < nz; k++) {
     		z = k * voxdim[2];
     		for (j = 0; j < ny; j++) {
         		y = j * voxdim[1];
         		for (i = 0; i < nx; i++) {
             			x = i * voxdim[0];
             			ndx = i + nx * ( j + ny * k );
             			imgh[ndx] = bb[0] + bb[1] * x + bb[2] * y + bb[3] * z;
         		}
     		}
  	}
 	norms (imgh, hnorms, mask0, nx*ny*nz, " imgh (at end of hcomp) "); 
}            /* End hcomp() */

/*---------------------------------------------------------------------
// void imgerr() - function to measure the error (=difference) between
// 2 images. Both images r supposed to have the same sizes. The error(s)
// in the vector err, r the sub1, sub2 & sub-infinity norms of the diff.
//--------------------------------------------------------------------*/

void imgerr (typ *img0, typ *img1, typ *imgg, typ *imgh, typ *imgd, short *mask0,
             short *mask1, int nx, int ny, int nz, typ *err, char *errstr)
{
	typ 	x, err1 = 0.0, err2 = 0.0, erri = 0.0, nr[3];   
	long 	int ndx, i, j, k, size, voxcnt = 0;

	size = nx * ny * nz;
	norms (img0,nr,mask0,size," img0 ");
	norms (imgg,nr,mask0,size," imgf ");
	norms (imgh,nr,mask0,size," imgh ");
	norms (img1,nr,mask0,size," img1 ");
	for (ndx = 0; ndx < size; ndx++) {
          	if (!( mask0[ndx] && mask1[ndx])) continue; 
        	imgd[ndx] = img0[ndx] - img1[ndx];       
          	voxcnt++;
          	x = fabs (img0[ndx] - img1[ndx]);
          	err1 += x;
          	err2 += x * x;
          	if (x > erri) erri = x;
  	}
	err[0] = err1; err[1] = sqrt (err2) / (float) voxcnt; err[2] = erri;
	printf ("\n %s The differences between the images are: %10.4f %10.4f %10.4f\n", errstr, err[0], err[1], err[2]);

}

void norms (typ *img, typ *nr, short *mask, int size, char *str)
{
	typ 	x;
	int 	ndx, cnt = 0;

	nr[0] = nr[1] = nr[2] = 0.0;
	for (ndx = 0; ndx < size; ndx++) {
     		if (!mask[ndx]) continue;
     		cnt++;
     		x = fabs (img[ndx]);
     		nr[0] += x;
     		nr[1] += x * x;
     		if (x > nr[2]) nr[2] = x;
  	}
	nr[0] = nr[0] / (typ) cnt;
  	nr[1] = sqrt (nr[1]) / (typ) cnt;
  	printf ("\nNorms of %s %10.4e %10.4e %10.4f   ", str, nr[0], nr[1], nr[2]); 
	printf ("Voxcnt = %d\n", cnt);
}

void inorms (short *mask, long size, char *str)
{
	int 	nr[2] = {0, 0}, j;
 
	for (j = 0; j < size; j++) {
     		nr[0] += mask[j];
     		if (mask[j] > nr[1]) nr[1] = mask[j];
 	}
	printf ("%s %d %d\n", str, nr[0], nr[1]);
}

typ imagerr (typ *img0,typ *img1,short *mask0,short *mask1,int nx, int ny,
             int nz, char *str)
{
	typ 	err[] = {0.0, 0.0, 0.0}, aux;
	int 	size, j;

	size = nx * ny * nz;
	for(j = 0; j < size; j++) {
/*     if( !( mask0[j] && mask1[j] ) ) continue;   */
     		aux = fabs(img0[j] - img1[j]);
     		err[0] += aux;
     		err[1] += (aux * aux);
     		if( aux > err[2] ) err[2] = aux;
 	}
	err[0] = err[0]/(typ)size;
 	err[1] = (sqrt (err[1])) / (typ) size;
	printf ("\nerrs between %s=%10.4e %10.4e %10.4e",str,err[0],err[1], err[2]);
 	return err[1];
}

