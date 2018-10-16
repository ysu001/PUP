/* 2003.05.07 - module to do Parzen Intensity Inhomogeneity Correction (PIIC).*/
/*$Header: /data/petsun4/data1/src_solaris/nonlin/RCS/piic.c,v 1.5 2004/11/03 20:48:19 rsachs Exp $*/
/*$Log: piic.c,v $
 * Revision 1.5  2004/11/03  20:48:19  rsachs
 * Aesthetics.
 *
 * Revision 1.4  2004/10/28  21:57:27  rsachs
 * Removed mtg++.
 *
 * Revision 1.3  2004/07/21  19:30:05  rsachs
 * Minor formatting/cosmetic changes.
 *
 * Revision 1.2  2003/12/05  20:10:21  rsachs
 * Reduced #piic arguments by 1 (No more iter).
 *
 * Revision 1.1  2003/07/02  18:41:55  rsachs
 * Initial revision
 *
 
**/
 
 #include <stdlib.h>
 #include <stdio.h>
 #include <math.h>
 typedef float typ; 
 #define TWOPI 6.283185307179586
 #define Debug

 extern void ludcmp (typ **a, int n, int *mdx, typ *d);
 extern void lubksb (typ **a, int n, int *mdx, typ *b); 
 extern void parzen_ (typ *img, int *nx, int *ny, int *nz, typ *voxdim, typ *fwhm);
 extern void imgpac_ (typ *img, int *nx, int *ny, int *nz, typ *imgp, int *nxp, int *nyp, int *nzp);
 extern void imgdap_ (typ *img, int *nx, int *ny, int *nz, typ *imgp, int *nxp, int *nyp, int *nzp);
 extern void relax (typ *img, int *imgdim, typ *voxdim, short *mask0, short *mask1,
                    typ *factor, typ *crit);
 
/*----------------------------------------------------------------------
 *  void piic () - function to do Parzen iic.
 *  Definition of Arguments:
 *  Input Arguments:
 *
 *    img0     - target image. (Not modified). 
 *    img1     - source image that need be corrected.  (Not modified).
 *    bb       - a 4 vector rendering the linear gradient.
 *    mask0    - mask for img0.
 *    mask1    - mask for img1.
 *    voxdim   - voxel dimensions vector. Has 3 components (x,y,z). 
 *    nx,ny,nz - nr. of voxels in x,y,z directions respectively.
 *    step     - used in function hcomp.       Usually, step = 1.
 *    factor   - used by function relax.
 *    crit     - used by function relax.
 *    fwhm     - full width half max (mm).
 * 
 *  Output Arguments:
 *
 *    imgh     - volume resulting from affine transformation.
 *    imgf     - volume resulting from Parzen filtering. 
 *    img2     - output image. Should look 'like' img0.
 *    imgd     - difference image.
 * 
 *  Argument count:
 *    10 typ * (6 imgs + 4 more typ * (bb,voxdim,factor,crit))
 *     1 typ   (fwhm)
 *     2 short *
 *     4 int   (nx,ny,nz,step)
 *  Total   = 17 arguments. 
 *
 *  Important Note: All argument images passed to this function are
 *  padded. And so are the masks.
 *--------------------------------------------------------------------*/
  void piic (typ *img0, typ *img1, typ *img2, typ *imgf, typ *imgh, typ 
             *imgd, typ *bb, short *mask0, short *mask1, float *voxdim, int nx, 
              int ny, int nz, int step, typ *factor, typ *crit, typ fwhm)
{
	static int j, k, size, ndx, M, imgdim[3], margin, mtg = 0;
	static typ err[3], Lx, Ly, Lz, nr[3], wrap_flg = 0, val, f0;

 	size = nx * ny * nz;              /* initialising the h & f volumes */
	for (ndx = 0; ndx < size; ndx++)         
		imgf[ndx] = imgh[ndx] = 1.0;

	imgerr (img0,img1,imgf,imgh,imgd,mask0,mask1,nx,ny,nz,err,"0 & 1");
	norms (imgd,nr,mask0,size," imgdi "); 

                            /*               img0 / (h * img1)    */
	inorms (mask0,(long)size,"\nnorms of mask0(b4 hcomp): ");
	inorms (mask1,(long)size,"\nnorms of mask1(b4 hcomp): ");
	norms (img0,nr,mask0,size,"img0(b4 hcomp): ");
	norms (img1,nr,mask0,size,"img1(b4 hcomp): ");
	norms (imgf,nr,mask0,size,"imgf(b4 hcomp): ");
	norms (imgh,nr,mask0,size,"imgh(b4 hcomp): ");
	hcomp (img0,img1,imgf,imgh,mask0,mask1,voxdim,bb,nx,ny,nz,step);
	for (ndx = 0; ndx < size; ndx++) {       /*  img0 / (h * img1)   */
		if (!(mask0[ndx] && mask1[ndx])) imgf[ndx] = 1.0;
        	else if (imgh[ndx] * img1[ndx] != 0.0)
                	imgf[ndx] = img0[ndx] / (imgh[ndx] * img1[ndx]);
        	else imgf[ndx] = 1.0;
   	}
	norms (imgf,nr,mask0,size," imgf (intermediate) ");
	imgdim[0] = nx; imgdim[1] = ny; imgdim[2] = nz;
	relax (imgf, imgdim, voxdim, mask0, mask1, factor, crit); 
	M = fwhm / .7223517245; margin = 2 * M + 1;
	parzen_ (imgf,&nx,&ny,&nz,voxdim,&fwhm);   
	for (k = 0; k < size; k++) img2[k] = imgf[k] * imgh[k] * img1[k]; 
	norms (img2,nr,mask0,size," img2 "); 
	imgerr (img0,img1,imgf,imgh,imgd,mask0,mask1,nx,ny,nz,err,"0&1 final");
	norms (imgd,nr,mask0,size," imgdf "); 
/*	mtg++; printf("\npiic: mtg = %d\n", mtg);          */
}                                          /*  End piic()   */

void relax (typ *imag, int *imgdim, float *voxdim, short *mask0, short *mask1, 
            typ *factor, typ *crit) {
	static	int		i, j, k, dimension;
	static	int		ix, iy, iz, kndex, jndex, index, nvox, iter;
	static	double		del, sum6, q, w[3], err2, err, err0, amean;

	dimension = imgdim[0]*imgdim[1]*imgdim[2];
/*********************************************/
/* compute within-mask mean and count voxels */
/*********************************************/
	for (amean = nvox = i = 0; i < dimension; i++) if ( mask0[i] && mask1[i]) {
		amean += imag[i];
		nvox++;
	}
	amean /= nvox;
	if (nvox == dimension) return;

/************************************************/
/* compute voxdim dependent directional weights */
/************************************************/
	for (q = k = 0; k < 3; k++) {
		w[k] = 1./(voxdim[k]*voxdim[k]);
		q += w[k];
	}
	for (k = 0; k < 3; k++) w[k] *= 0.5/q;

/**********************/
/* relax to criterion */
/**********************/
	iter = 0; err0 = 0.0; do {
		err2 = 0;
		for (iz = 1; iz < imgdim[2]-1; iz++) {kndex = iz*imgdim[0]*imgdim[1];
		for (iy = 1; iy < imgdim[1]-1; iy++) {jndex = iy*imgdim[0] + kndex;
		for (ix = 1; ix < imgdim[0]-1; ix++) {index = ix + jndex;
			if ( mask0[index] && mask1[index] ) continue;
			sum6 =   w[0]*(imag[index+1]+imag[index-1])
		    		+w[1]*(imag[index+imgdim[0]]+imag[index-imgdim[0]])
		    		+w[2]*(imag[index+imgdim[0]*imgdim[1]]+imag[index-imgdim[0]*imgdim[1]]);
			del = sum6 - imag[index];
			imag[index] += *factor*del;
			err2 += del*del;
		}}}
		err = sqrt (err2 / (dimension-nvox));
		iter++;
/*		printf ("iter %4d rm error %.6f Dlnerr/Diter %.6f\n", iter, err/amean, (err-err0)/err0);  */
		err0 = err;
	} while (err > *crit*amean);
}                                     /* End relax */
