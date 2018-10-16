#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <knee.h>

MOMENT blobs[MAXBLOB];
int    adj_reg[MAXBLOB], parents[MAXBLOB];

/***** feature processes *****/
/*
** name - find_blob
** use  - finds connected regions in segmented image 
** input - src: pointer to source image
**         dst: pointer to destination image
**         npixel: number of pixels in image
** output- returns in dst the map of the new connected regions with
**         the new coding. stat information is available in objects
**         array, and raw info is in blobs, adj_reg and parents arrays
**         returns as an integer the number of connected regions
** effect- modify various arrays as indicated in output
** note  - finds connected blobs by scanning the image using pattern:
**  c b a      p is the current pixel, a-d is the scan sequence,
**  d p e      each a-d are checked first for a numbered region 
** that would indicate connection to a previous region. If none then
** create a new one in the last else. adjacent regions are stored
** in the array adj_reg which is manipulated by connect_reg. 
** position e is not necessary but can save time in certain situtations.
** "a" is a special case since can create connections between regions.
** foreground is 8-connected, background is 4-connected.
** works for images that have a frame of background.
*/

int
find_blob(short int *src, short int *dst, int ncol, int npix)
{
    int i, nblob=0, nblob1, *lub;
    short int *bimg,*ba,*bb,*bc,*bd,*be;

/* allocate and copy original image */
    bimg = (short int *) malloc(npix*sizeof(short int));
    for (i=0; i<npix; ++i) *(bimg+i) = *(src+i);  
/* fill frame as bg. this does not add the area of bg correctly */
    nblob = init_moments(bimg,0,nblob,0,0,ncol);   
    for (i=0; i<=ncol; ++i) *(bimg+i) = nblob;  
/* set pointers */
    ba = bimg+2; bb = bimg+1;
    bc = bimg;   bd = bimg+ncol;
    be = bimg+ncol+2;  /* lead pixel, used to save some time and memory */

    for (i=ncol+1; i<npix; ++i) {
        if (bimg[i]) {      /* foreground */
            if (blobs[*ba].fgflg==1) { add_moments(bimg,i,*ba,ncol);
                if (blobs[*bd].fgflg==1 && *ba != *bd) connect_reg(*ba,*bd);
              else if (blobs[*bc].fgflg==1 && *ba != *bc) connect_reg(*ba,*bc);
            }
            else if (blobs[*bb].fgflg==1) add_moments(bimg,i,*bb,ncol);
            else if (blobs[*bc].fgflg==1) add_moments(bimg,i,*bc,ncol);
            else if (blobs[*bd].fgflg==1) add_moments(bimg,i,*bd,ncol);
            else nblob = init_moments(bimg,i,nblob,*bd,1,ncol);
        }
        else {    /* background */
            if (blobs[*bb].fgflg==0) { add_moments(bimg,i,*bb,ncol);
                if (blobs[*bd].fgflg==0 && *bb != *bd) connect_reg(*bb,*bd);
            }
            else if (blobs[*bd].fgflg==0) add_moments(bimg,i,*bd,ncol);
            else if (blobs[*ba].fgflg==0 && *be==0) 
		add_moments(bimg,i,*ba,ncol);
            else nblob = init_moments(bimg,i,nblob,*bd,0,ncol);
        }
        ++ba; ++bb; ++bc; ++bd; ++be;
    } 
    if (nblob >= MAXBLOB-1) {
	printf("Error: num raw blobs too large\n");
	return(-1);
    }

/* lub converts from raw blobs to final blobs */
    lub = (int *) malloc((nblob+1)*sizeof(int));
    /* keep holes with the blobs */
    /* nblob1 = calc_blobs(nblob, lub); */
    /* fill up holes in blobs */
    nblob1 = calc_blobs_fill(nblob, lub);
    /*
    for (i=0; i<nblob+1; i++) printf("lub %d: %d\n",i,lub[i]); 
    */
    for (i=0; i<npix; ++i) *(dst+i) = lub[*(bimg+i)];

    free(lub);
    free(bimg);
    return(nblob1);
}


/*
** name - find_blob_2d
** use - like find_blob except the search is done d,c,b,a (reverse order)
**       this may save time in certain types of image
*/
int
find_blob_2d(short int *src, short int *dst, int npixel, int ncol)
{
    int i, nblob=0, *lub;
    short int *bimg,*ba,*bb,*bc,*bd,*be;

/* allocate and copy original image */
    bimg = (short int *) malloc(npixel*sizeof(short int));
    for (i=0; i<npixel; ++i,++src) bimg[i] = *src;  
/* fill 1st row as bg. this does not add the area of bg correctly */
    nblob = init_moments(bimg,0,nblob,0,0,ncol);
    for (i=0; i<=DY; ++i) bimg[i] = nblob;  
/* set pointers */
    ba = bimg+2; bb = bimg+1;
    bc = bimg; bd = bimg+DY;
    be = bimg+DY+2;  /* lead pixel, used to save some time and memory */

    for (i=DY+1; i<npixel; ++i) {
        if (bimg[i]) {      /* foreground */
            if (blobs[*bd].fgflg==1) { add_moments(bimg,i,*bd, ncol);
                if (blobs[*ba].fgflg==1 && *bd != *ba) connect_reg(*bd,*ba);
            }
            else if (blobs[*bc].fgflg==1) { add_moments(bimg,i,*bc, ncol);
                if (blobs[*ba].fgflg==1 && *bc != *ba) connect_reg(*bc,*ba);
            }
            else if (blobs[*bb].fgflg==1) add_moments(bimg,i,*bb, ncol);
            else if (blobs[*ba].fgflg==1) add_moments(bimg,i,*ba, ncol);
            else nblob = init_moments(bimg,i,nblob,*bd,1,ncol);
        }
        else {    /* background */
            if (blobs[*bd].fgflg==0) { add_moments(bimg,i,*bd, ncol);
                if (blobs[*bb].fgflg==0 && *bd != *bb) connect_reg(*bd,*bb);
            }
            else if (blobs[*bb].fgflg==0) add_moments(bimg,i,*bb, ncol);
            else if (blobs[*ba].fgflg==0 && *be==0) add_moments(bimg,i,*ba, ncol);
            else nblob = init_moments(bimg,i,nblob,*bd,0,ncol);
        }
        ++ba; ++bb; ++bc; ++bd; ++be;
    }
    printf("num raw blobs = %d\n", nblob);

/* lub converts from raw blobs to final blobs */
    lub = (int *) malloc((nblob+1)*sizeof(int));
    nblob = calc_blobs(nblob, lub);
    for (i=0; i<npixel; ++i, ++dst) *dst = lub[bimg[i]];
    free(lub);
    free(bimg);
    return(nblob);
}

/*
** name - find_blob_rv
** use - reverse the 4,8 connectivity, foreground is 4, background is 8
*/
int
find_blob_rv(short int *src,short int *dst,int npixel,int ncol)
{
    int i, nblob=0, *lub;
    short int *bimg,*ba,*bb,*bc,*bd,*be;

/* allocate and copy original image */
    bimg = (short int *) malloc(npixel*sizeof(short int));
    for (i=0; i<npixel; ++i,++src) bimg[i] = *src;  
/* fill 1st row as bg. this does not add the area of bg correctly */
    nblob = init_moments(bimg,0,nblob,0,0,ncol);
    for (i=0; i<=DY; ++i) bimg[i] = nblob;  
/* set pointers */
    ba = bimg+2; bb = bimg+1;
    bc = bimg; bd = bimg+DY;
    be = bimg+DY+2;  /* lead pixel, used to save some time and memory */

    for (i=DY+1; i<npixel; ++i) {
        if (!bimg[i]) {      /* background */
            if (blobs[*bd].fgflg==0) { add_moments(bimg,i,*bd, ncol);
                if (blobs[*ba].fgflg==0 && *bd != *ba) connect_reg(*bd,*ba);
            }
            else if (blobs[*bc].fgflg==0) { add_moments(bimg,i,*bc, ncol);
                if (blobs[*ba].fgflg==0 && *bc != *ba) connect_reg(*bc,*ba);
            }
            else if (blobs[*bb].fgflg==0) add_moments(bimg,i,*bb, ncol);
            else if (blobs[*ba].fgflg==0) add_moments(bimg,i,*ba, ncol);
            else nblob = init_moments(bimg,i,nblob,*bd,0,ncol);
        }
        else {    /* foreground */
            if (blobs[*bd].fgflg==1) { add_moments(bimg,i,*bd, ncol);
                if (blobs[*bb].fgflg==1 && *bd != *bb) connect_reg(*bd,*bb);
            }
            else if (blobs[*bb].fgflg==1) add_moments(bimg,i,*bb, ncol);
            else if (blobs[*ba].fgflg==1 && *be==1) add_moments(bimg,i,*ba, ncol);
            else nblob = init_moments(bimg,i,nblob,*bd,1,ncol);
        }
        ++ba; ++bb; ++bc; ++bd; ++be;
    }
    printf("num raw blobs = %d\n", nblob);

/* lub converts from raw blobs to final blobs */
    lub = (int *) malloc((nblob+1)*sizeof(int));
    nblob = calc_blobs(nblob, lub);
    for (i=0; i<npixel; ++i, ++dst) *dst = lub[bimg[i]];
    free(lub);
    free(bimg);
    return(nblob);
}

/*
** name - find_blob_3d 
** use - finds connected blobs in 3d image using the scan sequence:
**  e f g               c b a   p is the point a-m is scan sequence
**  h i j  (top layer)  d p
**  k l m 
** foreground is 26-connected, background is 4-connected.
** direct extension of the basic find_blob
** works for images that have a frame of background.
*/
int
find_blob_3d(short int *src, short int *dst, int npixel, int ncol)
{
    int i, nblob=0, *lub;
    short int *bimg,*ba,*bb,*bc,*bd,*be,*bf,*bg,*bh,*bi,*bj,*bk,*bl,*bm;

/* allocate and copy original image */
    bimg = (short int *) malloc(npixel*sizeof(short int));
    for (i=0; i<npixel; ++i,++src) bimg[i] = *src;  
/* fill 1st plane as bg. this does not add the area of bg correctly */
    nblob = init_moments(bimg,0,nblob,0,0,ncol);
    for (i=0; i<=DZ; ++i) bimg[i] = nblob;  
/* set pointers */
    ba = bimg+IND(2,0,1); bb = bimg+IND(1,0,1); 
    bc = bimg+IND(0,0,1); bd = bimg+IND(0,1,1); 
    be = bimg+IND(0,0,0); bf = bimg+IND(1,0,0); bg = bimg+IND(2,0,0);
    bh = bimg+IND(0,1,0); bi = bimg+IND(1,1,0); bj = bimg+IND(2,1,0);
    bk = bimg+IND(0,1,0); bl = bimg+IND(1,1,0); bm = bimg+IND(2,1,0);

    for (i=DZ+DY+1; i<npixel; ++i) {
        if (bimg[i]) {      /* foreground */
            if (blobs[*ba].fgflg==1) { add_moments(bimg,i,*ba, ncol);
                if (blobs[*bc].fgflg==1 && *ba != *bc) connect_reg(*ba,*bc);
                if (blobs[*bd].fgflg==1 && *ba != *bd) connect_reg(*ba,*bd); 
                if (blobs[*be].fgflg==1 && *ba != *be) connect_reg(*ba,*be); 
                if (blobs[*bh].fgflg==1 && *ba != *bh) connect_reg(*ba,*bh); 
                if (blobs[*bk].fgflg==1 && *ba != *bk) connect_reg(*ba,*bk); 
                if (blobs[*bl].fgflg==1 && *ba != *bl) connect_reg(*ba,*bl); 
                if (blobs[*bm].fgflg==1 && *ba != *bm) connect_reg(*ba,*bm); 
            }
            else if (blobs[*bb].fgflg==1) { add_moments(bimg,i,*bb, ncol);
                if (blobs[*bk].fgflg==1 && *bb != *bk) connect_reg(*bb,*bk); 
                if (blobs[*bl].fgflg==1 && *bb != *bl) connect_reg(*bb,*bl); 
                if (blobs[*bm].fgflg==1 && *bb != *bm) connect_reg(*bb,*bm); 
            }
            else if (blobs[*bc].fgflg==1) { add_moments(bimg,i,*bc, ncol);
                if (blobs[*bg].fgflg==1 && *bc != *bg) connect_reg(*bc,*bg); 
                if (blobs[*bj].fgflg==1 && *bc != *bj) connect_reg(*bc,*bj); 
                if (blobs[*bk].fgflg==1 && *bc != *bk) connect_reg(*bc,*bk); 
                if (blobs[*bl].fgflg==1 && *bc != *bl) connect_reg(*bc,*bl); 
                if (blobs[*bm].fgflg==1 && *bc != *bm) connect_reg(*bc,*bm); 
            }
            else if (blobs[*bd].fgflg==1) { add_moments(bimg,i,*bd, ncol);
                if (blobs[*bg].fgflg==1 && *bd != *bg) connect_reg(*bd,*bg); 
                if (blobs[*bj].fgflg==1 && *bd != *bj) connect_reg(*bd,*bj); 
                if (blobs[*bm].fgflg==1 && *bd != *bm) connect_reg(*bd,*bm); 
            }
            else if (blobs[*be].fgflg==1) { add_moments(bimg,i,*be, ncol);
                if (blobs[*bg].fgflg==1 && *be != *bg) connect_reg(*be,*bg); 
                if (blobs[*bj].fgflg==1 && *be != *bj) connect_reg(*be,*bj); 
                if (blobs[*bk].fgflg==1 && *be != *bk) connect_reg(*be,*bk); 
                if (blobs[*bl].fgflg==1 && *be != *bl) connect_reg(*be,*bl); 
                if (blobs[*bm].fgflg==1 && *be != *bm) connect_reg(*be,*bm); 
            }
            else if (blobs[*bf].fgflg==1) { add_moments(bimg,i,*bf, ncol);
                if (blobs[*bk].fgflg==1 && *bf != *bk) connect_reg(*bf,*bk); 
                if (blobs[*bl].fgflg==1 && *bf != *bl) connect_reg(*bf,*bl); 
                if (blobs[*bm].fgflg==1 && *bf != *bm) connect_reg(*bf,*bm); 
            }
            else if (blobs[*bg].fgflg==1) { add_moments(bimg,i,*bg, ncol);
                if (blobs[*bh].fgflg==1 && *bg != *bh) connect_reg(*bg,*bh); 
                if (blobs[*bk].fgflg==1 && *bg != *bk) connect_reg(*bg,*bk); 
                if (blobs[*bl].fgflg==1 && *bg != *bl) connect_reg(*bg,*bl); 
                if (blobs[*bm].fgflg==1 && *bg != *bm) connect_reg(*bg,*bm); 
            }
            else if (blobs[*bh].fgflg==1) { add_moments(bimg,i,*bh, ncol);
                if (blobs[*bj].fgflg==1 && *bh != *bj) connect_reg(*bh,*bj); 
                if (blobs[*bm].fgflg==1 && *bh != *bm) connect_reg(*bh,*bm); 
            }
            else if (blobs[*bi].fgflg==1) add_moments(bimg,i,*bi, ncol);
            else if (blobs[*bj].fgflg==1) { add_moments(bimg,i,*bj, ncol);
                if (blobs[*bk].fgflg==1 && *bj != *bk) connect_reg(*bj,*bk); 
            }
            else if (blobs[*bk].fgflg==1) { add_moments(bimg,i,*bk, ncol);
                if (blobs[*bm].fgflg==1 && *bk != *bm) connect_reg(*bk,*bm); 
            }
            else if (blobs[*bl].fgflg==1) add_moments(bimg,i,*bl, ncol);
            else if (blobs[*bm].fgflg==1) add_moments(bimg,i,*bm, ncol);
            else nblob = init_moments(bimg,i,nblob,*bd,1,ncol);
        }
        else {    /* background */
            if (blobs[*bb].fgflg==0) { add_moments(bimg,i,*bb, ncol);
                if (blobs[*bd].fgflg==0 && *bb != *bd) connect_reg(*bb,*bd);
                if (blobs[*bi].fgflg==0 && *bb != *bi) connect_reg(*bb,*bi);
            }
            else if (blobs[*bd].fgflg==0) { add_moments(bimg,i,*bd, ncol);
                if (blobs[*bi].fgflg==0 && *bd != *bi) connect_reg(*bd,*bi);
            }
            else nblob = init_moments(bimg,i,nblob,*bd,0,ncol);
        }
        ++ba; ++bb; ++bc; ++bd; 
        ++be; ++bf; ++bg; 
        ++bh; ++bi; ++bj; 
        ++bk; ++bl; ++bm; 
    }
    printf("num raw blobs = %d\n", nblob);

/* lub converts from raw blobs to final blobs */
    lub = (int *) malloc((nblob+1)*sizeof(int));
    nblob = calc_blobs(nblob, lub);
    for (i=0; i<npixel; ++i, ++dst) *dst = lub[bimg[i]];
    free(lub);
    free(bimg);
    return(nblob);
}

/*
** name - init_moments
** use  - allocate memory for new blob in blobs array
** input - bimg: blob image array
**         i: first pixel of new blob
**         nblob: last blob number, will be incremented and returned
**         nprnt: number of parent of new blob
**         fgflg: 1=forground or 0=background flag
**	   ncol: number pixels in a row=number of columns
** output- returns new blob number
**         if end of memory, print error and return last nblob
** effect- fill the blob information to track new blob 
** note  - can be used for 3D but z has to be set 
*/
int
init_moments(short int *bimg,int i,int nblob,int nprnt,int fgflg,int ncol)
{
    MOMENT *blob;
    int x=i%ncol, y=i/ncol, z = 0;

    ++nblob; 
    if (nblob >= MAXBLOB) {
        /* printf("Error - number raw blobs too large\n"); */
        return(nblob-1);
    }
    bimg[i] = nblob; 
    adj_reg[nblob] = nblob;
    parents[nblob] = nprnt;
    blob = &blobs[nblob];
    blob->m000 = 1;
    blob->m100 = blob->xbeg = x;
    blob->m010 = blob->ybeg = y;
    blob->m001 = blob->zbeg = z;
    blob->m110 = (double) x*y;
    blob->m101 = (double) x*z;
    blob->m011 = (double) y*z;
    blob->m200 = (double) x*x;
    blob->m020 = (double) y*y;
    blob->m002 = (double) z*z;
    blob->fgflg = fgflg;
    return(nblob);
}

/*
** name - add_moments
** use  - add new pixel to existing blob and keep track of stats
** input - bimg: blob image array
**         i: new pixel of blob
**         nblob: blob that gets new pixel
** output- none
** effect- add the blob information to track stats
** note  - can be used for 3D but z has to be set 
*/
void
add_moments(short int *bimg, int i, int nblob, int ncol)
{
    MOMENT *blob;
    int x=i%ncol, y=i/ncol, z = 0;

    bimg[i] = nblob;
    blob = &blobs[nblob];
    blob->m000 += 1;
    blob->m100 += x;
    blob->m010 += y;
    blob->m001 += z;
    blob->m110 += (double) x*y;
    blob->m101 += (double) x*z;
    blob->m011 += (double) y*z;
    blob->m200 += (double) x*x;
    blob->m020 += (double) y*y;
    blob->m002 += (double) z*z;
}

/*
** name - connect_reg
** use  - make a note of regions/blobs that are adjacent/same
** input - n,m: the two blob numbers that are connected and 
**         belong to the same region and will be consolidated
** output- none
** effect- make the change in the adj_reg array to track this info
** note  - good for 2 and 3 dimensions
*/
void
connect_reg(int n, int m)
{
    int tmp;

    tmp = adj_reg[n];
    while (tmp != n) {   /* check if connected */
        if (tmp==m) return;
        tmp = adj_reg[tmp];
    }
    tmp = adj_reg[n];    /* no so connect */
    adj_reg[n] = adj_reg[m];
    adj_reg[m] = tmp;
}

/*
** name - calc_blobs
** use  - convert from raw blobs to real blobs
** input - nblob: final number of raw blobs
**         lub: table to point from raw blob to new blob number
** output- returns the number of real blobs
** effect- consolidates stats for real blobs in new allocated area
**         assigned to objects global array
** note  - good for 2 and 3 dimensions
**         works in three steps: consolidates raw info,
**         applies filters to eliminate noise, then tags on
**         daughter holes as a linked list to each object
*/
int
calc_blobs(int nblob, int *lub)
{
    int     i,j,tmp,next,nhole, nobj=0;
    MOMENT  *blob;
    PBLOB   obj;

/* consolidate blobs */
    for (i=1; i<=nblob; ++i) {
        if (adj_reg[i]==0) continue;
        blob = &blobs[i];
        lub[i] = i;
        next = adj_reg[i];
        while (next != i) {
            blob->m000 += blobs[next].m000;
            blob->m100 += blobs[next].m100;
            blob->m010 += blobs[next].m010;
            blob->m001 += blobs[next].m001;
            blob->m110 += blobs[next].m110;
            blob->m101 += blobs[next].m101;
            blob->m011 += blobs[next].m011;
            blob->m200 += blobs[next].m200;
            blob->m020 += blobs[next].m020;
            blob->m002 += blobs[next].m002;
            tmp = adj_reg[next];
            lub[next] = i;
            adj_reg[next] = 0;
            for (j=0;j<=nblob;++j) if (parents[j]==next) parents[j]=i;
            next = tmp;
        }
/* if foreground and greater than min area keep as an object */
        if (blob->fgflg==1 && blob->m000 >= AREA_MIN) {
            objects[nobj] = obj = alloc_blob();
            obj->xbeg = blob->xbeg;
            obj->ybeg = blob->ybeg;
            obj->zbeg = blob->zbeg;
            obj->m000  = blob->m000;
            obj->m100  = blob->m100;
            obj->m010  = blob->m010;
            obj->m001  = blob->m001;
            obj->m110  = blob->m110;
            obj->m101  = blob->m101;
            obj->m011  = blob->m011;
            obj->m200  = blob->m200;
            obj->m020  = blob->m020;
            obj->m002  = blob->m002;
            obj->nblob = i;
            obj->parent = parents[i];
            ++nobj;
            if (nobj >= MAXOBJ) return(nobj);
        }
    }

/* add holes to their parents */
    for (nobj=0, i=1; i<=nblob; ++i) {
/* continue only if real object */
        if (adj_reg[i]==0 || blobs[i].m000 < AREA_MIN ||
            blobs[i].fgflg==0) continue;
/* scan for holes that belong to object */
        for (nhole=0, j=1; j<=nblob; ++j) {
/* check for all conditions for hole */
            if (adj_reg[j]!=0 && blobs[j].m000 >= HOLE_MIN &&
                blobs[j].fgflg==0 && i!=j && parents[j]==i) {
                blob = &blobs[j];
                ++nhole;
                obj = alloc_blob();
                obj->xbeg = blob->xbeg;
                obj->ybeg = blob->ybeg;
                obj->zbeg = blob->zbeg;
                obj->m000  = blob->m000;
                obj->m100  = blob->m100;
                obj->m010  = blob->m010;
                obj->m001  = blob->m001;
                obj->m110  = blob->m110;
                obj->m101  = blob->m101;
                obj->m011  = blob->m011;
                obj->m200  = blob->m200;
                obj->m020  = blob->m020;
                obj->m002  = blob->m002;
                obj->nblob = j;
                obj->parent = parents[j];
                obj->next = objects[nobj]->hole;
                objects[nobj]->hole = obj;
            }
        }
        objects[nobj]->nholes = nhole;
        ++nobj;
    }
    return(nobj);
}


/*
** name - calc_blobs_fill
** use  - convert from raw blobs to real blobs filling holes
** input - nblob: final number of raw blobs
**         lub: table to point from raw blob to new blob number
** output- returns the number of real blobs
** effect- consolidates stats for real blobs in new allocated area
**         assigned to objects global array
** note  - good for 2 and 3 dimensions
**         works in three steps: consolidates raw info,
**         applies filters to eliminate noise, then fills holes
*/
int
calc_blobs_fill(int nblob, int *lub)
{
    int     i,j,k,tmp,next,nhole, nobj=0;
    MOMENT  *blob;
    PBLOB   obj;

/* consolidate blobs */
    for (i=1; i<=nblob; ++i) {
        if (adj_reg[i]==0) continue;
        blob = &blobs[i];
        lub[i] = i;
        next = adj_reg[i];
        while (next != i) {
            blob->m000 += blobs[next].m000;
            blob->m100 += blobs[next].m100;
            blob->m010 += blobs[next].m010;
            blob->m001 += blobs[next].m001;
            blob->m110 += blobs[next].m110;
            blob->m101 += blobs[next].m101;
            blob->m011 += blobs[next].m011;
            blob->m200 += blobs[next].m200;
            blob->m020 += blobs[next].m020;
            blob->m002 += blobs[next].m002;
            tmp = adj_reg[next];
            lub[next] = i;
            adj_reg[next] = 0;
            for (j=0;j<=nblob;++j) if (parents[j]==next) parents[j]=i;
            next = tmp;
        }

/* if foreground and greater than min area keep as an object */
        if (blob->fgflg==1 && blob->m000 >= AREA_MIN) {
            objects[nobj] = obj = alloc_blob();
            obj->xbeg = blob->xbeg;
            obj->ybeg = blob->ybeg;
            obj->zbeg = blob->zbeg;
            obj->m000  = blob->m000;
            obj->m100  = blob->m100;
            obj->m010  = blob->m010;
            obj->m001  = blob->m001;
            obj->m110  = blob->m110;
            obj->m101  = blob->m101;
            obj->m011  = blob->m011;
            obj->m200  = blob->m200;
            obj->m020  = blob->m020;
            obj->m002  = blob->m002;
            obj->nblob = i;
            obj->parent = parents[i];
            ++nobj;
            if (nobj >= MAXOBJ) return(nobj);
        }
    }

/* add holes to their parents */
    for (nobj=0, i=1; i<=nblob; ++i) {
/* continue only if real object */
        if (adj_reg[i]==0 || blobs[i].m000 < AREA_MIN ||
            blobs[i].fgflg==0) continue;
/*
printf("real obj %d area %d\n", i, blobs[i].m000);
*/
/* scan for holes that belong to object */
        for (nhole=0, j=1; j<=nblob; ++j) {
/* check for all conditions for hole */
            if (adj_reg[j]!=0 && blobs[j].m000 >= HOLE_MIN &&
                blobs[j].fgflg==0 && i!=j && parents[j]==i) {
                blob = &blobs[j];
                ++nhole;
                obj = objects[nobj];
                obj->m000  += blob->m000;
                obj->m100  += blob->m100;
                obj->m010  += blob->m010;
                obj->m001  += blob->m001;
                obj->m110  += blob->m110;
                obj->m101  += blob->m101;
                obj->m011  += blob->m011;
                obj->m200  += blob->m200;
                obj->m020  += blob->m020;
                obj->m002  += blob->m002;
		/* hole filler */
		lub[j] = i; 
                for (k=0;k<=nblob;k++) if (lub[k]==j) lub[k] = i;
		/*
		printf("    hole %d area %d xy %d %d \n", j, blob->m000, blob->xbeg, blob->ybeg);
		*/
            }
        }
        ++nobj;
    }
    return(nobj);
} 


/*
** name - alloc_blob
** use  - allocate memory for a blob
** input - none
** output- returns the address of allocated memory
** effect- allocate memory, error message if can't
** note  - good for 2 and 3 dimensions
*/
PBLOB
alloc_blob()
{
    PBLOB blob;

    if ((blob=(PBLOB)malloc(sizeof(BLOB))) != (PBLOB)NULL) {
        blob->hole = (PBLOB) NULL;
        blob->bay  = (PBLOB) NULL;
        blob->next = (PBLOB) NULL;
        blob->perim = (PPOINT) NULL;
        blob->con_hull = (PPOINT) NULL;
    }
    else printf("Error - blob memory allocation\n");
    return(blob);
}

/*
** centroid - calculate rotated system
*/
void
centroid(int iobj)
{
    double num,xcm,ycm,zcm,xycm,xzcm,yzcm,x2cm,y2cm,z2cm;
    double xcen,ycen,zcen,mx2,my2,mxy,a1,a2,c,mbov2,eg1,eg2,eccen;

    num = (double) objects[iobj]->m000;
    xcm = (double) objects[iobj]->m100;
    ycm = (double) objects[iobj]->m010;
    zcm = (double) objects[iobj]->m100;
    xycm = objects[iobj]->m110;
    xzcm = objects[iobj]->m101;
    yzcm = objects[iobj]->m011;
    x2cm = objects[iobj]->m200;
    y2cm = objects[iobj]->m020;
    z2cm = objects[iobj]->m002;

    xcen = xcm/num; ycen = ycm/num; zcen = zcm/num;
    mx2 = x2cm - (xcm*xcm/num);
    my2 = y2cm - (ycm*ycm/num);
    mxy = xycm - (xcm*ycm/num);
/*
    my2 *= PIXRAT*PIXRAT; mxy *= PIXRAT;
*/
    a1 = atan2(mxy, (mx2-my2)/2.0)/2.0;
    mbov2 = (mx2+my2)/2.0;
    c = mx2*my2 - mxy*mxy;
    eg1 = mbov2 + sqrt(mbov2*mbov2 - c);
    eg2 = mbov2 - sqrt(mbov2*mbov2 - c);
    eccen = eg2/eg1;
/*
    a1 = atan2(eg1-mx2, mxy);
    a2 = atan2(eg2-mx2, mxy);
*/

/* to find min and max in rotated system
    xmin = ymin = 512;
    xmax = ymax = 0;
loop over all perimeter points {
        x0 = PROTX(x,y,ca,sa);
        y0 = PROTY(x,y,ca,sa);
        if (x0 > xmax) xmax = x0;
        if (x0 < xmin) xmin = x0;
        if (y0 > ymax) ymax = y0;
        if (y0 < ymin) ymin = y0;
    }
*/
}

/*
** name - perim_trace 
** use  - trace the perimeter of a blob in 2D only
** input - npixel: number of pixels
**         src: source image
**         xbeg,ybeg: start location
**         xp,yp: arrays to store the value of the perimeter points
**         id: what is the id value of the object in the image
**         cw_flg: clockwise or counter cw (counter used for holes)
** output- dst: mask of perimeter
** return- returns the number of points on the perimeter 
**         the points are stored in the arrays xp and yp
** effect- print the circumference. other stats can be gathered but 
**         are already available from the blob routine
** note  - none
*/
int
perim_trace(short int *src, short int *dst, int xbeg, int ybeg,
	int *xp,int *yp,int id, int cw_flg, int ncol,int npix)
{
    int x8n(int d),y8n(int d);
    int i,j,xt,yt,x0,y0,dt,dbeg, npts=0;

/* cw - trace a white object; ccw - trace a black hole; */
    if (cw_flg)  j=7; 
    else { j=3; --xbeg; }

    xt = xbeg; yt = ybeg;
    for (i=0; i<8; ++i) {
        x0 = x8n(MOD8(j+i));
        y0 = y8n(MOD8(j+i));
        if (*(src+(yt+y0)*ncol+xt+x0) == id) {
            dt = dbeg = MOD8(j+i);
            j = MOD8(2*(dt/2)+7);
            break;
        }
    }
    do {
	if (npts>=MAXPTS) printf("Error - npts in perimeter too large\n");
	xp[npts] = xt;
	yp[npts] = yt;
        ++npts;
        xt += x0; yt += y0;
        for (i=0; i<8; ++i) {
            x0 = x8n(MOD8(j+i));
            y0 = y8n(MOD8(j+i));
            if (*(src+(yt+y0)*ncol+xt+x0) == id) {
                dt = MOD8(j+i);
                j = MOD8(2*(dt/2)+7);
                break;
            }
        }
    } while ( xt!=xbeg || yt!=ybeg || dt!=dbeg );

/* fill dst with perimeter */
    if (dst != (short int *)NULL) {
	for (i=0; i<npix; ++i)  *(dst+i)=0;
	for (i=0; i<npts; ++i) 
            *(dst + yp[i]*ncol + xp[i]) = NGRAY-1;
    }
    return(npts);

/* used for the hip
    hough_circ(src,xp,yp,circ,ncol,npix);
*/
}


/*
** name - x8n, y8n
** use  - return the shift in x and y needed for scan of 
**        8-connected neighbors
** input - direction code from 0 to 7
** output- returns the shift in x and y
** effect- none
** note  - none
*/
int
x8n(int d)
{
    switch(d) {
        case(0): return(0); 
        case(1): return(1); 
        case(2): return(1); 
        case(3): return(1); 
        case(4): return(0); 
        case(5): return(-1); 
        case(6): return(-1); 
        case(7): return(-1); 
        default: return(0);
    }
}

int
y8n(int d)
{
    switch(d) {
        case(0): return(-1); 
        case(1): return(-1); 
        case(2): return(0); 
        case(3): return(1); 
        case(4): return(1); 
        case(5): return(1); 
        case(6): return(0); 
        case(7): return(-1); 
        default: return(0);
    }
}


/*
** name - trace_line
** use  - give the pixels that lay on a line between point 1 and 2
** input - x1,y1: point 1
**         x2,y2: point 2
**         xo,yo: output array pointers
** output- return the pixels on line in arrays xo,yo
** effect- none
** note  - Bresenham's integer line algorithm from
**         Newman & Sproul, Interactive Computer Graphics, 1979
*/
void
trace_line(int x1,int y1,int x2,int y2,int *xo,int *yo)
{
    int i, e, dx, dy, ix=1, iy=1;
    int x=x1, y=y1;

    dx = ABSV(x2 - x1);
    dy = ABSV(y2 - y1);
    if (x2 < x1) ix = -1;
    if (y2 < y1) iy = -1;

    if (dx >= dy) {
        e = 2*dy - dx;
        for (i=0; i<=dx; ++i) {
            xo[i] = x; yo[i] = y;
            if (e > 0) { y += iy; e += 2*(dy - dx); }
            else e += 2*dy;
            x += ix;
        }
    }
    else {
        e = 2*dx - dy;
        for (i=0; i<=dy; ++i) {
            xo[i] = x; yo[i] = y;
            if (e > 0) { x += ix; e += 2*(dx - dy); }
            else e += 2*dx;
            y += iy;
        }
    }
}

void
draw_line(short int *dst,int x1,int y1,int x2,int y2,int ncol)
{
    int i, e, dx, dy, ix=1, iy=1;
    int x=x1, y=y1;

    dx = ABSV(x2 - x1);
    dy = ABSV(y2 - y1);
    if (x2 < x1) ix = -1;
    if (y2 < y1) iy = -1;

    if (dx >= dy) {
        e = 2*dy - dx;
        for (i=0; i<=dx; ++i) {
	    *(dst+ncol*y+x)=NGRAY-1;
            if (e > 0) { y += iy; e += 2*(dy - dx); }
            else e += 2*dy;
            x += ix;
        }
    }
    else {
        e = 2*dx - dy;
        for (i=0; i<=dy; ++i) {
	    *(dst+ncol*y+x)=NGRAY-1;
            if (e > 0) { x += ix; e += 2*(dx - dy); }
            else e += 2*dx;
            y += iy;
        }
    }
}

void
draw_circle(short int *dst,int xcen,int ycen,double rad,int ncol,int npix)
{
    int i;
    double cs,sn,n,xo,yo,xn,yn,xc,yc;

    xc = (double)xcen; yc = (double)ycen;
/* pick number of points to be 10% larger */
    n = 2.2*M_PI*rad;

    cs = cos(2.0*M_PI/n);
    sn = sin(2.0*M_PI/n);

    xo = rad;
    yo = 0.0;
    for (i=0; i<(int)n; ++i) {
        if ((int)(yo+yc)*ncol+(int)(xo+xc) > npix) continue;
	*(dst+(int)(yo+yc)*ncol+(int)(xo+xc))=NGRAY-1;
	xn = xo*cs - yo*sn;
	yn = yo*cs + xo*sn;
        xo = xn;
	yo = yn;
    }
}
	
/*
double 
fit_circle(int *xp, int *yp, int npts, double xcen,double ycen,double rad)
{
}
*/

/* 
** name - fit_line
** use  - make a linear fit to point collection and return variance
** input - xp,yp: arrays with the points
**         npts: number of points
** output- the variance of the fit
*/
double
fit_line(int *xp,int *yp,int npts)
{
    int i;
    double x,y,n;
    double del,a,b,var;
    double sx=0.0,sy=0.0,sxx=0.0,sxy=0.0,syy=0.0;

/* if only two points no error */
    if (npts==2) return(0.0);

    n = (double)npts;
/* check if to use x or y as base */
    if (ABSV(xp[npts-1]-xp[0])>=ABSV(yp[npts-1]-yp[0])) {
        for (i=0; i<npts; ++i) {
    	    x = (double)xp[i];
	    y = (double)yp[i];
	    sx += x;
	    sy += y;
	    sxx += x*x;
	    sxy += x*y;
	    syy += y*y;
	}
    }
    else {
        for (i=0; i<npts; ++i) {
    	    y = (double)xp[i];
	    x = (double)yp[i];
	    sx += x;
	    sy += y;
	    sxx += x*x;
	    sxy += x*y;
	    syy += y*y;
	}
    }

    del = n*sxx - sx*sx;
    if (del==0.0) printf("Error - div by zero in fit_line\n");
    a = (sxx*sy - sx*sxy)/del;
    b = (n*sxy - sx*sy)/del;
    var = (syy+n*a*a+sxx*b*b-2.0*(a*sy+b*sxy-a*b*sx))/(n-2);

    return(var);
}

/*
** name- find_line
** use - find straight lines in a general perimeter
** input - xp,yp: array with the points 
**         npts: number of points on the perimeter 
**	   dst: image of lines found
*/
int
find_line(short int *dst,int *xp,int *yp,int npts,
	int minlen,double maxvar,int ncol,int npix)
{
    int i=0,j,len,nline=0;
    double var,ang;

/* good assumption: begin points are not in the middle of a line */
/* scan over the perimeter points */
    while (i+minlen<npts) {
/* fit the points to a longer line */
	var = fit_line(xp+i,yp+i,minlen);
	if (var > maxvar) {
	    ++i; continue;
	}
	len = minlen;
	maxvar = MINV(maxvar,(1.3*var));
	while (var<=maxvar && i+len<npts) {
	    ++len;
	    var = fit_line(xp+i,yp+i,len);
	}
	--len;
/* save the line */
/* get a better angle from fit_line */
	ang=atan2((double)(yp[i+len-1]-yp[i]),
	    (double)(xp[i+len-1]-xp[i]));
	fitline[nline].indx = i;
	fitline[nline].len = len;
	fitline[nline].ang = ang;
	fitline[nline].var = var;
	fitline[nline].x1 = xp[i];
	fitline[nline].y1 = yp[i];
	fitline[nline].x2 = xp[i+len-1];
	fitline[nline].y2 = yp[i+len-1];
	printf("%d i %d len %d ang %f 1: %d %d 2: %d %d\n",
	nline,i,len,ang*180.0/M_PI,xp[i],yp[i],xp[i+len-1],yp[i+len-1]);
/*
	printf("    var %f var/n %e\n", var,var/len);
*/
	++nline;
        i += len;
    }

/* draw if needed */
    if (dst != (short int *)NULL) {
	for (i=0; i<npix; ++i)  *(dst+i)=0;
        for (i=0; i<nline; ++i) {
	    for (j=fitline[i].indx; j<fitline[i].indx+fitline[i].len; ++j) {
		*(dst + ncol*yp[j] + xp[j]) = NGRAY-1;
	    }
	}
    }
    return(nline);
}

/* 
** name - convex_hull
** use  - create the convex hull from a collection of points on the 
**	  perimeter of a blob
** input - dst: destination array for image output
**         xp,yp: arrays with perimeter points of shape
**         npts: number of points in xp,yp
**	   ncol,npix: # columns and pixels in dst
** note - Reference: Shin & Woo, Pattern Recognition V19 N6 P453, 1986.
*/
int
convex_hull(short int *dst,int *xp,int *yp,int npts,int ncol,int npix)
{
    int convex_deff(int *xp,int *yp,int *vp,int *zp,int nhull,int npts);
    int i,p=1,q=1;
    int *zp,*vp;
    int minx, miny, mini;
    int xmin,xmax,ymin,ymax;
    int left(int p,int a,int b,int *xp,int *yp);
    int right(int p,int a,int b,int *xp,int *yp);

    zp = (int *) malloc(npts*sizeof(int));
    vp = (int *) malloc((npts+1)*sizeof(int));

/* find the point with minimal x and y values */
    minx = xp[0];
    miny = yp[0];
    mini = 0;
    xmin=xmax=xp[0];
    ymin=ymax=yp[0];
    for (i=1; i<npts; ++i) {
	if ((yp[i]<miny) || (yp[i]==miny && xp[i]<minx)) {
	    minx = xp[i];
	    miny = yp[i];
	    mini = i;
	}
	if (xp[i]<xmin) xmin=xp[i];
	if (xp[i]>xmax) xmax=xp[i];
	if (yp[i]<ymin) ymin=yp[i];
	if (yp[i]>ymax) ymax=yp[i];
    }

/* set up the vp array indexed from the min point */
    vp[npts] = mini;
    for (i=0; i<npts; ++i) vp[i] = (mini+i)%npts;
    zp[0] = vp[0];
    zp[1] = vp[1];

/* loop to create the hull indexed by zp with p points */
    while (q < npts) {
	if (right(vp[q+1],zp[p-1],zp[p],xp,yp)) {
	    if (right(vp[q+1],vp[q-1],vp[q],xp,yp)) {
		++p; ++q; zp[p]=vp[q];
	    }
	    else while(!left(vp[q+1],zp[p-1],zp[p],xp,yp)) ++q;
	}
	else {
	    while(p>0 && !right(zp[p-1],zp[p],vp[q+1],xp,yp)) --p;
	    ++p; ++q; zp[p]=vp[q];
	}
    }

/* draw the convex hull lines  */
    if (dst != (short int *)NULL) {
	for (i=1; i<=p; ++i) draw_line(dst,xp[zp[i-1]],yp[zp[i-1]],
	xp[zp[i]],yp[zp[i]],ncol);
    }

/* clear and draw the convex hull points 
    for (i=0; i<npix; ++i)  *(dst+i)=0;
    for (i=0; i<=p; ++i){
	*(dst + ncol*yp[zp[i]] + xp[zp[i]]) = NGRAY-1;
	printf("i %d, zp %d, xp %d, yp %d\n",
	i,zp[i],xp[zp[i]],yp[zp[i]]);
    }
*/

    i = convex_deff(xp,yp,vp,zp,p,npts);
    free(vp); free(zp);
    return(i);
}

/*
** name - left, right
** use  - assist functions for finding the convex hull
*/
int 
left(int p,int a,int b,int *xp,int *yp)
{
    int lflg=(yp[p]-yp[a])*(xp[b]-xp[a])-(xp[p]-xp[a])*(yp[b]-yp[a]);

    if (lflg<0) return(1);
    else if (lflg>0) return(0);
    else if 
	((yp[p]-yp[b])*(yp[b]-yp[a])+(xp[p]-xp[b])*(xp[b]-xp[a])>0)
	return(1);
    else return(0);
}

int 
right(int p,int a,int b,int *xp,int *yp)
{
    int rflg=(yp[p]-yp[a])*(xp[b]-xp[a])-(xp[p]-xp[a])*(yp[b]-yp[a]);

    if (rflg>0) return(1);
    else if (rflg<0) return(0);
    else if 
	((yp[p]-yp[b])*(yp[b]-yp[a])+(xp[p]-xp[b])*(xp[b]-xp[a])>0)
	return(0);
    else return(1);
}


/* 
** name - convex_deff
** use  - evaluate the convex defficiencies using the convex hull
*/
int
convex_deff(int *xp,int *yp,int *vp,int *zp,int nhull,int npts)
{
    int right_shft(int p,int a,int b,int *xp,int *yp);
    int i,k,q,qp,circ,area,nopn;
    int ndeff=0;
    double ang;
    PBLOB deff;

/* calculate the convex defficiency */
    qp=0; q=1;
    for (i=0; i<npts; ++i) {
	if (right_shft(vp[i],zp[q-1],zp[q],xp,yp))
	{
	    /* can draw outline here */
	    if (qp==q) continue;
	    qp=q;
	    if ((circ = zp[q]-zp[q-1])<0) circ=npts+circ;
	    for (area=0, k=zp[q-1]; k<zp[q-1]+circ; ++k)
		area += xp[k%npts]*yp[(k+1)%npts]-xp[(k+1)%npts]*yp[k%npts];
	    area += xp[zp[q]]*yp[zp[q-1]] - xp[zp[q-1]]*yp[zp[q]];
	    nopn = MAXV(ABSV(xp[zp[q]]-xp[zp[q-1]]),
			ABSV(yp[zp[q]]-yp[zp[q-1]]));
	    circ += nopn;
	    area = ABSV(area/2)-circ/2;
/* if area > min save info on defficiency */
	    if (area > HOLE_MIN) {
		cdeff[ndeff].area = area;
		cdeff[ndeff].len = nopn;
		cdeff[ndeff].circ = circ;
		cdeff[ndeff].n1 = zp[q-1];
		cdeff[ndeff].x1 = xp[zp[q-1]];
		cdeff[ndeff].y1 = yp[zp[q-1]];
		cdeff[ndeff].n2 = zp[q];
		cdeff[ndeff].x2 = xp[zp[q]];
		cdeff[ndeff].y2 = yp[zp[q]];
	        ang = atan2((double)(yp[zp[q]]-yp[zp[q-1]]),
			    (double)(xp[zp[q]]-xp[zp[q-1]]));
		cdeff[ndeff].ang = ang;
		printf("deff %d, area %d, nopn %d, circ %d ang %f\n",
			ndeff,area,nopn,circ,ang*180.0/M_PI);
		printf("     edge pt 1 n:%d x:%d y:%d 2 n:%d x:%d y:%d\n",
		zp[q-1],xp[zp[q-1]],yp[zp[q-1]],zp[q],xp[zp[q]],yp[zp[q]]);
		++ndeff;
		if (ndeff>=MAXDEFF) printf("Error - too many con deffs\n");
	    }
	}
	else if (zp[q]==vp[i]) ++q;
    }
    return(ndeff);
}


int 
right_shft(int p,int a,int b,int *xp,int *yp)
{
    int dx,dy;

    if (ABSV(xp[b]-xp[a]) > ABSV(yp[b]-yp[a]))
    {
	dx=0; dy= -1;
	if (xp[b]>xp[a]) dy=1;
    }
    else
    { 
	dx=1; dy=0;
	if (yp[b]>yp[a]) dx= -1;
    }
    if ((yp[p]-yp[a]-dy)*(xp[b]-xp[a])>(xp[p]-xp[a]-dx)*(yp[b]-yp[a]))
	return(1);
    else return(0);
}


