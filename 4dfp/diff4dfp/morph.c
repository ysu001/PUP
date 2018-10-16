#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <knee.h>

/***** utilities *****/
/*
** name  - get_image
** use   - get image from file, raw uchars
** input - image: pointer to image array
**         filename: name of image file
**         size: number of elements
** output- error message only
** effect- fills array image with data
** note  - none
*/

void 
get_image(unsigned char *image,char *filename,int size)
{
    FILE *fp;

    if ((fp = fopen(filename, "r")) == NULL) {
        printf("\r\n Error opening file %s.\r\n", filename);
        return;
    }
    printf("Reading %s size %d\n",filename,size);
    fread(image, size, 1, fp);
    fclose(fp);
}

/*
** name  - put_image
** use   - put image to file, raw uchars, no headers
** input - image: pointer to image array
**         filename: name of image file
**         size: number of elements
** output- error msg only
** effect- creates file with array image data
** note  - none
*/
void 
put_image(unsigned char *image,char *filename,int size)
{
    FILE *fp;

    if ((fp = fopen(filename, "w")) == NULL) {
        printf("\r\n Error opening file %s.\r\n", filename);
        return;
    }
    printf("Writing %s size %d\n",filename,size);
    fwrite(image, size, 1, fp);
    fclose(fp);
}

/*
** name  - dicom_hdr
** use   - extract header information from dicom file
**	   assume the header is 496 bytes long
*/
int
dicom_hdr(short int *src,int nhdr,int *ncol,int *nrow,int *npix)
{
    int i;
/*
    for (i=0; i<nhdr/2; ++i) printf("%d %d \n",i,src[i]);
*/
    *ncol = src[223 + (nhdr-496)/2];
    *nrow = src[218 + (nhdr-496)/2];
    *npix = (*ncol)*(*nrow);
    return(2*(*npix)+nhdr);
}

/*
** name  - inv_chk_img
** use   - invert and check image
** input - src: pointer to source array
**	   npix: number of pix
** output- returns changes in src
*/
void
inv_chk_img(short int *src, int npix)
{
    int i;
    short int tmp;

    for (i=0; i<npix; ++i) {
	tmp = NGRAY - src[i];
	if (tmp>=NGRAY) src[i]=NGRAY-1;
	else if (tmp<0) src[i]=0;
	else src[i]=tmp;
    }

/* old conversion 
    for (i=nhdr; i<size; i+=2) {
	dst[(i-hdr)/2]= 4096 - (256*src[i]+src[i+1]);
    }
*/
}

/* 
** name - compress4
** use - reduce a file by a factor of 4
**       divides the ncol and nrow by two and saves these values.
** note - overwrites the input file
*/
void 
compress4(short int *src,int *ncol,int *nrow)
{
    int i,j,ncolo;

    ncolo = *ncol;
    *ncol /= 2;
    *nrow /= 2;
    for (i=0; i< *nrow; ++i) {
    for (j=0; j< *ncol; ++j) {
	*(src+i*(*ncol)+j) = (
	*(src+(2*i+1)*ncolo+2*j) + *(src+(2*i+1)*ncolo+2*j+1) +
	*(src+(2*i)*ncolo+2*j) + *(src+(2*i)*ncolo+2*j+1) )/4;
    }
    }
}

/*
** name  - vax_to_uchar
** use   - convert from vax format to raw uchars, no header
** input - src: pointer to source array
**         dst: pointer to destination array
** output- none
** effect- creates new array in raw, uchar format 256x256
** note  - this function deals only with 256x256 images. 
**         vax format has a 4K header and each pixel is stored
**         in two bytes. the header is removed and if the high
**         order bit is set the uchar is set to 0xff. otherwise
**         the low order byte is used.
*/
void
vax_to_uchar(unsigned char *src, unsigned char *dst)
{
    int i;

    for (i=4096; i<135168; i+=2) {
        if (src[i+1])  dst[(i-4096)/2] = 0xff;
        else dst[(i-4096)/2] = src[i];
    }
}

/*
** int_to_char - convert from integer to string
** name  - int_to_char
** use   - convert from integer to string 
** input - num: integer to convert
**         num_text: char pointer to result
** output- none
** effect- fills char array with result
** note  - none
*/
void
int_to_char(int num, char *num_text)
{
    int i=0, j=0, dig;
    char rnum[20];

/* convert one digit at a time */
    if (num==0) { rnum[0] = '0'; i=1; }
    while (num != 0) {
        dig = num % 10;
        num /= 10;
        rnum[i] = (char) (dig + 48);
        ++i;
    }
    rnum[i] = 0x00;
/* invert the order of rnum for final result */
    while (i>0) {
        --i;
        num_text[j] = rnum[i];
        ++j;
    }
    num_text[j] = 0x00;
}


/*
** name - prt_hex_2D
** use  - make a hex 2D print outs
** input - blk: pointer to array 
**         x,y: center location
** output- none
** effect - print out in hex
** note - none
*/
void
prt_hex_2D(unsigned char *blk, int x, int y)
{
    int i,j,indx;

    printf("2D print centered on location %d %d\n",x,y);
    for (i=y-9; i<y+9; ++i) {
        for (j=x-9; j<x+9; ++j) {
            indx = i*DY + j;
            if (indx < 0 || indx >= DZ) printf("** ");
            else if (i==y && j==x) printf("%2x<",*(blk+indx));
            else printf("%2x ",*(blk+indx));
        }
        printf("\n");
    }
}


/***** point processes *****/

/*
** name  - segment_img
** use   - segment an image into forground and background
** input - npixel: size of image
**         src: pointer to source image
**         dst: pointer to destination image
**         lo: low limit of segment
**         hi: high limit of segment
** output- none
** effect- fills dst image with segmented image
** note  - done via look up table to save time
*/
void
segment_img(short int *src, short int *dst, int lo, int hi, int npix)
{
    short int lut[NGRAY];
    int i;

    for (i=0; i<NGRAY; ++i) {
        if (i < lo || i > hi) lut[i] = 0;
        else lut[i] = NGRAY-1;
    }
    for (i=0; i<npix; ++i) {
        *dst = lut[*src];
        dst++; src++;
    }
}

void
segment_img_bin(short int *src, short int *dst, int n, int npix)
{
    int i;

    /* dilate case */
    if (n < 9) {
	for (i=0; i<npix; i++) {
            if (*src < 9 && *src >= n) *dst = 1;
	    else *dst = (int) *src/9;
            dst++; src++;
	}
    }
    /* erode case */
    else {
	for (i=0; i<npix; i++) {
            if (*src >=9 && *src < n) *dst = 0;
	    else *dst = (int) *src/9;
            dst++; src++;
	}
    }
}

/*
** name - hist_img 
** use - histogram an image
** input - npixel: number of pixel in image
**         img: pointer to image
**         bins: pointer to bin array, size 256
** output - values in bin array
** effect - print out of max,avg,var
**          and of all bin numbers
*/
void
hist_img(short int *img, int npix)
{
    int i, imax, max=0,ufbin=0,ofbin=0,high=0;
    int bins[NGRAY];
    double sum=0.0, sum2=0.0, avg, var;

/* init and sum values */
    for (i=0; i<NGRAY; ++i) bins[i] = 0;
    for (i=0; i<npix; ++i) {
	if (*img < 0) ++ufbin;
	else if (*img >= NGRAY) ++ofbin;
        else bins[(int) *img]++;
        ++img;
    }

/* calculate stats */
    for (i=0; i<NGRAY; ++i) {
        sum += bins[i]*i;
        sum2 += bins[i]*i*i;
        if (bins[i] > max){ imax=i; max=bins[i]; }
        if (bins[i]!=0) high=i;
    }
    avg = sum/npix;
    var = sum2/npix - avg*avg; 
    printf("hist max bin[%d]=%d avg=%10.3f var=%10.3f under=%d over=%d\n",
	 imax, max, avg, var, ufbin, ofbin);

    draw_hist(bins,10000);
/* print values of bin 
    for (i=0; i<=high; ++i) {
	printf("b[%d]=%d ",i,bins[i]);
	if (i%4==3) printf("\n");
    }
*/
}

/*
** name - draw_hist 
** use - create a visual picture of histogram
** input - bin: pointer to filled bins arrray
** 	   max: max height in graph, if =0 auto calc
** output- into dst image
** effect - none
** note - adds ticks every 32 pixels
*/
void
draw_hist(int *bin, int max)
{
    int i,j,sum,amax,ncol=256,nrow=256,*sbin;
    short int *hist;

    sbin = (int *) malloc(ncol*sizeof(int));
    hist = (short int *) malloc(ncol*nrow*sizeof(short int));

/* rebin and find max value */
    amax=0;
    for (i=0; i<NGRAY; i+=(NGRAY/ncol)) {
	sum=0;
	for (j=0; j<(NGRAY/ncol); ++j) sum += *(bin+i+j);
	*(sbin+i*ncol/NGRAY)=sum;
	if (sum>amax) amax=sum;
    }
    if (max!=0) amax=max;

/* fill in the histogram */
    for (i=0; i<ncol; ++i) {
        for (j=0; j<nrow; ++j) {
            if (j <= (*(sbin+i)*(nrow-1)/amax)) *(hist+i+ncol*j) = NGRAY-1;
            else *(hist+i+ncol*j) = 0;
        }
    }

/* draw in ticks ,#ticks,length,bot,top
    for (i=1; i<8; ++i) {                   
        for (j=0; j<5; ++j) {               
            *(dst+i*32+256*(255-j)) = 0xff; 
            *(dst+i*32+256*j) = 0x00;       
        }
    }
*/
    put_image((unsigned char *)hist,"hist.img",2*ncol*nrow);
    free(hist); free(sbin);
}

/*
** name - profile_img
** use - get a row or column profile of image
** input - src: pointer to source imagee
**         x,y,z: coordinates for profile
**         bins: pointer to bins array, size 255
** output - in bins array
** effect - none
** note - a neg number in x,y, or z indicates the axis to scan
*/
void
profile_img(unsigned char *src, int x, int y, int z, int *bins)
{
    int i;

    if (x<0) {   
        for (i=0; i<NGRAY; ++i) bins[i] = *(src + IND(i,y,z));
    }
    if (y<0) {   
        for (i=0; i<NGRAY; ++i) bins[i] = *(src + IND(x,i,z));
    }
    if (z<0) {   
        for (i=0; i<NGRAY; ++i) bins[i] = *(src + IND(x,y,i));
    }
}


/***** area processes *****/

/*
** name - conv3x3
** use  - 3x3 convolution of image
** input - npixel: number of pixels in image
**         src: pointer  to source image
**         dst: pointer to destination image 
**         conv: 9 element convolution  matrix
**         norm: normalization factor
** output- none
** effect- fill dst matrix by convolution  result 
** note  - no convolution is performed on outer edge
**         all calculation is in integers
**         assume DY has line length
*/
void
conv3x3(short int *src,short int *dst,int *conv,int norm,
	int ncol, int npix)
{
    short int *i0,*i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;
    int i,res,c0,c1,c2,c3,c4,c5,c6,c7,c8;

/* reset image origin to avoid image frame and fill with zero */
    for (i=0; i<ncol+1; ++i) {
	*dst=0; ++dst;
    }
    src += ncol+1;

/* save computation time */
    i0 = src-ncol-1; i1 = src-ncol; i2 = src-ncol+1;
    i3 = src-1;      i4 = src;      i5 = src+1;
    i6 = src+ncol-1; i7 = src+ncol; i8 = src+ncol+1;

    c0 = conv[0]; c1 = conv[1]; c2 = conv[2];
    c3 = conv[3]; c4 = conv[4]; c5 = conv[5];
    c6 = conv[6]; c7 = conv[7]; c8 = conv[8];

/* modify npixel to avoid edges */
    npix -= 2*(ncol+1);
/* loop over pixels */
    for (i=0; i<npix; ++i) {
        res = ( *i0*c0 + *i1*c1 + *i2*c2 +
                *i3*c3 + *i4*c4 + *i5*c5 +
                *i6*c6 + *i7*c7 + *i8*c8)/norm;
        res = ABSV((int) res);
        if (res > NGRAY-1) *dst = NGRAY-1;
        else *dst = res;
        dst++; i0++; i1++; i2++; i3++; i4++; i5++; i6++; i7++; i8++;
    }
/* fill edge with zero */
    for (i=0; i<ncol+1; ++i) {
	*dst=0; ++dst;
    }
}



/*
** name - conv3x3x3
** use  - 3x3x3 convolution of image
** input - npixel: number of pixels in image
**         src: pointer  to source image
**         dst: pointer to destination image 
**         conv: 27 element convolution  matrix
**         norm: normalization factor
** output- none
** effect- fill dst matrix by convolution  result 
** note  - no convolution is performed on outer edge
**         all calculation is in integers
**         assume DY has line length
*/
void
conv3x3x3(short int *src,short int *dst,int *conv,int norm,
	int ncol, int npix, int nvol)
{
    short int *i00,*i01,*i02,*i03,*i04,*i05,*i06,*i07,*i08;
    short int *i09,*i10,*i11,*i12,*i13,*i14,*i15,*i16,*i17;
    short int *i18,*i19,*i20,*i21,*i22,*i23,*i24,*i25,*i26;
    int i,res;

/* reset image origin to avoid image frame and fill with zero */
    for (i=0; i<npix+ncol+1; ++i) {
	*dst=0; ++dst;
    }
    src += npix+ncol+1;

/* save computation time */
    i00 = src-npix-ncol-1; i01 = src-npix-ncol; i02 = src-npix-ncol+1;
    i03 = src-npix-1;      i04 = src-npix;      i05 = src-npix+1;
    i06 = src-npix+ncol-1; i07 = src-npix+ncol; i08 = src-npix+ncol+1;

    i09 = src-ncol-1; i10 = src-ncol; i11 = src-ncol+1;
    i12 = src-1;      i13 = src;      i14 = src+1;
    i15 = src+ncol-1; i16 = src+ncol; i17 = src+ncol+1;

    i18 = src+npix-ncol-1; i19 = src+npix-ncol; i20 = src+npix-ncol+1;
    i21 = src+npix-1;      i22 = src+npix;      i23 = src+npix+1;
    i24 = src+npix+ncol-1; i25 = src+npix+ncol; i26 = src+npix+ncol+1;

/* modify npixel to avoid edges */
    nvol -= 2*(npix+ncol+1);
/* loop over pixels */
    for (i=0; i<nvol; ++i) {

        res = ( *i00* conv[0] + *i01* conv[1] + *i02* conv[2] +
                *i03* conv[3] + *i04* conv[4] + *i05* conv[5] +
                *i06* conv[6] + *i07* conv[7] + *i08* conv[8] +
        	*i09* conv[9] + *i10*conv[10] + *i11*conv[11] +
                *i12*conv[12] + *i13*conv[13] + *i14*conv[14] +
                *i15*conv[15] + *i16*conv[16] + *i17*conv[17] +
        	*i18*conv[18] + *i19*conv[19] + *i20*conv[20] +
                *i21*conv[21] + *i22*conv[22] + *i23*conv[23] +
                *i24*conv[24] + *i25*conv[25] + *i26*conv[26]  )/norm;

	/*
        res = ABSV((int) res);
	*/
        if (res > NGRAY-1) *dst = NGRAY-1;
        else *dst = res;
        dst++; 
	i00++; i01++; i02++; i03++; i04++; i05++; i06++; i07++; i08++;
	i09++; i10++; i11++; i12++; i13++; i14++; i15++; i16++; i17++;
	i18++; i19++; i20++; i21++; i22++; i23++; i24++; i25++; i26++;
    }
/* fill edge with zero */
    for (i=0; i<npix+ncol+1; ++i) {
	*dst=0; ++dst;
    }
}


void
sobelx(short int *src,short int *dst,int ncol,int npix)
{
    short int *i0,*i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;
    int i,res;

/* reset image origin to avoid image frame and fill with zero */
    for (i=0; i<ncol+1; ++i) {
	*dst=0; ++dst;
    }
    src += ncol+1;

/* save computation time */
    i0 = src-ncol-1; i1 = src-ncol; i2 = src-ncol+1;
    i3 = src-1;      i4 = src;      i5 = src+1;
    i6 = src+ncol-1; i7 = src+ncol; i8 = src+ncol+1;

/* modify npixel to avoid edges */
    npix -= 2*(ncol+1);
/* loop over pixels */
    for (i=0; i<npix; ++i) {
        res = -(*i0) + *i2 - *i3*2 + *i5*2 - *i6  + *i8;
        res = ABSV((int) res);
        if (res > NGRAY-1) *dst = NGRAY-1;
        else *dst = res;
        dst++; i0++; i2++; i3++; i5++; i6++; i8++;
    }
/* fill edge with zero */
    for (i=0; i<ncol+1; ++i) {
	*dst=0; ++dst;
    }
}

void
sobely(short int *src,short int *dst,int ncol,int npix)
{
    short int *i0,*i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;
    int i,res;

/* reset image origin to avoid image frame and fill with zero */
    for (i=0; i<ncol+1; ++i) {
	*dst=0; ++dst;
    }
    src += ncol+1;

/* save computation time */
    i0 = src-ncol-1; i1 = src-ncol; i2 = src-ncol+1;
    i3 = src-1;      i4 = src;      i5 = src+1;
    i6 = src+ncol-1; i7 = src+ncol; i8 = src+ncol+1;

/* modify npixel to avoid edges */
    npix -= 2*(ncol+1);
/* loop over pixels */
    for (i=0; i<npix; ++i) {
        res = -(*i0) - *i1*2 - *i2 + *i6 + *i7*2 + *i8;
        res = ABSV((int) res);
        if (res > NGRAY-1) *dst = NGRAY-1;
        else *dst = res;
        dst++; i0++; i1++; i2++; i6++; i7++; i8++;
    }
/* fill edge with zero */
    for (i=0; i<ncol+1; ++i) {
	*dst=0; ++dst;
    }
}

void
conv_de(short int *src,short int *dst,int ncol,int npix)
{
    short int *i0,*i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;
    int i,res;

/* reset image origin to avoid image frame and fill with zero */
    for (i=0; i<ncol+1; ++i) {
	*dst=0; ++dst;
    }
    src += ncol+1;

/* save computation time */
    i0 = src-ncol-1; i1 = src-ncol; i2 = src-ncol+1;
    i3 = src-1;      i4 = src;      i5 = src+1;
    i6 = src+ncol-1; i7 = src+ncol; i8 = src+ncol+1;

/* modify npixel to avoid edges */
    npix -= 2*(ncol+1);
/* loop over pixels */
    for (i=0; i<npix; ++i) {
        res = ( *i0 + *i1 + *i2 +
                *i3 + *i4*9 + *i5 +
                *i6 + *i7 + *i8)/NGRAY-1;
        res = ABSV((int) res);
        if (res > NGRAY-1) *dst = NGRAY-1;
        else *dst = res;
        dst++; i0++; i1++; i2++; i3++; i4++; i5++; i6++; i7++; i8++;
    }
/* fill edge with zero */
    for (i=0; i<ncol+1; ++i) {
	*dst=0; ++dst;
    }
}

/* assumes a binary input 0 or 1 */
void
conv_de_bin(short int *src,short int *dst,int ncol,int npix)
{
    short int *i0,*i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;
    int i,res;

/* reset image origin to avoid image frame and fill with zero */
    for (i=0; i<ncol+1; ++i) {
	*dst=0; ++dst;
    }
    src += ncol+1;

/* save computation time */
    i0 = src-ncol-1; i1 = src-ncol; i2 = src-ncol+1;
    i3 = src-1;      i4 = src;      i5 = src+1;
    i6 = src+ncol-1; i7 = src+ncol; i8 = src+ncol+1;

/* modify npixel to avoid edges */
    npix -= 2*(ncol+1);
/* loop over pixels */
    for (i=0; i<npix; ++i) {
        res = ( *i0 + *i1 + *i2 +
                *i3 + *i4*9 + *i5 +
                *i6 + *i7 + *i8);
        *dst = res;
        dst++; i0++; i1++; i2++; i3++; i4++; i5++; i6++; i7++; i8++;
    }
/* fill edge with zero */
    for (i=0; i<ncol+1; ++i) {
	*dst=0; ++dst;
    }
}

/*
** name  - dilate_erode
** use   - performs morphological primitives on image such as
**         dilate, fill, erode, delete
** input - npixel: size of image
**         src: pointer to source image
**         dst: pointer to destination image
**         n: number of 8-neighbors. use the following code
**         dilate(n)  - 0 to 1 if #n=1 is >=n  (also called grow)
**         fill (8-n) - 0 to 1 if #n=0 is <=n
**         erode(17-n) - 1 to 0 if #n=0 is >=n
**         delete(9+n) - 1 to 0 if #n=1 is <=n
** output- none
** effect- fills dst with operation on src
** note  - n ranges from 1 (coarse) to 7(fine).
*/
void
dilate_erode(short int *src,short int *dst,int n,int ncol,int npix)
{
    void conv_de(short int *src,short int *dst,int ncol,int npix);
    int i;
    short int *tmp;

    tmp = (short int *) malloc(npix*sizeof(short int));
    conv_de(src, tmp,ncol,npix);
    segment_img(tmp, dst, n, NGRAY-1, npix);
    free(tmp);
}

/* expect a binary 0 or 1 input */
/* n<9 dilate; n>=9 dilate */
void
dilate_erode_bin(short int *src,short int *dst,int n,int ncol,int npix)
{
    void conv_de_bin(short int *src,short int *dst,int ncol,int npix);
    int i;
    short int *tmp;

    tmp = (short int *) malloc(npix*sizeof(short int));
    conv_de_bin(src, tmp,ncol,npix);
    segment_img_bin(tmp, dst, n, npix);
    free(tmp);
}

/*
** name - open_img
** use  - erode then dilate, leaves less area, smooths convexities
** note - see comments in dilate_erode for details
*/
void
open_img(short int *src, short int *dst, int n, int ncol, int npix)
{
    short int *tmp;

    tmp = (short int*) malloc(npix*sizeof(short int));
    dilate_erode(src, tmp, 17-n, ncol, npix);
    dilate_erode(tmp, dst, n, ncol, npix);
    free(tmp);
}

/*
** name - close_img
** use  - dilate then erode, leaves more area, fills concavities
** note - see comments in dilate_erode for details
*/
void
close_img(short int *src, short int *dst, int n, int ncol, int npix)
{
    short int *tmp;

    tmp = (short int*) malloc(npix*sizeof(short int));
    dilate_erode(src, tmp, n, ncol, npix);
    dilate_erode(tmp, dst, 17-n, ncol, npix);
    free(tmp);
}

/*
** name - dilate_erode_3D
** use  - 3D version of dilate_erode
** note - see commentss in dilate_erode for details
*/
void
dilate_erode_3D(int npixel, unsigned char *src, int n, unsigned char *dst)
{
    int i, prim[27];
    unsigned char *tmp;

    for (i=0; i<27; ++i) prim[i]=1; prim[13]=9;
    tmp = (unsigned char *) malloc(npixel*sizeof(unsigned char));
    for (i=0; i<npixel; ++i) tmp[i] = 0x00;
/*
    segment_img(npixel, tmp, n, 255, dst);
*/
    free(tmp);
}


/*
** name - avgpix
** use  - get the average of pixels in image under a mask
** input - src: pointer to image to be averaged
**         msk: pointer to mask image
**         dst: pointer to resulting average image
**         min,max: range over which to average
**         span: area to average over +/- span
** output - result in dst
** effect - none
** note - none
*/
void
avgpix(unsigned char *src, unsigned char *msk, unsigned char *dst,
	int xmin,int xmax,int xspan,int ymin,int ymax,int yspan,int zmin,int zmax,int zspan)
{
    int csum, cnum, psum[65], pnum[65];
    int sum, num, i,j,k, x,y,z;
    int xsmin, xsmax, ysmin, ysmax, zsmin, zsmax;

    for (z=zmin; z<zmax; ++z) 
      for (y=ymin; y<ymax; ++y) 
        for (x=xmin; x<xmax; ++x) {
          zsmin=z-zspan; zsmax=z+zspan;
          ysmin=y-yspan; ysmax=y+yspan;
          xsmin=x-xspan; xsmax=x+xspan;
          sum = 0; num = 0;
          for (k=zsmin; k<=zsmax; ++k) 
            for (j=ysmin; j<=ysmax; ++j) 
              for (i=xsmin; i<=xsmax; ++i) {
                if (*(msk+k*DZ+j*DY+i)) {sum += *(src+k*DZ+j*DY+i); ++num;}
          }
          if (num==0) *(dst+z*DZ+y*DY+x) = 0x00;
          else *(dst+z*DZ+y*DY+x) = (unsigned char) (sum/num);
    }
}

/***** frame processes *****/

/*
** mag_img - get magnitude of combining two images
** name  - mag_img
** use   - get magnitude of combining two images as vectors
** input - npixel: size of image
**         srcx: pointer to first source image
**         srcy: pointer to second source image
**         dst: pointer to destination image
** output- none
** effect- fills dst = sqrt( srcx^2 + srcy^2 )
** note  - none
*/
void
mag_img(short int *srcx,short int *srcy,short int *dst, int npix)
{
    int i, res;
    for (i=0; i<npix; ++i) {
/*
        res = (int) sqrt((double) (*srcx)* (double) (*srcx)+
			 (double) (*srcy)* (double) (*srcy));
*/
	res = *srcx + *srcy;
        if (res > NGRAY-1) *dst = NGRAY-1;
        else *dst = res;
        dst++; srcx++; srcy++;
    }
}

/*
** name  - linear_img
** use   - combine two images in linear way
** input - npixel: size of image
**         src0: pointer to first source image
**         src1: pointer to second source image
**         dst: pointer to destination image
**         c0,c1,c2: constants for combination
**         ov_flag: signal how to deal with underflow
** output- none
** effect- fills dst = c0*src0 + c1*src1 +c2
** note  - none
*/
void
linear_img(short int *src0, short int *src1, short int *dst, 
	   int c0, int c1, int c2, int ov_flag, int npix)
{
    int i, res;

    for (i=0; i<npix; ++i) {
        res = (*src0)*c0 + (*src1)*c1 + c2;
        if (res < 0) { 
            if (ov_flag) res = 0;
            else res = ABSV(res);
        }
        if (res > NGRAY-1) res = NGRAY-1;
        *dst = res;
        dst++; src0++; src1++;
    }
}

/*
** name  - scale_img
** use   - scale the pixels in an image
** input - npixel: size of image
**         src: pointer to source image
**         dst: pointer to destination image
**         c0,c1: constants for scaling
**         ov_flag: signal how to deal with underflow
** output- none
** effect- fills dst = c0*src0 + c1
** note  - none
*/
void
scale_img(short int *src, short int *dst,
	  int c0, int c1, int ov_flag, int npix)
{
    short int lut[NGRAY];
    int i, res;

    for (i=0; i<NGRAY; ++i) {
        res = i*c0 + c1;
        if (res < 0) { 
            if (ov_flag) res = 0;
            else res = ABSV(res);
        }
        if (res > NGRAY-1) res = NGRAY-1;
        lut[i] = res;
    }
    for (i=0; i<npix; ++i) {
        *dst = lut[*src];
        dst++; src++;
    }
}

/*
** name  - grad_img
** use   - perform a gradient of image
** input - npixel: size of image
**	   scale:scale the final answer
**         src: pointer to source image
**         dst: pointer to destination image
** output- none
** effect- fills dst 
** note  - none
*/
void
grad_img(short int *src,short int *dst,int scale,int ncol,int npix)
{
    short int *imgx,*imgy,*imgm;

    imgx = (short int*) malloc(npix*sizeof(short int));
    imgy = (short int*) malloc(npix*sizeof(short int));
    imgm = (short int*) malloc(npix*sizeof(short int));

    sobelx(src,imgx,ncol,npix);
    sobely(src,imgy,ncol,npix);
    mag_img(imgx,imgy,imgm,npix);
    scale_img(imgm,dst,scale,1,1,npix);
    free(imgx); free(imgy); free(imgm);
}


/*
** name  - and_img
** use   - and two images pixel by pixel
** input - npixel: size of image
**         src0: pointer to first source image
**         src1: pointer to second source image
**         dst: pointer to destination image
** output- none
** effect- fills dst = src0 and src1
** note  - none
*/
void
and_img(short int *src0,short int *src1,short int *dst,int npix)
{
    int i;
    for (i=0; i<npix; ++i) {
        *dst = *src0 & *src1;
        src0++; src1++; dst++;
    }
}

/*
** name  - xor_img
** use   - xor two images pixel by pixel
** input - npixel: size of image
**         src0: pointer to first source image
**         src1: pointer to second source image
**         dst: pointer to destination image
** output- none
** effect- fills dst = src0 xor src1
** note  - none
*/
void
xor_img(int npixel, unsigned char *src0, unsigned char *src1, unsigned char *dst)
{
    int i;
    for (i=0; i<npixel; ++i) {
        *dst = *src0 ^ *src1;
        src0++; src1++; dst++;
    }
}


/*
** name  - or_img
** use   - or two images pixel by pixel
** input - npixel: size of image
**         src0: pointer to first source image
**         src1: pointer to second source image
**         dst: pointer to destination image
** output- none
** effect- fills dst = src0 or src1
** note  - none
*/
void
or_img(int npixel, unsigned char *src0, unsigned char *src1, unsigned char *dst)
{
    int i;
    for (i=0; i<npixel; ++i) {
        *dst = *src0 | *src1;
        src0++; src1++; dst++;
    }
}

/*
** name  - not_img
** use   - not an image pixel by pixel
** input - npixel: size of image
**         src: pointer to source image
**         dst: pointer to destination image
** output- none
** effect- fills dst = not src
** note  - none
*/
void
not_img(short int *src, short int *dst, int npix)
{
    short int lut[NGRAY];
    int i;

    lut[0] = NGRAY-1;
    for (i=1; i<NGRAY; ++i) lut[i] = 0;
    for (i=0; i<npix; ++i) {
        *dst = lut[*src];
        dst++; src++;
    }
}


/***** geometrical processes *****/

/* 
** name  - zoom_img
** use   - zoom magnify an image by some power of 2
** input - fac2: zoom factor, use power of 2
**         src: pointer to source image
**         dst: pointer to destination image
**         x,y: center point of zoom in src
** output- none
** effect- fills dst by zoomed version of src
** note  - assumes 256x256 images
*/
void
zoom_img(int fac2, unsigned char *src, int x, int y, unsigned char *dst)
{
    int i,j,dx,dy;

    for (i=0; i<DY; ++i) {
        for (j=0; j<DY; ++j) {
            dx = (j-128)/fac2;
            dy = (i-128)/fac2;
            if (y+dy<0 || y+dy>255 || x+dx<0 || x+dx>255) 
                *(dst+i*DY+j) = 0;
            else *(dst+i*DY+j) = *(src+(y+dy)*DY+x+dx);
        }
    }
}
