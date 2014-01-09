#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>
#include <librms.h>
#include <nrutil.h>
#include <RSF.h>
	
IMAGE_4dfp *read_image_4dfp(char *filespc, IMAGE_4dfp *image)
{
	int n, m1, m2, isbig;
	char imgroot[MAXL], imgfile[MAXL];
	FILE *imgfp=NULL;
	
	if (!image) {
		image=(IMAGE_4dfp *)malloc((size_t)sizeof(IMAGE_4dfp));
		if (!image) errm("read_image_4dfp");
		getroot(filespc, imgroot);
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		Getifh(imgfile, &(image->ifh));
		m2=image->ifh.number_of_bytes_per_pixel*image->ifh.matrix_size[0]*image->ifh.matrix_size[1]*image->ifh.matrix_size[2]*image->ifh.matrix_size[3];		
	} else {
		m1=image->ifh.number_of_bytes_per_pixel*image->ifh.matrix_size[0]*image->ifh.matrix_size[1]*image->ifh.matrix_size[2]*image->ifh.matrix_size[3];
		getroot(filespc, imgroot);
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		Getifh(imgfile, &(image->ifh));
		m2=image->ifh.number_of_bytes_per_pixel*image->ifh.matrix_size[0]*image->ifh.matrix_size[1]*image->ifh.matrix_size[2]*image->ifh.matrix_size[3];		
		if (m1!=m2) { /*free old memory if the image size does not match*/
			free((void *) image->image);
		}	 
	}
	isbig = strcmp (image->ifh.imagedata_byte_order, "littleendian");
	if (!(imgfp = fopen (imgfile, "rb"))) errr ("read_image_4dfp", imgfile);
	m2=m2/image->ifh.number_of_bytes_per_pixel;
	image->image=(float *)malloc((size_t)(m2*sizeof(float))); /*Allocate memory for reading the image*/
	if (!image->image) errm("read_image_4dfp");
	if (eread(image->image, m2, isbig, imgfp)) errr("read_image_4dfp", imgfile);
	fclose(imgfp);
	return image;	
}

int write_image_4dfp(char *filespc, IMAGE_4dfp *image)
{
	int isbig, vdim, status = 0;
	char imgroot[MAXL], imgfile[MAXL], command[MAXL];
	FILE *imgfp;
	char control = '\0';
	
	isbig = strcmp (image->ifh.imagedata_byte_order, "littleendian");
	if (!control) control = (isbig) ? 'b' : 'l';
	vdim=image->ifh.matrix_size[0]*image->ifh.matrix_size[1]*image->ifh.matrix_size[2]*image->ifh.matrix_size[3];
	getroot(filespc, imgroot);
	sprintf(imgfile, "%s.4dfp.img", imgroot);
	if (!(imgfp = fopen (imgfile, "wb"))) errw("write_image_4dfp", imgfile);
	if (ewrite(image->image, vdim, control, imgfp)) errw("write_image_4dfp", imgfile);
	
	/* ifh hdr rec */
	if (fclose(imgfp)) errw("write_image_4dfp", imgfile);
	if (Writeifh("write_image_4dfp", imgfile, &image->ifh, control)) errw("write_image_4dfp", imgroot);
	sprintf(command, "ifh2hdr %s", imgroot);
	status |= system (command);
	return status;
}
	
float getroimean_4dfp(IMAGE_4dfp *in_image, IMAGE_4dfp *Mask, int val)
{
	float *fptr1, *fptr2, roimean, lo, hi;
	int i, n, nv;
	
	/*need to implement a input checking mechanism*/
	
	roimean=0;
	n=0;
	fptr1 = in_image->image;
	fptr2 = Mask->image;
	lo=(float)val - 0.5;
	hi=(float)val + 0.5;
	nv = in_image->ifh.matrix_size[0]*in_image->ifh.matrix_size[1]*in_image->ifh.matrix_size[2]*in_image->ifh.matrix_size[3];
	for (i=0; i<nv; i++)
	{
		if (*fptr2<hi && *fptr2>lo) {
			n++;
			roimean+=*fptr1;
		}
		fptr1++;
		fptr2++;
	}
	if (n>0) {
		roimean=roimean/(float)n;
	} else {
		printf("Warning: ROI has no voxels\n");
	}
	return roimean;
}


IMAGE_4dfp *copy_image_4dfp(IMAGE_4dfp *in_image, IMAGE_4dfp *out_image)
{
	int i, m1, m2;
	float *fptr1, *fptr2;
	
	m2=in_image->ifh.number_of_bytes_per_pixel*in_image->ifh.matrix_size[0]*in_image->ifh.matrix_size[1]*in_image->ifh.matrix_size[2]*in_image->ifh.matrix_size[3];
	if (!out_image) {
		out_image=(IMAGE_4dfp *)malloc((size_t)sizeof(IMAGE_4dfp));
		if (!out_image) errm("copy_image_4dfp");
		out_image->image=(float *)malloc((size_t)m2);
		if (!out_image->image) errm("copy_image_4dfp");	
	} else {
		m1=out_image->ifh.number_of_bytes_per_pixel*out_image->ifh.matrix_size[0]*out_image->ifh.matrix_size[1]*out_image->ifh.matrix_size[2]*out_image->ifh.matrix_size[3];
		if (m1!=m2) {
			free((void *) out_image->image);
			out_image->image=(float *)malloc((size_t)m2);
			if (!out_image->image) errm("copy_image_4dfp");
		}
		 
	}
	strcpy(out_image->ifh.interfile, in_image->ifh.interfile);
	strcpy(out_image->ifh.version_of_keys, in_image->ifh.version_of_keys);
	strcpy(out_image->ifh.conversion_program, in_image->ifh.conversion_program);
	strcpy(out_image->ifh.name_of_data_file, in_image->ifh.name_of_data_file);
	/* needs to modify out_image->name_of_data_file to avoid duplication of file names*/
	strcpy(out_image->ifh.number_format, in_image->ifh.number_format);
	strcpy(out_image->ifh.imagedata_byte_order, in_image->ifh.imagedata_byte_order);
	out_image->ifh.number_of_bytes_per_pixel=in_image->ifh.number_of_bytes_per_pixel;
	out_image->ifh.number_of_dimensions=in_image->ifh.number_of_dimensions;
	for (i=0; i<3; i++)
	{
		out_image->ifh.matrix_size[i]=in_image->ifh.matrix_size[i];
		out_image->ifh.scaling_factor[i]=in_image->ifh.scaling_factor[i];
		out_image->ifh.mmppix[i]=in_image->ifh.mmppix[i];
		out_image->ifh.center[i]=in_image->ifh.center[i];
	}
	out_image->ifh.matrix_size[3]=in_image->ifh.matrix_size[3];
	out_image->ifh.scaling_factor[3]=in_image->ifh.scaling_factor[3];
	out_image->ifh.orientation=in_image->ifh.orientation;
	m2=m2/in_image->ifh.number_of_bytes_per_pixel;
	fptr1=in_image->image;
	fptr2=out_image->image;
	for (i=0; i<m2; i++)
	{
		*fptr2=*fptr1;
		fptr1++;
		fptr2++;
	}
	return out_image;
}

IMAGE_4dfp *extract_roi_4dfp(IMAGE_4dfp *image, int MaskVal, IMAGE_4dfp *roimask)
{
	int i, m1, m2;
	float *fptr1, *fptr2, lo, hi;
	
	m2=image->ifh.number_of_bytes_per_pixel*image->ifh.matrix_size[0]*image->ifh.matrix_size[1]*image->ifh.matrix_size[2]*image->ifh.matrix_size[3];
	if (!roimask) {
		roimask=(IMAGE_4dfp *)malloc((size_t)sizeof(IMAGE_4dfp));
		if (!roimask) errm("extract_roi_4dfp");
		roimask->image=(float *)malloc((size_t)m2);
		if (!roimask->image) errm("extract_roi_4dfp");	
	} else {
		m1=roimask->ifh.number_of_bytes_per_pixel*roimask->ifh.matrix_size[0]*roimask->ifh.matrix_size[1]*roimask->ifh.matrix_size[2]*roimask->ifh.matrix_size[3];
		if (m1!=m2) {
			free((void *) roimask->image);
			roimask->image=(float *)malloc((size_t)m2);
			if (!roimask->image) errm("extract_roi_4dfp");
		}
		 
	}
	strcpy(roimask->ifh.interfile, image->ifh.interfile);
	strcpy(roimask->ifh.version_of_keys, image->ifh.version_of_keys);
	strcpy(roimask->ifh.conversion_program, image->ifh.conversion_program);
	strcpy(roimask->ifh.name_of_data_file, image->ifh.name_of_data_file);
	/* needs to modify roimask->name_of_data_file to avoid duplication of file names*/
	strcpy(roimask->ifh.number_format, image->ifh.number_format);
	strcpy(roimask->ifh.imagedata_byte_order, image->ifh.imagedata_byte_order);
	roimask->ifh.number_of_bytes_per_pixel=image->ifh.number_of_bytes_per_pixel;
	roimask->ifh.number_of_dimensions=image->ifh.number_of_dimensions;
	for (i=0; i<3; i++)
	{
		roimask->ifh.matrix_size[i]=image->ifh.matrix_size[i];
		roimask->ifh.scaling_factor[i]=image->ifh.scaling_factor[i];
		roimask->ifh.mmppix[i]=image->ifh.mmppix[i];
		roimask->ifh.center[i]=image->ifh.center[i];
	}
	roimask->ifh.matrix_size[3]=image->ifh.matrix_size[3];
	roimask->ifh.scaling_factor[3]=image->ifh.scaling_factor[3];
	roimask->ifh.orientation=image->ifh.orientation;
	fptr1=image->image;
	fptr2=roimask->image;
	m2=m2/image->ifh.number_of_bytes_per_pixel;
	lo = (float)MaskVal - 0.5;
	hi = (float)MaskVal + 0.5;
	for (i=0; i<m2; i++)
	{
		if (*fptr1<hi && *fptr1>lo) {
			*fptr2=1;
		} else {
			*fptr2=0;
		}
		fptr1++;
		fptr2++;
	}
	return roimask;
}

void getroismean_4dfp(IMAGE_4dfp *in_image, IMAGE_4dfp *Mask, int n, float *a)
{
	float *fptr1, *fptr2, *p;
	int i, nv, idx;
	
	p=vector(1,n);
	for (i=1; i<=n; i++) {a[i]=0.0; p[i]=0.;}
	nv = in_image->ifh.matrix_size[0]*in_image->ifh.matrix_size[1]*in_image->ifh.matrix_size[2]*in_image->ifh.matrix_size[3];
	fptr1=in_image->image;
	fptr2=Mask->image;
	for (i=0; i<nv; i++) {
		idx=(int)rint((double)(*fptr2))+1;
		a[idx]+=*fptr1;
		p[idx]+=1.;
		if (idx>n) printf("something is wrong.\n");
		fptr1++;
		fptr2++;
	}
	for (i=1; i<=n; i++) {a[i]=a[i]/p[i];}
}

void getroismean2_4dfp(IMAGE_4dfp *in_image, IMAGE_4dfp *Mask, IMAGE_4dfp *PETMask, int n, float *a, float *p)
{
	float *fptr1, *fptr2, *fptr3, m;
	int i, nv, idx, nm;
	
	
	for (i=1; i<=n; i++) {a[i]=0.0; p[i]=0.;}
	nv = in_image->ifh.matrix_size[0]*in_image->ifh.matrix_size[1]*in_image->ifh.matrix_size[2]*in_image->ifh.matrix_size[3];
	fptr1=in_image->image;
	fptr2=Mask->image;
	fptr3=PETMask->image;
	nm=0;
	m=0;
	for (i=0; i<nv; i++) {
		if (*fptr3>0){
			idx=(int)rint((double)(*fptr2))+1; /*belong to ROI # idx*/
			a[idx]+=*fptr1; /*sum for ROI #idx + voxel value*/
			p[idx]+=1.; /* voxel number for ROI #idx increase by 1 */
			if (idx>n) printf("something is wrong.\n");
			if (*fptr2 >0 ){
				nm++;
				m+=*fptr1;
			}
		}
		fptr1++;
		fptr2++;
		fptr3++;
	}
	m=m/(float)nm*.8;
	for (i=1; i<=n; i++) {
		if (p[i]>0) {
			a[i]=a[i]/p[i];
			/*printf("%f\t%f\n",p[i], a[i]);*/
		} else { /* to protect cases that ROI #i is completely out of PET field of view */
			a[i]=m; /* in such case, the ROI mean is replaced with mean value of all non-unknown regions multiplied by 0.8 */
		}
	}
}

void getroismean3_4dfp(IMAGE_4dfp *in_image, IMAGE_4dfp *Mask, IMAGE_4dfp *PETMask, int n, float *a, float *p)
{
/* get roimean excluding non-valid voxels */
	float *fptr1, *fptr2, *fptr3, m;
	int i, nv, idx, nm;
	
	
	for (i=1; i<=n; i++) {a[i]=0.0; p[i]=0.;}
	nv = in_image->ifh.matrix_size[0]*in_image->ifh.matrix_size[1]*in_image->ifh.matrix_size[2]*in_image->ifh.matrix_size[3];
	fptr1=in_image->image;
	fptr2=Mask->image;
	fptr3=PETMask->image;
	nm=0;
	m=0;
	for (i=0; i<nv; i++) {
		if (*fptr3>0 && isnormal (*fptr1) && *fptr1 != (float) 1.e-37 && *fptr1 !=0){
			idx=(int)rint((double)(*fptr2))+1; /*belong to ROI # idx*/
			a[idx]+=*fptr1; /*sum for ROI #idx + voxel value*/
			p[idx]+=1.; /* voxel number for ROI #idx increase by 1 */
			if (idx>n) printf("something is wrong.\n");
			if (*fptr2 >0 ){
				nm++;
				m+=*fptr1;
			}
		}
		fptr1++;
		fptr2++;
		fptr3++;
	}
	m=m/(float)nm*.8;
	for (i=1; i<=n; i++) {
		if (p[i]>0) {
			a[i]=a[i]/p[i];
			/*printf("%f\t%f\n",p[i], a[i]);*/
		} else { /* to protect cases that ROI #i is completely out of PET field of view */
			a[i]=0; /* in such case, the ROI mean is replaced with mean value of all non-unknown regions multiplied by 0.8 */
		}
	}
}


RSFMat *calrsfmat(IMAGE_4dfp *image, ROIList *rois, float fhalf)
{
	RSFMat *rsfmat = NULL;
	IMAGE_4dfp *tmp_img=NULL;
	int val, nx, ny, nz, i, j, nv;
	float cmppix[3], *fptr, l;
	IFH tmp_ifh;
	clock_t start, end; 
	
	
	rsfmat = (RSFMat *) malloc((size_t)sizeof(RSFMat));
	if (!rsfmat) errm("calrsfmat");
	rsfmat->n = rois->NumberOfRegions;
	rsfmat->mat = matrix(1,rsfmat->n,1,rsfmat->n);
	nx = image->ifh.matrix_size[0];
	ny = image->ifh.matrix_size[1];
	nz = image->ifh.matrix_size[2];
	cmppix[0]=image->ifh.scaling_factor[0]/10.;
	cmppix[1]=image->ifh.scaling_factor[1]/10.;
	cmppix[2]=image->ifh.scaling_factor[2]/10.;
	
	/* Converting roi values to consecutive numbers*/
	start = clock();
	nv = image->ifh.matrix_size[0]*image->ifh.matrix_size[1]*image->ifh.matrix_size[2]*image->ifh.matrix_size[3];
	
	for (j=0; j<rois->NumberOfRegions; j++) {
		val = rois->List[j].MaskVal;
		fptr = image->image;
		for (i=0; i<nv; i++) {
			if (fabs(*fptr-(float)val)<1e-4) {
				*fptr=(float)j;
			}
			fptr++;
		}
	}
	fptr=image->image;
	l=(float)(rois->NumberOfRegions-.5);
	for (i=0; i<nv; i++) {
		if (*fptr>l) *fptr=0;
		fptr++;
	}
	
	end = clock();
	printf("CPU TIME USED for roi preprocessing is %f\n",((double) (end - start)) / CLOCKS_PER_SEC);
	
	/* Calculate RSF Matrix */
	for (i=1;i<=rsfmat->n;i++) {
		printf("%s\n", rois->List[i-1].Name);
		tmp_img=extract_roi_4dfp(image, i-1, tmp_img); 
		start = clock(); 
		gauss3d(tmp_img->image, &nx, &ny, &nz, cmppix, &fhalf);
		end = clock();
		printf("CPU TIME USED for gauss3d is %f\n",((double) (end - start)) / CLOCKS_PER_SEC);
		start = clock();
		getroismean_4dfp(tmp_img, image, rsfmat->n, rsfmat->mat[i]);
		/*for (j=1; j<=rsfmat->n; j++) {
			printf("%e\n",rsfmat->mat[i][j]);
		}*/		
		end = clock();
		printf("CPU TIME USED for getroimeans is %f\n",((double) (end - start)) / CLOCKS_PER_SEC);
	}
	
	
/*	
	for (i=1;i<=rsfmat->n;i++) {
		printf("%s\n", rois->List[i-1].Name);
		val = rois->List[i-1].MaskVal;
		tmp_img=extract_roi_4dfp(image, val, tmp_img); 
		start = clock(); 
		gauss3d(tmp_img->image, &nx, &ny, &nz, cmppix, &fhalf);
		end = clock();
		printf("CPU TIME USED for gauss3d is %f\n",((double) (end - start)) / CLOCKS_PER_SEC);
		start = clock();
		for (j=1; j<=rsfmat->n; j++) {
			rsfmat->mat[i][j]=getroimean_4dfp(tmp_img, image, rois->List[j-1].MaskVal);
			printf("%e\n",rsfmat->mat[i][j]);
		}
		end = clock();
		printf("CPU TIME USED for getroimean is %f\n",((double) (end - start)) / CLOCKS_PER_SEC);
	}
*/
	return rsfmat;
}

RSFROIS *Preprocess_RSF(IMAGE_4dfp *FS_Mask, IMAGE_4dfp *Head_Mask, IMAGE_4dfp *PET_Mask, RSFROIS *rsfrois, char *FSLUT)
/*
Combine Freesurfer defined masks, MR based head mask, and PET field of view based mask into a RSF Mask. It makes
sure that meaningful ROIs are only defined on those area within PET field of view.
*/
{
	int i, j, k, nv, tmp, xb, yb, zb, x, y, z, FSVal[16384], nfsroi, nroi, fsflag[16384], oflag[64], fsidx, val, maxfsval;
	float *fptr1, *fptr2, *fptr3, *fptr4, l;
	char **FSRegions = NULL;
	char line[512], otherstr[256];
	FILE *fp = NULL;
	
	printf("0\n");
	
	/* allocate memory and initiate rsfrois */
	nv=FS_Mask->ifh.number_of_bytes_per_pixel*FS_Mask->ifh.matrix_size[0]*FS_Mask->ifh.matrix_size[1]*FS_Mask->ifh.matrix_size[2]*FS_Mask->ifh.matrix_size[3];
	if (!rsfrois) {
		rsfrois=(RSFROIS *)malloc((size_t)sizeof(RSFROIS));
		if (!rsfrois) errm("Preprocess_RSF");
		rsfrois->RSFMask=(IMAGE_4dfp *)malloc((size_t)sizeof(IMAGE_4dfp));
		if (!rsfrois->RSFMask) errm("Preprocess_RSF");
		rsfrois->RSFMask->image=(float *)malloc((size_t)nv);
		if (!rsfrois->RSFMask->image) errm("Preprocess_RSF");	
		rsfrois->rois=(ROIList *)malloc((size_t)sizeof(ROIList));
		if (!rsfrois->rois) errm("Preprocess_RSF");
	} else {
		nrerror("rsfrois should be a NULL pointer passing to Preprocess_RSF");
	}
	nv = nv/FS_Mask->ifh.number_of_bytes_per_pixel;
	strcpy(rsfrois->RSFMask->ifh.interfile, FS_Mask->ifh.interfile);
	strcpy(rsfrois->RSFMask->ifh.version_of_keys, FS_Mask->ifh.version_of_keys);
	strcpy(rsfrois->RSFMask->ifh.conversion_program, FS_Mask->ifh.conversion_program);
	strcpy(rsfrois->RSFMask->ifh.name_of_data_file, "RSFMask");
	
	strcpy(rsfrois->RSFMask->ifh.number_format, FS_Mask->ifh.number_format);
	strcpy(rsfrois->RSFMask->ifh.imagedata_byte_order, FS_Mask->ifh.imagedata_byte_order);
	rsfrois->RSFMask->ifh.number_of_bytes_per_pixel=FS_Mask->ifh.number_of_bytes_per_pixel;
	rsfrois->RSFMask->ifh.number_of_dimensions=FS_Mask->ifh.number_of_dimensions;
	for (i=0; i<3; i++)
	{
		rsfrois->RSFMask->ifh.matrix_size[i]=FS_Mask->ifh.matrix_size[i];
		rsfrois->RSFMask->ifh.scaling_factor[i]=FS_Mask->ifh.scaling_factor[i];
		rsfrois->RSFMask->ifh.mmppix[i]=FS_Mask->ifh.mmppix[i];
		rsfrois->RSFMask->ifh.center[i]=FS_Mask->ifh.center[i];
	}
	rsfrois->RSFMask->ifh.matrix_size[3]=FS_Mask->ifh.matrix_size[3];
	rsfrois->RSFMask->ifh.scaling_factor[3]=FS_Mask->ifh.scaling_factor[3];
	rsfrois->RSFMask->ifh.orientation=FS_Mask->ifh.orientation;
	
	printf("1\n");
	/* load in freesurfer region definition*/
	FSRegions = (char **)malloc((size_t)(16384*sizeof(char *)));
	if (!FSRegions) errm("Preprocess_RSF");
	fp = fopen(FSLUT,"r");
	if (!fp) errr("Preprocess_RSF",FSLUT);
	maxfsval=0;
	for (i=0; i<16384; i++) fsflag[i]=0;
	while (fgets(line, 512, fp)){
		if (line[0]>='0'&&line[0]<='9') {
			sscanf(line, "%d%s%*s", &maxfsval, otherstr);
			FSRegions[maxfsval]=(char *)malloc((size_t)(256*sizeof(char)));
			if (!FSRegions[maxfsval]) errm("Preprocess_RSF");
			strcpy(FSRegions[maxfsval], otherstr);
			fsflag[maxfsval]=1;
			printf("%s\t%d\n",FSRegions[maxfsval],maxfsval);
			/*nfsroi++;*/
		}
	}
	fclose(fp);
	printf("maxfsval=%d\n",maxfsval);
	
	/* combine FS_Mask, Head_Mask and PET_Mask to form RSFMask*/
	xb=FS_Mask->ifh.matrix_size[0]/4;
	yb=FS_Mask->ifh.matrix_size[1]/4;
	zb=FS_Mask->ifh.matrix_size[2]/4;
	fptr1=FS_Mask->image;
	fptr2=Head_Mask->image;
	fptr3=rsfrois->RSFMask->image;
	fptr4=PET_Mask->image;
	l=0;
	for (i=0; i<64; i++) oflag[i]=0; /*set all "other region" flags to be 0, no "other region" is defined at this point*/
	for (z=0;z<FS_Mask->ifh.matrix_size[2];z++) {
		for (y=0; y<FS_Mask->ifh.matrix_size[1]; y++) {
			for (x=0; x<FS_Mask->ifh.matrix_size[0]; x++) {
				tmp=x/xb+y/yb*4+z/zb*16;
				if (*fptr1>0 && *fptr4>0) { /*for voxels defined in FS_Mask and within PET field of view*/
					if (l<*fptr1) l=*fptr1;
					*fptr3=*fptr1;
					fsidx=(int)rint((double)*fptr1);
					if (fsflag[fsidx]<=0) { /*the fs region is undefined*/
						if (*fptr2>0) { /*within Head Mask*/
							*fptr3=(float)(tmp+16384);
							oflag[tmp]=1;
						} else {/*outside of Head Mask*/
							*fptr3=0;
							fsflag[0]=2;
						}
					} else { /* fs region is defined */
						fsflag[fsidx]=2;
					}
				} else if (*fptr1<=0&&*fptr2>0&&*fptr4>0) { /*outside of fs defined region, within Head Mask*/
					*fptr3=(float)(tmp+16384);
					oflag[tmp]=1;
				} else {
					*fptr3=0;
					fsflag[0]=2;
				}
				fptr1++;
				fptr2++;
				fptr3++;
				fptr4++;
			}
		}
	}
	nroi=0;
	maxfsval=(int)rint((double)l);
	printf("maxfsval=%d\n",maxfsval);
	printf("fsflag[%d]=%d\n",1999,fsflag[1999]);
	for (i=0; i<=maxfsval; i++) {
		if (fsflag[i]>1) {nroi++; printf("%s\t%d\n",FSRegions[i], i);}
	}
	printf("nroi=%d\n",nroi);
	for (i=0; i<64; i++) {
		if (oflag[i]>0) nroi++;
	}
	rsfrois->rois->NumberOfRegions=nroi;
	rsfrois->rois->List=(ROI_Info *)malloc((size_t)(nroi*sizeof(ROI_Info)));
	if (!rsfrois->rois->List) errm("Preprocess_RSF");
	k=0;
	for (i=0; i<=maxfsval; i++) {
		if (fsflag[i]>1) {
			strcpy(rsfrois->rois->List[k].Name, FSRegions[i]);
			rsfrois->rois->List[k].MaskVal=i;
			k++;
		}
	}
	for (i=0; i<64; i++) {
		if (oflag[i]>0) {
			sprintf(otherstr,"other%d",i+16384);
			strcpy(rsfrois->rois->List[k].Name, otherstr);
			rsfrois->rois->List[k].MaskVal=i+16384;
			k++;
		}
	}
	printf("%d\n",rsfrois->rois->NumberOfRegions);
	for (j=0; j<rsfrois->rois->NumberOfRegions; j++) {
		val = rsfrois->rois->List[j].MaskVal;
		fptr3 = rsfrois->RSFMask->image;
		for (i=0; i<nv; i++) {
			if (fabs(*fptr3-(float)val)<1e-4) {
				*fptr3=(float)j;
			}
			fptr3++;
		}
	}
	fptr3=rsfrois->RSFMask->image;
	l=(float)(rsfrois->rois->NumberOfRegions)-.5;
	for (i=0; i<nv; i++) {
		if (*fptr3>l) *fptr3=0;
		fptr3++;
	}
	for (i=0; i<nroi; i++) {rsfrois->rois->List[i].MaskVal=i; rsfrois->rois->List[i].NVoxels=0;}
	fptr3=rsfrois->RSFMask->image;
	for (i=0; i<nv; i++) {
		val=(int)rint((double)*fptr3);
		rsfrois->rois->List[val].NVoxels++;
		fptr3++;
	}
	return rsfrois;
}

RSFROIS *Preprocess_RSF2(IMAGE_4dfp *FS_Mask, IMAGE_4dfp *Head_Mask, RSFROIS *rsfrois, char *FSLUT)
/*
Combine Freesurfer defined masks, MR based head maskinto a RSF Mask.
*/
{
	int i, j, k, nv, tmp, xb, yb, zb, x, y, z, FSVal[16384], nfsroi, nroi, fsflag[16384], oflag[64], fsidx, val, maxfsval;
	float *fptr1, *fptr2, *fptr3, l;
	char **FSRegions = NULL;
	char line[512], otherstr[256];
	FILE *fp = NULL;
	
	
	/* allocate memory and initiate rsfrois */
	nv=FS_Mask->ifh.number_of_bytes_per_pixel*FS_Mask->ifh.matrix_size[0]*FS_Mask->ifh.matrix_size[1]*FS_Mask->ifh.matrix_size[2]*FS_Mask->ifh.matrix_size[3];
	if (!rsfrois) {
		rsfrois=(RSFROIS *)malloc((size_t)sizeof(RSFROIS));
		if (!rsfrois) errm("Preprocess_RSF");
		rsfrois->RSFMask=(IMAGE_4dfp *)malloc((size_t)sizeof(IMAGE_4dfp));
		if (!rsfrois->RSFMask) errm("Preprocess_RSF");
		rsfrois->RSFMask->image=(float *)malloc((size_t)nv);
		if (!rsfrois->RSFMask->image) errm("Preprocess_RSF");	
		rsfrois->rois=(ROIList *)malloc((size_t)sizeof(ROIList));
		if (!rsfrois->rois) errm("Preprocess_RSF");
	} else {
		nrerror("rsfrois should be a NULL pointer passing to Preprocess_RSF");
	}
	nv = nv/FS_Mask->ifh.number_of_bytes_per_pixel;
	strcpy(rsfrois->RSFMask->ifh.interfile, FS_Mask->ifh.interfile);
	strcpy(rsfrois->RSFMask->ifh.version_of_keys, FS_Mask->ifh.version_of_keys);
	strcpy(rsfrois->RSFMask->ifh.conversion_program, FS_Mask->ifh.conversion_program);
	strcpy(rsfrois->RSFMask->ifh.name_of_data_file, "RSFMask");
	
	strcpy(rsfrois->RSFMask->ifh.number_format, FS_Mask->ifh.number_format);
	strcpy(rsfrois->RSFMask->ifh.imagedata_byte_order, FS_Mask->ifh.imagedata_byte_order);
	rsfrois->RSFMask->ifh.number_of_bytes_per_pixel=FS_Mask->ifh.number_of_bytes_per_pixel;
	rsfrois->RSFMask->ifh.number_of_dimensions=FS_Mask->ifh.number_of_dimensions;
	for (i=0; i<3; i++)
	{
		rsfrois->RSFMask->ifh.matrix_size[i]=FS_Mask->ifh.matrix_size[i];
		rsfrois->RSFMask->ifh.scaling_factor[i]=FS_Mask->ifh.scaling_factor[i];
		rsfrois->RSFMask->ifh.mmppix[i]=FS_Mask->ifh.mmppix[i];
		rsfrois->RSFMask->ifh.center[i]=FS_Mask->ifh.center[i];
	}
	rsfrois->RSFMask->ifh.matrix_size[3]=FS_Mask->ifh.matrix_size[3];
	rsfrois->RSFMask->ifh.scaling_factor[3]=FS_Mask->ifh.scaling_factor[3];
	rsfrois->RSFMask->ifh.orientation=FS_Mask->ifh.orientation;
	
	
	/* load in freesurfer region definition*/
	FSRegions = (char **)malloc((size_t)(16384*sizeof(char *)));
	if (!FSRegions) errm("Preprocess_RSF");
	fp = fopen(FSLUT,"r");
	if (!fp) errr("Preprocess_RSF",FSLUT);
	maxfsval=0;
	for (i=0; i<16384; i++) fsflag[i]=0;
	while (fgets(line, 512, fp)){
		if (line[0]>='0'&&line[0]<='9') {
			sscanf(line, "%d%s%*s", &maxfsval, otherstr);
			FSRegions[maxfsval]=(char *)malloc((size_t)(256*sizeof(char)));
			if (!FSRegions[maxfsval]) errm("Preprocess_RSF");
			strcpy(FSRegions[maxfsval], otherstr);
			fsflag[maxfsval]=1;
			/*printf("%s\t%d\n",FSRegions[maxfsval],maxfsval);*/
			/*nfsroi++;*/
		}
	}
	fclose(fp);
	printf("maxfsval=%d\n",maxfsval);
	
	/* combine FS_Mask, Head_Mask to form RSFMask*/
	xb=FS_Mask->ifh.matrix_size[0]/4;
	yb=FS_Mask->ifh.matrix_size[1]/4;
	zb=FS_Mask->ifh.matrix_size[2]/4;
	fptr1=FS_Mask->image;
	fptr2=Head_Mask->image;
	fptr3=rsfrois->RSFMask->image;
	l=0;
	for (i=0; i<64; i++) oflag[i]=0; /*set all "other region" flags to be 0, no "other region" is defined at this point*/
	for (z=0;z<FS_Mask->ifh.matrix_size[2];z++) {
		for (y=0; y<FS_Mask->ifh.matrix_size[1]; y++) {
			for (x=0; x<FS_Mask->ifh.matrix_size[0]; x++) {
				tmp=x/xb+y/yb*4+z/zb*16;
				if (*fptr1>0 ) { /*for voxels defined in FS_Mask*/
					if (l<*fptr1) l=*fptr1;
					*fptr3=*fptr1;
					fsidx=(int)rint((double)*fptr1);
					if (fsflag[fsidx]<=0) { /*the fs region is undefined*/
						if (*fptr2>0) { /*within Head Mask*/
							*fptr3=(float)(tmp+16384);
							oflag[tmp]=1;
						} else {/*outside of Head Mask*/
							*fptr3=0;
							fsflag[0]=2;
						}
					} else { /* fs region is defined */
						fsflag[fsidx]=2;
					}
				} else if (*fptr1<=0&&*fptr2>0) { /*outside of fs defined region, within Head Mask*/
					*fptr3=(float)(tmp+16384);
					oflag[tmp]=1;
				} else {
					*fptr3=0;
					fsflag[0]=2;
				}
				fptr1++;
				fptr2++;
				fptr3++;
			}
		}
	}
	nroi=0;
	maxfsval=(int)rint((double)l);
	printf("maxfsval=%d\n",maxfsval);
	printf("fsflag[%d]=%d\n",1999,fsflag[1999]);
	for (i=0; i<=maxfsval; i++) {
		if (fsflag[i]>1) {nroi++; printf("%s\t%d\n",FSRegions[i], i);}
	}
	printf("nroi=%d\n",nroi);
	for (i=0; i<64; i++) {
		if (oflag[i]>0) nroi++;
	}
	rsfrois->rois->NumberOfRegions=nroi;
	rsfrois->rois->List=(ROI_Info *)malloc((size_t)(nroi*sizeof(ROI_Info)));
	if (!rsfrois->rois->List) errm("Preprocess_RSF");
	k=0;
	for (i=0; i<=maxfsval; i++) {
		if (fsflag[i]>1) {
			strcpy(rsfrois->rois->List[k].Name, FSRegions[i]);
			rsfrois->rois->List[k].MaskVal=i;
			k++;
		}
	}
	for (i=0; i<64; i++) {
		if (oflag[i]>0) {
			sprintf(otherstr,"other%d",i+16384);
			strcpy(rsfrois->rois->List[k].Name, otherstr);
			rsfrois->rois->List[k].MaskVal=i+16384;
			k++;
		}
	}
	printf("%d\n",rsfrois->rois->NumberOfRegions);
	for (j=0; j<rsfrois->rois->NumberOfRegions; j++) {
		val = rsfrois->rois->List[j].MaskVal;
		fptr3 = rsfrois->RSFMask->image;
		for (i=0; i<nv; i++) {
			if (fabs(*fptr3-(float)val)<1e-4) {
				*fptr3=(float)j;
			}
			fptr3++;
		}
	}
	fptr3=rsfrois->RSFMask->image;
	l=(float)(rsfrois->rois->NumberOfRegions)-.5;
	for (i=0; i<nv; i++) {
		if (*fptr3>l) *fptr3=0;
		fptr3++;
	}
	for (i=0; i<nroi; i++) {rsfrois->rois->List[i].MaskVal=i; rsfrois->rois->List[i].NVoxels=0;}
	fptr3=rsfrois->RSFMask->image;
	for (i=0; i<nv; i++) {
		val=(int)rint((double)*fptr3);
		rsfrois->rois->List[val].NVoxels++;
		fptr3++;
	}
	return rsfrois;
}

float *RSFPVC(RSFMat *RSFMat, float *roimean, int iters)
/* 
RSFMat is a RSFMat structure that contains the matrix size n and a floating point matrix of n by n;
roimean is a vector that counts from 1 to n;
iters is the number of iterations to be performed using the iterative PVC algorithm;
The return value is a pointer that points to a vector that counts from 1 to n 
*/
{
	int i, j, k, NROI;
	float *r, *val, *m;
	
	NROI=RSFMat->n;
	r=vector(1,NROI);
	val=vector(1,NROI);
	m=vector(1,NROI);
	for (i=1; i<=NROI; i++) {m[i]=roimean[i];} /* Initial estimation */
	for (k=0; k<iters; k++)
	{
		for (i=1; i<=NROI; i++){
			val[i]=0.;
			for (j=1; j<=NROI; j++) {
				val[i]+=m[j]*RSFMat->mat[j][i]; /* Blurred regional value assuming current estimation */
			}
			r[i]=fabs(roimean[i]/val[i]); /* ratio between reblurred value and observed value */
			r[i]= (r[i]<1.2)?r[i]:1.2; /* constraints for the ratio to be applied to avoid blow up */
			r[i]= (r[i]>0.8)?r[i]:0.8;
		} 
		for (i=1; i<=NROI; i++) {
			m[i]=m[i]*r[i];	/* Apply correction for current iteration */
		}
	}
	
	free_vector(r, 1, NROI);
	free_vector(val, 1, NROI);
	return m;
}

void readtac(char *fn, float *tac, float *frd, float *st, float *nv, int *nframes)
{
	FILE *fp;
	char line[MAXL], dummy[MAXL];
	float fdum, *fptr1, *fptr2, *fptr3, *fptr4;
	int i;
	
	if (!(fp=fopen(fn, "r"))) errr("readtac","fn");
	fgets(line, 512, fp);
	sscanf(line, "%s%s%s%s%s%f", dummy, dummy, dummy, dummy, dummy, nv);
	fptr3=tac;
	fptr2=frd;
	fptr1=st;
	*nframes=0;
	while(fgets(line, 512, fp))
	{
		sscanf(line, "%f%f%f%f", &fdum, fptr1, fptr2, fptr3);
		fptr1++;
		fptr2++;
		fptr3++;
		(*nframes)++;
	}
	fclose(fp);
} 
