/**********************************************************************************
 Generation of whole brain ROI from FREESURFER segmentation results.
  
 USAGE:
 fs2brain FSfile fsROIlist
 FSfile is the 4dfp file of the freesurfer segmentation output, typically 
 wmparc001.4dfp.img; fsROIlist is three column text file that specifies the ROI 
 name and value corresponding to the FreeSurferColorLUT.txt file that contains all 
 the defined regions in the FSfile. The third column contains the voxel counts for
 each ROI.
 
 The output file has the name fsbrain, which is a mask of the regions belong to the
 brain, all other voxels are set to zero. 

 Yi Su, 12/24/2013
**********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>
#include <librms.h>
#include <nrutil.h>
#include <RSF.h>

int main(int argc, char **argv)
{
	int i, j, nvroi[16384], m2, vallist[16384], nroi, nx, ny, nz, tmp;
	float *fptr1, *fptr2, *fptr3, cmppix[3], fhalf, pvc[16384];
	IMAGE_4dfp *FSdata=NULL, *PETMask=NULL, *mask=NULL;
	FILE *fp;
	char **fsregions = NULL;
	char optf[MAXL], tmpstr[MAXL], line[512];

	/* read FSdata*/
	FSdata=read_image_4dfp(argv[1], FSdata);
	
	
	/* allocate memory for mask */
	m2=FSdata->ifh.number_of_bytes_per_pixel*FSdata->ifh.matrix_size[0]*FSdata->ifh.matrix_size[1]*FSdata->ifh.matrix_size[2]*FSdata->ifh.matrix_size[3];
	mask=(IMAGE_4dfp *)malloc((size_t)sizeof(IMAGE_4dfp));
	if (!mask) errm("fs2brain");
	mask->image=(float *)malloc((size_t)m2);
	if (!mask->image) errm("fs2brain");	
	m2=m2/FSdata->ifh.number_of_bytes_per_pixel;
	
	/* create binary mask that defines brain tissue against CSF and background */
	fptr1=FSdata->image;
	fptr2=mask->image;
	tmp=0;
	for (i=0; i<m2; i++)
	{
		if ( (int)(*fptr1)==0 || (int)(*fptr1)==4  || (int)(*fptr1)==5  || 
		    (int)(*fptr1)==14 || (int)(*fptr1)==15 || (int)(*fptr1)==30 ||
		    (int)(*fptr1)==43 || (int)(*fptr1)==44 || (int)(*fptr1)==62 ||
		    (int)(*fptr1)==72 || (int)(*fptr1)==73 || (int)(*fptr1)==74 ||
		    (int)(*fptr1)==75 || (int)(*fptr1)==76 || (int)(*fptr1)==24 ||
		    (int)(*fptr1)==98 || ((int)(*fptr1)>=118 && (int)(*fptr1)<=158) ||
		    (int)(*fptr1)==165 || (int)(*fptr1)==166 || (int)(*fptr1)==167 ||
		    (int)(*fptr1)==168 || ((int)(*fptr1)>=169 && (int)(*fptr1)<=215) ||
		    (int)(*fptr1)==217 || (int)(*fptr1)==221 || (int)(*fptr1)== 224 ||
		    ((int)(*fptr1)>=256 && (int)(*fptr1)<=499) || 
		    ((int)(*fptr1)>=559 && (int)(*fptr1)<=701) ||
		    ((int)(*fptr1)>=704 && (int)(*fptr1)<=999)) {
		    *fptr2=0;
		    }
		else {
		    *fptr2=1;
		    tmp++;
		    }
		fptr1++;
		fptr2++;
	}
	printf("%d\n",tmp);
	printf("Done...Binarize FS mask\n");

	strcpy(mask->ifh.interfile, FSdata->ifh.interfile);
	strcpy(mask->ifh.version_of_keys, FSdata->ifh.version_of_keys);
	strcpy(mask->ifh.conversion_program, FSdata->ifh.conversion_program);
	strcpy(mask->ifh.number_format, FSdata->ifh.number_format);
	strcpy(mask->ifh.imagedata_byte_order, FSdata->ifh.imagedata_byte_order);
	mask->ifh.number_of_bytes_per_pixel=FSdata->ifh.number_of_bytes_per_pixel;
	mask->ifh.number_of_dimensions=FSdata->ifh.number_of_dimensions;
	for (i=0; i<3; i++)
	{
		mask->ifh.matrix_size[i]=FSdata->ifh.matrix_size[i];
		mask->ifh.scaling_factor[i]=FSdata->ifh.scaling_factor[i];
		mask->ifh.mmppix[i]=FSdata->ifh.mmppix[i];
		mask->ifh.center[i]=FSdata->ifh.center[i];
	}
	mask->ifh.matrix_size[3]=FSdata->ifh.matrix_size[3];
	mask->ifh.scaling_factor[3]=FSdata->ifh.scaling_factor[3];
	mask->ifh.orientation=FSdata->ifh.orientation;

	strcpy(mask->ifh.name_of_data_file, "BrainMask");
	
	write_image_4dfp("BrainMask",mask);
	free(mask->image);
	free(mask);
	
	return 0;
}
