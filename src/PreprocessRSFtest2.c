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

int main(int argc, char **argv)
{
	int i, j, n;
	char FS4dfp[MAXL], FSLUT[MAXL], Head4dfp[MAXL], RSF4dfp[MAXL], ROIfile[MAXL];
	IMAGE_4dfp *FS_Mask = NULL;
	IMAGE_4dfp *Head_Mask = NULL;
	
	RSFROIS *rsfrois = NULL;
	FILE *fp = NULL;
	
	strcpy(FS4dfp, argv[1]);
	strcpy(Head4dfp, argv[2]);
	
	strcpy(RSF4dfp, argv[3]);
	strcpy(ROIfile, argv[4]);
	strcpy(FSLUT, argv[5]);
	
	FS_Mask=read_image_4dfp(FS4dfp, FS_Mask);
	Head_Mask=read_image_4dfp(Head4dfp, Head_Mask);

	rsfrois=Preprocess_RSF2(FS_Mask, Head_Mask, rsfrois, FSLUT);
	strcpy(rsfrois->RSFMask->ifh.name_of_data_file, RSF4dfp);
	write_image_4dfp(RSF4dfp, rsfrois->RSFMask);
	
	fp=fopen(ROIfile, "w");
	if (!fp) errw("PreprocessRSFtest2",ROIfile);
	for (i=0;i<rsfrois->rois->NumberOfRegions;i++) {
		fprintf(fp,"%s\t%d\t%d\n",rsfrois->rois->List[i].Name, rsfrois->rois->List[i].MaskVal, rsfrois->rois->List[i].NVoxels);
	}
	fclose(fp);
}
