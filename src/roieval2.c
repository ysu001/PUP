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
	int i,  NROI, f;
	char Frame4dfp[MAXL], RSFMask4dfp[MAXL],  PETMask4dfp[MAXL], ROIfile[MAXL], PETID[MAXL], f1[MAXL], f2[MAXL], line[512];
	float *fptr, *roimean, *rsfroimean, *roinv;
	IMAGE_4dfp *frame = NULL;
	IMAGE_4dfp *RSFMask = NULL;
	IMAGE_4dfp *PETMask = NULL;
	ROIList rois;
	ROI_Info *rip;
	FILE *fp;
	
	strcpy(Frame4dfp, argv[1]);
	strcpy(RSFMask4dfp, argv[2]);
	strcpy(PETMask4dfp, argv[3]);
	strcpy(ROIfile, argv[4]);
	NROI=atoi(argv[5]);
	strcpy(PETID, argv[6]);
	f = atoi(argv[7]);
	
	/* read PET frame, RSFMask, PETMask,  and rois */
	frame=read_image_4dfp(Frame4dfp, frame);
	RSFMask=read_image_4dfp(RSFMask4dfp, RSFMask);
	PETMask=read_image_4dfp(PETMask4dfp, PETMask);
	
	
	if (!(fp =fopen(ROIfile, "r"))) errr("roieval2", ROIfile);
	rois.List = (ROI_Info *)malloc((size_t)(NROI*sizeof(ROI_Info)));
	if (!rois.List) errm("roieval2");
	rois.NumberOfRegions=NROI;
	i=0;
	rip=rois.List;
	while (fgets(line, 512, fp))
	{
		sscanf(line, "%s%d%d", rip->Name, &(rip->MaskVal), &(rip->NVoxels));
		printf("%s\t%d\t%d\n", rip->Name, rip->MaskVal, rip->NVoxels);
		rip++;
		i++;	
	}
	if (NROI!=i) nrerror("ROIfile and number of regions does not match");
	fclose(fp);
	
	/* Get raw roi means */
	roimean=vector(1,NROI);
	roinv=vector(1,NROI);
	getroismean3_4dfp(frame, RSFMask, PETMask, NROI, roimean, roinv);
	
	
	/* Write output file */
	sprintf(f1, "%s_ROI2_f%d", PETID, f);
	if (!(fp =fopen(f1, "w"))) errw("roieval2", f1);
	fprintf(fp, "%-35s %10s %16s\n", "Structure_Name", "NVoxels", "Mean Intensity");

	for (i=0; i<NROI; i++){
		fprintf(fp, "%-35s %10d %16.4f\n", rois.List[i].Name, (int)roinv[i+1], roimean[i+1]);
	}
	fclose(fp);
}
