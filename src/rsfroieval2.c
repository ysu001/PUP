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
	int i, Niter, NROI, f;
	char Frame4dfp[MAXL], RSFMask4dfp[MAXL], RSFMatfile[MAXL], PETMask4dfp[MAXL], ROIfile[MAXL], PETID[MAXL], f1[MAXL], f2[MAXL], line[512];
	float *fptr, *roimean, *rsfroimean, *roinv;
	IMAGE_4dfp *frame = NULL;
	IMAGE_4dfp *RSFMask = NULL;
	IMAGE_4dfp *PETMask = NULL;
	RSFMat *rsfmat;
	ROIList rois;
	ROI_Info *rip;
	FILE *fp, *fp2;
	
	strcpy(Frame4dfp, argv[1]);
	strcpy(RSFMask4dfp, argv[2]);
	strcpy(RSFMatfile, argv[3]);
	strcpy(PETMask4dfp, argv[4]);
	strcpy(ROIfile, argv[5]);
	NROI=atoi(argv[6]);
	Niter = atoi(argv[7]);
	strcpy(PETID, argv[8]);
	f = atoi(argv[9]);
	
	/* read PET frame, RSFMask, PETMask, RSFMat and rois */
	frame=read_image_4dfp(Frame4dfp, frame);
	RSFMask=read_image_4dfp(RSFMask4dfp, RSFMask);
	PETMask=read_image_4dfp(PETMask4dfp, PETMask);
	
	rsfmat = (RSFMat *) malloc((size_t)sizeof(RSFMat));
	if (!rsfmat) errm("rsfroieval");
	rsfmat->n = NROI;
	rsfmat->mat = matrix(1,NROI,1,NROI);
	
	fptr=rsfmat->mat[1];
	fptr++;
	if (!(fp =fopen(RSFMatfile, "r"))) errr("rsfroieval", RSFMatfile);
	while (fgets(line, 512, fp))
	{
		sscanf(line,"%e",fptr);
		fptr++;
	}
	fclose(fp);
	
	if (!(fp =fopen(ROIfile, "r"))) errr("rsfroieval", ROIfile);
	rois.List = (ROI_Info *)malloc((size_t)(NROI*sizeof(ROI_Info)));
	if (!rois.List) errm("rsfroieval");
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
	getroismean2_4dfp(frame, RSFMask, PETMask, NROI, roimean, roinv);
	
	/* Perform RSFPVC */
	rsfroimean=RSFPVC(rsfmat, roimean, Niter);
	
	/* Write output file */
	sprintf(f1, "%s_ROI2_f%d", PETID, f);
	sprintf(f2, "%s_RSF_ROI2_f%d", PETID, f);
	if (!(fp =fopen(f1, "w"))) errw("rsfroieval", f1);
	if (!(fp2 =fopen(f2, "w"))) errw("rsfroieval", f2);
	fprintf(fp, "%-35s %10s %16s\n", "Structure_Name", "NVoxels", "Mean Intensity");
	fprintf(fp2, "%-35s %10s %16s\n", "Structure_Name", "NVoxels", "Mean Intensity");

	for (i=0; i<NROI; i++){
		fprintf(fp, "%-35s %10d %16.4f\n", rois.List[i].Name, (int)roinv[i+1], roimean[i+1]);
		fprintf(fp2, "%-35s %10d %16.4f\n", rois.List[i].Name, (int)roinv[i+1], rsfroimean[i+1]);
	}
	fclose(fp);
	fclose(fp2);
}
