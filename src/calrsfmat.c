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


int main (int argc, char **argv)
{
	char imgfile[MAXL], roifile[MAXL], outfile[MAXL], line[512];
	int i, j;
	float fwhm;
	RSFMat *rsfmat =  NULL;
	ROIList rois;
	ROI_Info *rip;
	IMAGE_4dfp *image = NULL;
	FILE *fp;
	
	strcpy(imgfile, argv[1]);
	strcpy(roifile, argv[2]);
	strcpy(outfile, argv[3]);
	fwhm=atof(argv[4]);
	
	image=read_image_4dfp(imgfile, image);
	
	if (!(fp =fopen(roifile, "r"))) errr("calrsfmat_main", roifile);
	rois.List = (ROI_Info *)malloc((size_t)(256*sizeof(ROI_Info)));
	if (!rois.List) errm("calrsfmat_main");
	i=0;
	rip=rois.List;
	while (fgets(line, 512, fp))
	{
		sscanf(line, "%s%d%d", rip->Name, &(rip->MaskVal), &(rip->NVoxels));
		printf("%s\t%d\t%d\n", rip->Name, rip->MaskVal, rip->NVoxels);
		rip++;
		i++;	
	}
	rois.NumberOfRegions=i;
	rsfmat=calrsfmat(image, &rois, 4.412712/fwhm);
	fclose(fp);
	
	if (!(fp =fopen(outfile, "w"))) errr("calrsfmat_main", outfile);
	for (i=1; i<=rois.NumberOfRegions; i++)
	{
		for (j=1; j<=rois.NumberOfRegions; j++)
		{
			fprintf(fp, "%e\n", rsfmat->mat[i][j]);
		}
	}
	fclose(fp);
			
}
