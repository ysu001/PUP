#ifndef MAXL
#define MAXL 256
#endif

typedef struct {
	int n;
	float **mat;
	} RSFMat;

typedef struct {
	IFH ifh;
	float *image;
	} IMAGE_4dfp;
	
typedef struct {
	char Name[MAXL];
	int NVoxels;
	int MaskVal;
	} ROI_Info;
	
typedef struct {
	int NumberOfRegions;
	ROI_Info *List;
	} ROIList;

typedef struct {
	IMAGE_4dfp *RSFMask;
	ROIList *rois;
	} RSFROIS;
		
IMAGE_4dfp *read_image_4dfp(char *filespc, IMAGE_4dfp *image);
int write_image_4dfp(char *filespc, IMAGE_4dfp *image);
float getroimean_4dfp(IMAGE_4dfp *in_image, IMAGE_4dfp *Mask, int val);
IMAGE_4dfp *copy_image_4dfp(IMAGE_4dfp *in_image, IMAGE_4dfp *out_image);
IMAGE_4dfp *extract_roi_4dfp(IMAGE_4dfp *image, int MaskVal, IMAGE_4dfp *roimask);
void getroismean_4dfp(IMAGE_4dfp *in_image, IMAGE_4dfp *Mask, int n, float *a);
void getroismean2_4dfp(IMAGE_4dfp *in_image, IMAGE_4dfp *Mask, IMAGE_4dfp *PETMask, int n, float *a, float *p);
void getroismean3_4dfp(IMAGE_4dfp *in_image, IMAGE_4dfp *Mask, IMAGE_4dfp *PETMask, int n, float *a, float *p);
RSFROIS *Preprocess_RSF(IMAGE_4dfp *FS_Mask, IMAGE_4dfp *Head_Mask, IMAGE_4dfp *PET_Mask, RSFROIS *rsfrois, char *FSLUT);
RSFROIS *Preprocess_RSF2(IMAGE_4dfp *FS_Mask, IMAGE_4dfp *Head_Mask, RSFROIS *rsfrois, char *FSLUT);
RSFMat *calrsfmat(IMAGE_4dfp *image, ROIList *rois, float fhalf);
float *RSFPVC(RSFMat *RSFMat, float *roimean, int iters);
IMAGE_4dfp *patlakvol(IMAGE_4dfp *petdata, IMAGE_4dfp *mask, float *t, float *aif, float *frd, int nframe, int stframe, IMAGE_4dfp *patlakout);
float loganREF(float *REF, float *ROI, float k2, float *t,  float *frd, int nframes, int stframe, int endframe, float *intercept, float *slope, float *Rsq);
bool file_exists(const char * filename);
void readtac(char *fn, float *tac, float *frd, float *st, float *nv, int *nframes);

