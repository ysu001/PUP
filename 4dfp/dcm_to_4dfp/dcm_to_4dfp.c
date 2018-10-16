/*
          Copyright (C) 1998-2013 Washington University

*/
/* Copyright marker.  Copyright will be inserted above.  Do not remove */
/*
**		     Electronic Radiology Laboratory
**		   Mallinckrodt Institute of Radiology
**		Washington University School of Medicine
**
** Module Name(s):	main
** Author, Date:	Stephen M. Moore, 9-Nov-1998
** Intent:		A set of DICOM images are placed in
**			analyze volume(s) (Analyze 7.5 hdr + img files)
**
*/

/******************************************************************************/

/* 04-27-2006 Modified by Mohana Ramaratnam 
**   to create only 4dfp file, corrected  
**   checkOrientation method and adapted for porting to Sparc on x86
*/

extern const char *gittag, *gitrev;
static char program[] = "dcm_to_4dfp";


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/fcntl.h>
#include <math.h>

#include "config.h"

#include "Getifh.h"
#include "endianio.h"
#include "rec.h"

#include "dicom.h"
#include "ctnthread.h"
#include "condition.h"
#include "lst.h"
#include "utility.h"
#include "dicom_objects.h"
#include "ifh.h"
#include "dcm_to_analyze.h"

static CTNBOOLEAN verbose = FALSE;
static CTNBOOLEAN useSeqFileName = FALSE;
static CTNBOOLEAN useFolder = FALSE;
static CTNBOOLEAN outputFileName = FALSE;


/* we want different linker behavior when building unit tests */
#ifdef BUILD_UNIT_TEST
#define STATIC
#define MAIN testmain
#else
#define STATIC static
#define MAIN main
#endif


static const int BUF_SIZE = 1024;
static const char *extension_4dfp = "4dfp";


/****************************************************************************/
/* Global Variables                                                         */

int 	*imageDim;
float	*voxelSize;
int 	 orientation = 2;  /* Default is Transverse */

/****************************************************************************/
typedef struct {
  void* reserved[2];
  char text[1024];
} SELECTION_ITEM;

/****************************************************************************/
static void
usageerror()
{

  printf ("%s %s (commit %s)\n", program, gittag, gitrev);
  printf ("Usage:\t%s [-b base] [-d gggg eeee] [-f] [-g] [-u] file(s)\n",
	  program);
  printf ("Slice Spacing Options: [-c] [-t <flt> or S or T]\n");
  printf ("Slice Position Options: [-X] [-Y] [-Z]\n");

  printf (" e.g.,\t%s *\n", program);
  printf ("   or,\t%s -b ID101 -f -g -u *IMA\n", program);
  printf ("   or,\t%s -d 0008 0030 -t 4.98 -g *.dcm\n", program);
  printf ("   or,\t%s -b P0089 -t T -g mydir/* \n", program);

printf ("Options:\n");

printf("\t[-b base] Output base filename follows the -b.\n");
printf("\t[-c]	    Slice Spacing: By Image Position (0020 0032).\n");
printf("\t[-d gggg eeee] Divide series by group and element number. \n");
printf("\t	 ** Default will divide volumes using ID series time (0008 0031). \n");
printf("\t[-f]	    Directories will be created, and dicom files will be moved. \n");
printf("\t[-g]	    Add image name, XYZ relative position, and number to rec file.\n\n");

printf("\t[-q]       Slice Spacing: Do not compute by Image Position.\n");
printf("\t[-r]       Rescale: Use the rescale slope and intercept fields.\n");
printf("\t[-t <flt>] Slice Spacing: Use input value.[-t <flt>]\n");
printf("\t[-t T]     Slice Spacing: Use Slice Thickness 0018 0050.[-t T]\n");
printf("\t[-t S]     Slice Spacing: Use Slice Spacing 0018 0088 [-t S](** Default)\n");

printf("\t[-u]	Output files named using sequence tag 0018 0024 plus number. \n\n");

printf("\t	4dfp Coordinant System is determined by Image Position (0020 0032). \n");
printf("\t	Multivolume and BOLD images are ordered by REL Image Number (0020 0013). \n");
printf("\t[-X]	Sagittal:	image positions will be ordered low to high\n");
printf("\t[-Y]	Coronal:	image positions will be high to low\n");
printf("\t[-Z]	Transverse:	image positions will be high to low\n");
printf("\t	** Default is transverse ordered by REL Image Number (0020 0013).\n");
printf("\t[-@ <b|l>]\toutput big or little endian (default CPU endian)\n");
/**printf ("\t-v  Debugging Statements: verbose mode\n");**/

exit (1);

}

/****************************************************************************/
STATIC CONDITION
parseParams(DCM_OBJECT** obj, const char* sliceThickness, PARAMS* p) {
  char imagePosition[DICOM_DS_LENGTH*3 + 32];
  char imageOrientation[DICOM_DS_LENGTH*6 + 32];
  char series_num_txt[DICOM_IS_LENGTH+1];
  char txt[DICOM_IS_LENGTH+1];
  char slope[DICOM_DS_LENGTH+1], intercept[DICOM_DS_LENGTH+1];

   /* tag  representation "descrip" mult length d->NULL */
  DCM_ELEMENT e[] = {
    { DCM_IMGROWS, DCM_US, "", 1, sizeof(p->rows), NULL},
    { DCM_IMGCOLUMNS, DCM_US, "", 1, sizeof(p->columns), NULL},
    { DCM_IMGBITSALLOCATED, DCM_US, "", 1, sizeof(p->bitsAllocated), NULL},
    { DCM_IMGBITSSTORED, DCM_US, "", 1, sizeof(p->bitsStored), NULL},
    { DCM_IMGHIGHBIT, DCM_US, "", 1, sizeof(p->highBit), NULL},
    { DCM_IMGPIXELREPRESENTATION, DCM_US, "", 1, sizeof(p->pixelRepresentation), NULL},
    { DCM_RELSERIESINSTANCEUID, DCM_UI, "", 1, sizeof(p->seriesUID), NULL},
    { DCM_RELIMAGENUMBER, DCM_IS, "", 1, sizeof(txt), NULL},
    { DCM_IMGPHOTOMETRICINTERP, DCM_CS, "", 1, sizeof(p->photometricInterpretation), NULL}
  };
  DCM_ELEMENT e0[] = {
    { DCM_RELSERIESNUMBER, DCM_IS, "", 1, sizeof(series_num_txt), NULL}}; /* (0020 0011) */
  DCM_ELEMENT e1[] = {
    { DCM_IMGPIXELSPACING, DCM_DS, "", 1, sizeof(p->pixelSpacing), NULL},/* (0028 0030) */
    { DCM_ACQSLICETHICKNESS, DCM_DS, "", 1, sizeof(p->sliceThickness), NULL },/* (0018 0050) */
    { DCM_RELIMAGEPOSITIONPATIENT, DCM_DS, "", 1, sizeof(imagePosition), NULL},/* (0020 0032) */
    { DCM_RELIMAGEORIENTATIONPATIENT, DCM_DS, "", 1, sizeof(imageOrientation), NULL}/* (0020 0037) */
  };
  DCM_ELEMENT e2[] = {
    { DCM_IMGPLANARCONFIGURATION, DCM_US, "", 1, sizeof(p->planarConfiguration), NULL}};
  DCM_ELEMENT e3[] = {
    { DCM_IDSTUDYDATE, DCM_DS, "", 1, sizeof(p->studyDate), NULL}};
  DCM_ELEMENT e4[] = {
    { DCM_IDSERIESDESCR, DCM_DS, "", 1, sizeof(p->seriesDescription), NULL}};
  DCM_ELEMENT e5[] = {
    { DCM_ACQPROTOCOLNAME, DCM_DS, "", 1, sizeof(p->protocolName), NULL }};
  DCM_ELEMENT e6[] = {
    { DCM_PATNAME, DCM_DS, "", 1, sizeof(p->patientName), NULL}};
  DCM_ELEMENT e7[] = {
    { DCM_PATID, DCM_DS, "", 1, sizeof(p->patientID), NULL}};
  DCM_ELEMENT e8[] = {
    { DCM_IMGLARGESTIMAGEPIXELVALUE, DCM_US, "", 1, sizeof(p->imgLargestPix), NULL},  /* 0028 0107 */
    { DCM_IMGSMALLESTIMAGEPIXELVALUE, DCM_US, "", 1, sizeof(p->imgSmallestPix), NULL} /* 0028 0106 */
  };
  DCM_ELEMENT e9[] = {
    { DCM_ACQSLICESPACING, DCM_DS, "", 1, sizeof(p->sliceSpacing), NULL}};	/* 0018 0088 */
  DCM_ELEMENT e10[] = {
    { DCM_ACQSEQUENCENAME, DCM_DS, "", 1, sizeof(p->acqSequenceName), NULL}};	/* 0018 0024 */
  DCM_ELEMENT e11[] = {
    { DCM_IDMANUFACTURER, DCM_DS, "", 1, sizeof(p->idManufacturer), NULL}};	/* 0008 0070 */
  DCM_ELEMENT e12[] = {
    { DCM_IDMANUFACTURERMODEL, DCM_DS, "", 1, sizeof(p->idModel), NULL}};	/* 0008 1090 */
  DCM_ELEMENT e13[] = {
    { DCM_IDMODALITY, DCM_DS, "", 1, sizeof(p->idModality), NULL}};		/* 0008 0060 */
  DCM_ELEMENT e14[] = {
    { DCM_IDINSTITUTIONNAME, DCM_DS, "", 1, sizeof(p->idInstitution), NULL}};	/* 0008 0080 */
  DCM_ELEMENT e15[] = {
    { DCM_IMGRESCALEINTERCEPT, DCM_DS, "", 1, sizeof(intercept), NULL}};	/* 0008 0080 */
  DCM_ELEMENT e16[] = {
    { DCM_IMGRESCALESLOPE, DCM_DS, "", 1, sizeof(slope), NULL}};	/* 0008 0080 */



  CONDITION cond;

  e[0].d.us = &p->rows;				/*  DCM_IMGROWS		*/
  e[1].d.us = &p->columns;			/*  DCM_IMGCOLUMNS	*/
  e[2].d.us = &p->bitsAllocated;		/*  DCM_IMGBITSALLOCATED*/
  e[3].d.us = &p->bitsStored;			/*  DCM_IMGBITSSTORED	*/
  e[4].d.us = &p->highBit;			/*  DCM_IMGHIGHBIT	*/
  e[5].d.us = &p->pixelRepresentation;		/*  DCM_IMGPIXELREPRESENTATION */
  e[6].d.string = p->seriesUID;			/*  DCM_RELSERIESINSTANCEUID */
  e[7].d.string = txt;				/*  DCM_RELIMAGENUMBER	*/
  e[8].d.string = p->photometricInterpretation;	/*  DCM_IMGPHOTOMETRICINTERP */

  e0[0].d.string = series_num_txt;	 /* DCM_RELSERIESNUMBER	*/

  e1[0].d.string = p->pixelSpacing;	/*  DCM_IMGPIXELSPACING		*/
  e1[1].d.string = p->sliceThickness;	/*  DCM_ACQSLICETHICKNESS	*/
  e1[2].d.string = imagePosition;	/*  DCM_RELIMAGEPOSITIONPATIENT	*/
  e1[3].d.string = imageOrientation;     /* DCM_RELIMAGEORIENTATIONPATIENT */

  e2[0].d.us = &p->planarConfiguration;
  e3[0].d.string =  p->studyDate;		/*  DCM_IDSTUDYDATE	*/
  e4[0].d.string =  p->seriesDescription;	/*  DCM_IDSERIESDESCR	*/
  e5[0].d.string =  p->protocolName;		/*  DCM_ACQPROTOCOLNAME	*/
  e6[0].d.string =  p->patientName;		/*  DCM_PATNAME		*/
  e7[0].d.string =  p->patientID;		/*  DCM_PATID		*/

  e8[0].d.us =  &p->imgLargestPix;	/*  DCM_IMGLARGESTIMAGEPIXELVALUE  */
  e8[1].d.us =  &p->imgSmallestPix;	/*  DCM_IMGSMALLESTIMAGEPIXELVALUE */

  e9[0].d.string = p->sliceSpacing;     /*  DCM_ACQSLICESPACING, DCM_DS */
  e10[0].d.string =  p->acqSequenceName;	/*  DCM_ACQSEQUENCENAME	    */
  e11[0].d.string =  p->idManufacturer;		/*  DCM_IDMANUFACTURER	    */
  e12[0].d.string =  p->idModel;		/*  DCM_IDMANUFACTURERMODEL */
  e13[0].d.string =  p->idModality;		/*  DCM_IDMODALITY	    */
  e14[0].d.string =  p->idInstitution;		/*  DCM_IDINSTITUTIONNAME   */
  e15[0].d.string =  intercept;		/*  DCM_IMGRESCALEINTERCEPT   */
  e16[0].d.string =  slope;		/*  DCM_IMGRESCALESLOPE   */

  cond = DCM_ParseObject(obj, e, (int)DIM_OF(e), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    printf ("Problem Parsing e in parseParams\n");
    COND_DumpConditions();
    return 1;
  }
  cond = DCM_ParseObject(obj, e0, (int)DIM_OF(e0), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
     printf ("Problem Parsing e0 Rel Ser Num in parseParams\n");
     (void)COND_PopCondition(TRUE);
  }
  cond = DCM_ParseObject(obj, e1, (int)DIM_OF(e1), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
     printf ("Problem Parsing e1 in parseParams\n");
     printf ("Check DCM_IMGPIXELSPACING 0028 0030\n");		/* p->pixelSpacing */
     printf ("Check DCM_ACQSLICETHICKNESS 0018 0050\n");	/* p->sliceThickness */
     printf ("Check DCM_RELIMAGEPOSITIONPATIENT 0020 0032\n");	/* imagePosition */
     printf ("Check DCM_RELIMAGEORIENTATIONPATIENT 0020 0037\n");	/* imageOrientation */
     /*(void)COND_PopCondition(TRUE);*/
     COND_DumpConditions();
     return 1;
  }
  {
    char *c;
    c = strtok(imagePosition, "\\");
    p->positionX = atof(c);
    c = strtok(NULL, "\\");
    p->positionY = atof(c);
    c = strtok(NULL, "\\");
    p->positionZ = atof(c);
 
   
    /* printf("Img Ori STring %s\n",imageOrientation); */
    c = strtok(imageOrientation, "\\");
    p->imageOrientation[0] = atof(c);
    c = strtok(NULL, "\\");
    p->imageOrientation[1] = atof(c);
    c = strtok(NULL, "\\");
    p->imageOrientation[2] = atof(c);
    c = strtok(NULL, "\\");
    p->imageOrientation[3] = atof(c);
    c = strtok(NULL, "\\");
    p->imageOrientation[4] = atof(c);
    c = strtok(NULL, "\\");
    p->imageOrientation[5] = atof(c);
   /* printf("Orientation Components %f %f %f %f %f %f\n",p->imageOrientation[0], p->imageOrientation[1], p->imageOrientation[2], p->imageOrientation[3],p->imageOrientation[4],p->imageOrientation[5]); */
  }

  if (sliceThickness != NULL && *sliceThickness != '\0'  ) { /* -t option used. */
     strcpy(p->sliceThickness, sliceThickness);
  }
  cond = DCM_ParseObject(obj, e2, (int)DIM_OF(e2), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
    p->planarConfiguration = 0;
  }
  cond = DCM_ParseObject(obj, e3, (int)DIM_OF(e3), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
  cond = DCM_ParseObject(obj, e4, (int)DIM_OF(e4), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
  cond = DCM_ParseObject(obj, e5, (int)DIM_OF(e5), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
  cond = DCM_ParseObject(obj, e6, (int)DIM_OF(e6), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
  cond = DCM_ParseObject(obj, e7, (int)DIM_OF(e7), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
	/* e8  is calculated, not parsed */

  cond = DCM_ParseObject(obj, e9, (int)DIM_OF(e9), NULL, 0, NULL);
      if (cond != DCM_NORMAL)
        (void)COND_PopCondition(TRUE);

  cond = DCM_ParseObject(obj, e10, (int)DIM_OF(e10), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
  cond = DCM_ParseObject(obj, e11, (int)DIM_OF(e11), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
   cond = DCM_ParseObject(obj, e12, (int)DIM_OF(e12), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
  cond = DCM_ParseObject(obj, e13, (int)DIM_OF(e13), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
  cond = DCM_ParseObject(obj, e14, (int)DIM_OF(e14), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
  cond = DCM_ParseObject(obj, e15, (int)DIM_OF(e15), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }
  cond = DCM_ParseObject(obj, e16, (int)DIM_OF(e16), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    (void)COND_PopCondition(TRUE);
  }

  p->instanceNumber = atoi(txt); /* REL image number 0020 0013 used to order images */
  p->seriesNumber = atoi(series_num_txt); /*  REL series number 0020 0011 */

  p->rescaleIntercept = atof(intercept);
  p->rescaleSlope = atof(slope);
  if (p->rescaleSlope == 0.0f)
  {
	p->rescaleSlope = 1.0f;
  }
  
  return 0;
}

/******************************************************************************/
#define BINARY_CHECK(a, b, str, f) if ((a) != (b)) { \
  fprintf(stderr, "Values do not match (%d %d) for %s in %s\n", (a), (b), (str), (f)); \
  exit(3); \
}

/******************************************************************************/
#define STRING_CHECK(a, b, str, f) if (strcmp((a), (b)) != 0) { \
  fprintf(stderr, "Values do not match (%s %s) for %s in %s\n", (a), (b), (str), (f)); \
  exit(3); \
}

/******************************************************************************/

STATIC
CONDITION testSelectionText(DCM_OBJECT** obj,
			    DCM_TAG tag, const char* txt)
{
  char buf[1024];
  CONDITION cond;
  DCM_ELEMENT e;

  memset(&e, 0, sizeof(e));
  e.tag = tag;
  e.d.string = buf;
  e.length = sizeof(buf) - 1;


  cond = DCM_ParseObject(obj, &e, 1, 0, 0, NULL);
  if (cond != DCM_NORMAL) {
    COND_DumpConditions();
    exit(1);
  }

  if (strncmp(txt, buf, strlen(txt)) == 0)
    return 0;

  return 1;
}


/**************************************************************************
 * testRangeText(obj, tag, lowerRange, upperRange)
 * Returns 0 if the indicated attribute is within the given range,
 * 1 otherwise
 **************************************************************************/
STATIC CONDITION
testRangeText(DCM_OBJECT** obj, DCM_TAG tag,
	      const char* lowerRange, const char* upperRange) {
  char buf[1024];
  DCM_ELEMENT e;
  long lowerVal, upperVal, testVal;

  memset(&e, 0, sizeof(e));
  e.tag = tag;
  e.d.string = buf;
  e.length = sizeof(buf) - 1;

  if (DCM_NORMAL != DCM_ParseObject(obj, &e, 1, 0, 0, NULL)) {
    COND_DumpConditions();
    exit(1);
  }

  testVal = atol(buf);
  lowerVal = atol(lowerRange);
  upperVal = atol(upperRange);

  return (lowerVal <= testVal && testVal <= upperVal) ? 0 : 1;
}

/******************************************************************************/
static void
readParamsForOneImage(const char* fileName,
		      DCM_TAG selectionTag, const char* selectionText,
		      DCM_TAG rangeTag, const char* lowerRange,
		      const char* upperRange,const char* sliceThickness,
		      LST_HEAD* l)
{
  DCM_OBJECT* obj;
  CONDITION cond;
  PARAMS *p;
  unsigned long options = DCM_ORDERLITTLEENDIAN;

  p = calloc(1, sizeof(*p));
  if (p == NULL) {
    fprintf(stderr, "Could not malloc structure for parameters\n");
    exit(9);
  }

  cond = DCM_OpenFile(fileName, DCM_ORDERLITTLEENDIAN, &obj);
  if (cond != DCM_NORMAL) {

    options &= ~DCM_FILEFORMATMASK;
    options |= DCM_PART10FILE;
    options |= DCM_ACCEPTVRMISMATCH;

    if (verbose)printf("readParamsForOneImage options = %lx\n",options);

    cond = DCM_OpenFile(fileName, options, &obj);

    if (cond != DCM_NORMAL) {
      COND_DumpConditions();
      exit(2);
    }
    COND_PopCondition(TRUE);
    p->part10Flag = 1;

  };

  if (selectionText != 0) {
    cond = testSelectionText(&obj, selectionTag, selectionText);
    if (cond != 0) {
      free(p);
      DCM_CloseObject(&obj);
    }
  }

  if (lowerRange != 0) {
    cond = testRangeText(&obj, rangeTag, lowerRange, upperRange);
    if (cond != 0) {
      free(p);
      DCM_CloseObject(&obj);
    }
  }

  cond = parseParams(&obj, sliceThickness, p);

  if (cond != 0) {
        printf ("Error in readParamsForOneImage %s\n",fileName);
        free(p);
        DCM_CloseObject(&obj);
  }


  strcpy(p->fileName, fileName);

  if (LST_Enqueue(l, p) != LST_NORMAL) {
       fprintf(stderr, "Could not enqueue list element\n");
       exit(10);
  }

  DCM_CloseObject(&obj);

}

/******************************************************************************/
static void scanDirectory(const char* dirName,
	DCM_TAG selectionTag, const char* selectionText,
	DCM_TAG rangeTag, const char* lowerRange, const char* upperRange,
			    const char* sliceThickness, LST_HEAD* l)
{
  LST_HEAD* localList = 0;
  CONDITION cond;
  UTL_FILEITEM* item;

  cond = UTL_ScanDirectory(dirName, &localList);
  if (cond != UTL_NORMAL) {
    COND_DumpConditions();
    return;
  }

  while ((item = (UTL_FILEITEM*)LST_Dequeue(&localList)) != NULL) {
    char path[1024];
    if (strcmp(item->path, ".") == 0)
      continue;
    if (strcmp(item->path, "..") == 0)
      continue;

    sprintf(path, "%s/%s", dirName, item->path);
    readParamsForOneImage(path, selectionTag, selectionText,
			  rangeTag, lowerRange, upperRange,
			  sliceThickness, l);
  }
}

/******************************************************************************/
static void readImageParams(int argc, char** argv,
	DCM_TAG selectionTag, const char* selectionText,
	DCM_TAG rangeTag, const char* lowerRange, const char* upperRange,
	const char* sliceThickness, LST_HEAD* l)
{
  DCM_OBJECT* obj;
  CONDITION cond;
  int count;
  PARAMS *p;
  unsigned long options;

  if (argc < 1)
    return;
  if (argc == 1) {
    if (UTL_IsDirectory(*argv)) {

      scanDirectory(*argv, selectionTag, selectionText,
		    rangeTag, lowerRange, upperRange,
		    sliceThickness, l);
      return;
    }
  }

  count = argc;
  while (count-- > 0) {
    if (verbose)
      printf("%s\n", *argv);
    if (strcmp(*argv, ".") == 0)
      continue;
    if (strcmp(*argv, "..") == 0)
      continue;

    p = calloc(1, sizeof(*p)); /* p is pointer to PARAMS */
    if (p == NULL) {
      fprintf(stderr, "Could not malloc structure for parameters\n");
      exit(9);
    }

    cond = DCM_OpenFile(*argv, DCM_ORDERLITTLEENDIAN, &obj);
    if (cond != DCM_NORMAL) {

      if (verbose){
          printf("DCM_ORDERLITTLEENDIAN=%d in readImageParams \n",DCM_ORDERLITTLEENDIAN);
          printf("FILEFORMATMASK=%x PART10FILE=%x ACCEPTVRMISMATCH=%x \n",\
                  DCM_FILEFORMATMASK,DCM_PART10FILE,DCM_ACCEPTVRMISMATCH);
      }
      options = DCM_ORDERLITTLEENDIAN;			/* 0x02 */
      options &= ~DCM_FILEFORMATMASK;			/* 0x80 */
         if (verbose)printf("options = %lx\n",options);	/* 0 */

      options |= DCM_PART10FILE;			/* 0x80 */
         if (verbose)printf("options = %lx\n",options);	/* 80 */

      options |= DCM_ACCEPTVRMISMATCH;			/* 0x4000 */
         if (verbose)printf("options = %lx\n",options);	/* 4080 */
         if (verbose)printf("NOW DCM_OpenFile *argv, options, &obj\n",options);

      cond = DCM_OpenFile(*argv, options, &obj);
      if (cond != DCM_NORMAL) {
        if (verbose)printf("Problem reading DCM_PART10FILE in readImageParams\n");

	COND_DumpConditions();
	exit(2);
      }

      COND_PopCondition(TRUE);
      p->part10Flag = 1;

    };

    if (selectionText != 0) {
	cond = testSelectionText(&obj, selectionTag, selectionText);
	if (cond != 0) {
	    free(p);
	    DCM_CloseObject(&obj);
	    argv++;
	    continue;
	}
    } else if (lowerRange != 0) {
	cond = testRangeText(&obj, rangeTag, lowerRange, upperRange);
	if (cond != 0) {
	    free(p);
	    DCM_CloseObject(&obj);
	    argv++;
	    continue;
	}
    }

    
    cond = parseParams(&obj, sliceThickness, p);
    if (cond != 0) {
            printf("parseParams Error in: %s\n",*argv);
	    free(p);
	    DCM_CloseObject(&obj);
	    argv++;
	    continue;
    }
    

    strcpy(p->fileName, *argv);

    if (LST_Enqueue(l, p) != LST_NORMAL) {
    	fprintf(stderr, "Could not enqueue list element\n");
    	exit(10);
    }

    DCM_CloseObject(&obj);
    argv++;
  }
}



/******************************************************************************/
static void setDimensions(char* base, LST_HEAD** l, char* spacingFlag)
{
	  PARAMS 	*p;
	  PARAMS 	p1;
	  int 		volCount = 1;
	  /*Volume number is now calculated in setSliceVol*/
	  char*		newbase;

	  memset(&p1, 0, sizeof(p1));
	  p = LST_Head(l);
	  LST_Position(l, p);
	  while(p != NULL) {
	     if (p1.rows == 0)
		      p1 = *p;

	     BINARY_CHECK(p1.rows, p->rows, "ROWS", p->fileName);
	     BINARY_CHECK(p1.columns, p->columns, "COLUMNS", p->fileName);
	     BINARY_CHECK(p1.bitsAllocated, p->bitsAllocated, "BITS ALLOCATED", p->fileName);
	     BINARY_CHECK(p1.bitsStored, p->bitsStored, "BITS STORED", p->fileName);
	     BINARY_CHECK(p1.highBit, p->highBit, "HIGH BIT", p->fileName);
	     BINARY_CHECK(p1.pixelRepresentation, p->pixelRepresentation, "PIXEL REP", p->fileName);
	     if (strcmp(spacingFlag, "T") == 0)
	        STRING_CHECK(p1.sliceThickness, p->sliceThickness, "T THICKNESS 0018 0050", p->fileName);
	     if (strcmp(spacingFlag, "S") == 0)
	        STRING_CHECK(p1.sliceSpacing, p->sliceSpacing, "S SPACING 0018 0088", p->fileName);

	     if (verbose)printf("setDimensions:: X %f %f\tY %f %f\tZ %f %f\n", \
	     p1.positionX, p->positionX, p1.positionY, p->positionY, p1.positionZ, p->positionZ );

	     p = LST_Next(l);
	  }
	  p = LST_Head(l);

	  

	  /* num dimensions -- for 4dfp = 4 */
	  imageDim[0] = p1.columns;	    /* xdim */
	  imageDim[1] = p1.rows;		    /* ydim */
	  imageDim[2] = LST_Count(l)/volCount; /* zdim: (nodes in LST)/# in seq */
	  					    /* imageDim[2] will be set in setSliceVol */
	  imageDim[3] = volCount;		    /* volumes: see setSliceVol */

	  /* In pixel spacing, first number is spacing between rows. Second is */
	  /* spacing between columns */
	  {
	    char str[DICOM_DS_LENGTH * 2 + 2];
	    char *c;
	    strcpy(str, p1.pixelSpacing);	/* tag 0028 0030 is divided into x and y */
	    c = strtok(str, "\\");
	    voxelSize[1] = atof(c);	/* ypix (column) pixel sizes (mm) */
	    c = strtok(NULL, "\\");
	    voxelSize[0] = atof(c);	/* xpix (row) pixel sizes (mm) */
	  }
		 				/* zpix pixel sizes (mm)dsr->dime.pixdim[3] */
	  					/* -t sets sliceThickness number in parseParams */

	  if (strcmp(spacingFlag, "I") == 0) {		/* Use -t <flt> for pixdim[3] */
	     voxelSize[2] = atof(p1.sliceThickness);
	     /*if (verbose)printf("I Input sliceThickness is %f\n",atof(p1.sliceThickness));*/

	  }else{					/* Check for selected tag or use default */

	     if (strcmp(spacingFlag, "T") == 0) {	/* Use -t T (0018 0050) */
	           voxelSize[2] = atof(p1.sliceThickness);
	           /*if (verbose)printf("T Flag: Use p1.sliceThickness for pixdim[3]\n");*/

	     }else if (p1.sliceSpacing[0] != '\0' || strcmp(spacingFlag, "S") == 0) {
	           voxelSize[2] = atof(p1.sliceSpacing);
	           /*if (verbose)printf("S Flag: SliceSpacing = %f\n",atof(p1.sliceSpacing));*/

	     }else if (p1.sliceThickness[0] != '\0') {
	          voxelSize[2] = atof(p1.sliceThickness);
	          /*if (verbose)printf("Default: SliceSpacing NULL and -t not used \n");*/

	     }else{
	          spacingFlag = "C"; /* p1.sliceThickness Is NULL so set Flag to compute */
	          /*if (verbose)printf("Compute: sliceSpacing NULL & sliceThickness NULL & -t Not Used \n");*/
	     }
	  }

	  if (verbose){
	      printf("\nsliceSpacing is %f\n",atof(p1.sliceSpacing));
	      printf("sliceThickness is %f\n",atof(p1.sliceThickness));
	      printf("pixdim[3]=%f\n", voxelSize[2]);
	      printf("Done setDimensions\n");
	  }

}



/******************************************************************************/
int xOrderhl(PARAMS *p1, PARAMS *p2)
{
  if (p1->positionX == p2->positionX) {
    return 0;
  }

  if (p1->positionX < p2->positionX)
    return 1;
  else
    return -1;
}

/******************************************************************************/
int xOrderlh(PARAMS *p1, PARAMS *p2)
{
  if (p1->positionX == p2->positionX) {
    return 0;
  }

  if (p1->positionX > p2->positionX)
    return 1;
  else
    return -1;
}

/******************************************************************************/
int yOrderhl(PARAMS *p1, PARAMS *p2)
{
  if (p1->positionY == p2->positionY) {
    return 0;
  }

  if (p1->positionY < p2->positionY)
    return 1;
  else
    return -1;
}

/******************************************************************************/
int yOrderlh(PARAMS *p1, PARAMS *p2)
{
  if (p1->positionY == p2->positionY) {
    return 0;
  }

  if (p1->positionY > p2->positionY)
    return 1;
  else
    return -1;
}

/******************************************************************************/
int zOrderhl(PARAMS *p1, PARAMS *p2)
{
  if (p1->positionZ == p2->positionZ) {
    return 0;
  }

  if (p1->positionZ < p2->positionZ)
    return 1;
  else
    return -1;
}

/******************************************************************************/
int zOrderlh(PARAMS *p1, PARAMS *p2)
{
  if (p1->positionZ == p2->positionZ) {
    return 0;
  }

  if (p1->positionZ > p2->positionZ)
    return 1;
  else
    return -1;
}

/******************************************************************************/
int instanceNumberOrder(PARAMS *p1, PARAMS *p2)
{
  if (p1->instanceNumber == p2->instanceNumber) {
    return 0;
  }

  if (p1->instanceNumber > p2->instanceNumber)
    return 1;
  else
    return -1;
}

/******************************************************************************/
int instanceNumberOrderhl(PARAMS *p1, PARAMS *p2)
{
  if (p1->instanceNumber == p2->instanceNumber) {
    return 0;
  }

  if (p1->instanceNumber < p2->instanceNumber)
    return 1;
  else
    return -1;
}

/******************************************************************************/
static void orderOnInstance(LST_HEAD** l, CTNBOOLEAN sortOnInstanceNumber)
{
  PARAMS *p;

  if (sortOnInstanceNumber) {
    LST_Sort(l, sizeof(*p), instanceNumberOrder);
    return;
  }

}

/****************************************************************************/
/*
**  FUNCTION:
**	orderXYZSlices
**  FUNCTIONAL DESCRIPTION:
**	The purpose of this function is to process an image volume that has
**	been sorted by REL Image Number, so that the number of volumes and
**	orientation could be determined. The image will be sorted again using
**	the REL image position to verify that the image is in the Analyze
**	Coordinate System. Multivolume or fmri data will be ordered by Image
**	Number. This ordering may be reversed with "xyx"Direction = TRUE (case XYX).
**  FORMAL PARAMETERS:
**	LST_HEAD** l
**	CTNBOOLEAN xDirection yDirection zDirection
**
**  RETURN VALUE:
**
**  SIDE EFFECTS:
*/

static void orderXYZSlices(LST_HEAD** l, CTNBOOLEAN xDirection,
	        CTNBOOLEAN yDirection, CTNBOOLEAN zDirection)
{
  PARAMS *p;
  PARAMS p1;

  int	VolCount =	imageDim[3];

  p = LST_Head(l);
  (void)LST_Position(l, p);
  p1 = *p;

  p = LST_Next(l);
  
  if (verbose) {
    printf("xDirection = %d yDirection=%d zDirection=%d\n",xDirection,yDirection,zDirection);
    printf("orderXYZSlices::Volcount = %d imageDim[3] = %d  orientation = %d \n",VolCount, imageDim[3], orientation);}
  if (VolCount == 1){
      if (orientation == 2){     /** Transverse default is Z Image Position low to high **/
         if (zDirection) { /**  Command line: -z set **/
             printf("\nReordering Transverse Z-axis high to low\n");
             LST_Sort(l, sizeof(*p), zOrderhl);
         } else {
             LST_Sort(l, sizeof(*p), zOrderlh);
         }
      }
      if (orientation == 3){     /** Coronal default is Y Image Position high to low **/
         if (yDirection) { /**  Command line: -y set **/
             printf("\nReordering Coronal Y-axis\n");
             LST_Sort(l, sizeof(*p), yOrderhl);
         } else {
             LST_Sort(l, sizeof(*p), yOrderlh);
         }
      }
      if (orientation == 4){     /** Sagittal default is X Image Position high to low **/
         if (xDirection) { /**  Command line: -x set **/
             printf("\nReordering Sagittal X-axis low to high\n");
             LST_Sort(l, sizeof(*p), xOrderlh);
         } else {
             LST_Sort(l, sizeof(*p), xOrderhl);
         }
      }
  } else {
      if (zDirection||yDirection||xDirection) { /** multivolume or fmri **/
	     /** Reverse REL Image Numbers to read high to low **/
             printf(" + \nReverse Ordering Images Across All Volumes\n");
             LST_Sort(l, sizeof(*p), instanceNumberOrderhl);
      } else {
             printf(" + ");	/** Default is low to high REL Image number**/
      }
  }

}
/*****************************************************************************/
static void writeImg(const char* base, LST_HEAD** l, char control, float* minVal, float* maxVal, CTNBOOLEAN useRescale)
{
  CONDITION cond;
  FILE *fd;
  void *fpixels = NULL;	/* pixel pointer */
  float* p1;
  PARAMS *p;
  DCM_OBJECT* obj;
  int pixelCount, cnt;
  int byteSize=2;
  int slicetot, volcount;
  char fileName[1024];

  sprintf(fileName, "%s.%s.%s", base, extension_4dfp, "img");
  if (!(fd = fopen(fileName, "wb")))  errw (program, fileName);
  
  printf("Writing %s\n",fileName);
  p = LST_Head(l);
  *minVal = 65535; *maxVal = -65536;

  (void)LST_Position(l, p);
  while (p != NULL) {

     unsigned long options = DCM_ORDERLITTLEENDIAN;
     if (p->part10Flag){
    	options &= ~DCM_FILEFORMATMASK;
	options |= DCM_PART10FILE;
	options |= DCM_ACCEPTVRMISMATCH;
     }

     if (verbose){
       printf("writeImg options = %lx\n",options);
       printf(" %s\n", p->fileName);
     }

     cond = DCM_OpenFile(p->fileName, options, &obj);
     if (cond != DCM_NORMAL) {
       COND_DumpConditions();
       exit(6);
     }

     pixelCount=p->rows * p->columns;

     fpixels = extractPixels(&obj, p);
     
     p1 = fpixels; cnt=pixelCount;
	 if (useRescale == TRUE && (p->rescaleSlope != 1.0f || p->rescaleIntercept != 0.0f))
	 {
		if (cnt) *minVal = *maxVal = p1[0] * p->rescaleSlope + p->rescaleIntercept;//dont segfault if no image data
		while(cnt-- > 0) {
			*p1 = *p1  * p->rescaleSlope + p->rescaleIntercept;//rescale
			if (*p1 > *maxVal)
				*maxVal = *p1;
			if (*p1 < *minVal)
				*minVal = *p1;
			p1++;
		}
	 } else {
		if (cnt) *minVal = *maxVal = p1[0];//dont segfault if no image data
		while(cnt-- > 0) {
			if (*p1 > *maxVal)
				*maxVal = *p1;
			if (*p1 < *minVal)
				*minVal = *p1;
			p1++;
		}
	 }

     if (verbose)printf("4dfp write p->fileName= %s\n",p->fileName);
     if (ewrite (fpixels,pixelCount,control,fd)) errw(program,fileName);

     free(fpixels);

     DCM_CloseObject(&obj);
     
     p = LST_Next(l);

  } /* end of while */

  fclose(fd); /* done writing img file */
}

/******************************************************************************/
static void setSliceVol(LST_HEAD** l)
{

  PARAMS 	*p, *p1;
  int 		slicenum = 	1;
  int 		volcount = 	1;
  int 		fmriAdjust;
  float 	Xpos, Ypos, Zpos;
  float 	deltaX=0, deltaY=0, deltaZ=0;

  /* This is where the number of volumes is determined.  */
  /* Correct orientation must be known. The volcount is  */
  /* incrimented if there is a change in image position  */
  /* AND is the same as the position of the first slice. */
  /* Images are assumed to be ordered by REL Image number.*/

  p = LST_Head(l);
  (void)LST_Position(l, p);
  p1 = LST_Next(l);

  Xpos = p1->positionX;
  Ypos = p1->positionY;
  Zpos = p1->positionZ;
  p1 = LST_Next(l);
  slicenum++;

  while (p1 != NULL && slicenum < imageDim[2]) {
    if(orientation == 4 && Xpos == p1->positionX && deltaX != 0) volcount++;
    if(orientation == 3 && Ypos == p1->positionY && deltaY != 0) volcount++;
    if(orientation == 2 && Zpos == p1->positionZ && deltaZ != 0) volcount++;

    deltaX = p1->positionX - p->positionX;
    deltaY = p1->positionY - p->positionY;
    deltaZ = p1->positionZ - p->positionZ;

    p = p1;
    slicenum++;
    p1 = LST_Next(l);
  }

  fmriAdjust = LST_Count(l);
  if (volcount == fmriAdjust) volcount = 1; /* fmri: Every image has same location */

  imageDim[3] = volcount;
  imageDim[2] = slicenum/volcount;

  if (verbose){
  	printf("\nThe List Count in setSliceVol= %d\n",fmriAdjust);
  	printf("Volcount in setSliceVol= %d\n",volcount);
  	printf("Slicenum in setSliceVol= %d\n",slicenum);
  	printf("imageDim[2] in setSliceVol= %d\n",slicenum/volcount);
  }

}

/******************************************************************************/
static void setSliceThickness(LST_HEAD** l)
{
  PARAMS 	*p, *p1;

  float 	meanSpacing = 	0.;
  float 	sdSpacing = 	0.;
  double 	sum = 		0.;
  double 	sumDSquared = 	0.;
  float 	divisor = 	1.;
  int 		slicenum;
  int 		volcount = 	1;
  float 	deltaX, deltaY, deltaZ;
  float 	delta;

  /* Now find the mean spacing between slices */

  p = LST_Head(l);
  (void)LST_Position(l, p);
  p1 = LST_Next(l);

  slicenum = 1;
  while (p1 != NULL && slicenum < imageDim[2]) {
    deltaX = p1->positionX - p->positionX;
    deltaY = p1->positionY - p->positionY;
    deltaZ = p1->positionZ - p->positionZ;
    delta = (deltaX*deltaX) + (deltaY*deltaY) + (deltaZ*deltaZ);
    delta = sqrt(delta);
    sum += delta;
    p = p1;

    slicenum++;
    p1 = LST_Next(l);
  }
  
  meanSpacing = sum / (imageDim[2] - 1);

  if (LST_Count(l) > 2) {
    p = LST_Head(l);
    (void) LST_Position(l, p);
    p1 = LST_Next(l);

    slicenum = 1;
    while (p1 != NULL && slicenum < imageDim[2]) {
      deltaX = p1->positionX - p->positionX;
      deltaY = p1->positionY - p->positionY;
      deltaZ = p1->positionZ - p->positionZ;
      delta = (deltaX*deltaX) + (deltaY*deltaY) + (deltaZ*deltaZ);
      delta = sqrt(delta);
      delta -= meanSpacing;
      sumDSquared += (delta*delta);

      p = p1;

      slicenum++;
      p1 = LST_Next(l);
    }

    divisor = imageDim[2] - 2;
    sdSpacing = sqrt(sumDSquared / divisor);

    if (sdSpacing > (meanSpacing/10.)) {
     fprintf(stderr, "SD of pixel spacing is large compared to mean spacing. \n");
     fprintf(stderr, "sd = %f, mean = %f\n", sdSpacing, meanSpacing);
     exit(1);
    }

    if (verbose) printf("meanSpacing = %f sd = %f\n", meanSpacing, sdSpacing);

    voxelSize[2] = meanSpacing;
  } else {
    fprintf(stderr, "Not enough slices to compute slice spacing \n");
    exit(1);
  }
}


/****************************************************************************/
static void writeIFH(const char *base, LST_HEAD** l,  char control)
{
  char fileName[BUF_SIZE];
  int chars_written;

  chars_written = snprintf(fileName, BUF_SIZE, "%s.%s.img",
			   base, extension_4dfp);
  assert(chars_written < BUF_SIZE);

  if (writeifhe(program, fileName, imageDim, voxelSize, orientation, control))
    errw (program, fileName);

}

/****************************************************************************/
static void
writeAnalyzeHeader(const char *base, const float *minVal, const float *maxVal)
{
  char command[BUF_SIZE];
  int chars_written;

  chars_written = snprintf(command, BUF_SIZE, "ifh2hdr %s -r%fto%f",
			   base, *minVal, *maxVal);
  assert(chars_written < BUF_SIZE);

  if (system(command))
    errw(program,command);
}

/******************************************************************************/

static void
writeREC(const char *base, LST_HEAD** l, const char* text,
	 int argc, char **argv,
	 CTNBOOLEAN zorder, CTNBOOLEAN xorder, CTNBOOLEAN yorder,
	 char control) {
  FILE 		*f;
  char 		fileName[1024];
  char 		imgName[1024];
  PARAMS 	*p;
  char 		command[1024];
  int 		sliceCounter = 1;
  int 		k;
  char 		extention[6];

  if (useFolder) {
     mkdir(base, 0777);
  }

  sprintf(fileName, "%s.%s.img", base, extension_4dfp);
  printf("Writing %s.rec\n",fileName);
  startrecle (fileName, argc, argv, gittag, control);
  printrec (((CPU_is_bigendian ()) ? !(control == 'l' || control == 'L') : (control == 'b' || control == 'B')) ? "bigendian\n" : "littleendian\n");

  p = LST_Head(l);
  (void)LST_Position(l, p);

  sprintf(command, "%s text:        \t%s\n", program, text);
  printrec(command);

  if(orientation == 2) printrec("Orientation:             \tTransverse\n");
  else if(orientation == 3) printrec("Orientation:             \tCoronal\n");
  else if(orientation == 4) printrec("Orientation:             \tSagittal\n");
  else if(orientation > 4){
    sprintf (command,"Series Division Text = %s\n",text);
    printrec(command);
  }

  sprintf (command, "ID Sequence Name:        \t%s\n",	p->acqSequenceName); printrec(command);
  sprintf (command, "ACQ Protocol Name:        \t%s\n",	p->protocolName); printrec(command);
  sprintf (command, "ACQ Sequence Description:\t%s\n\n",	p->seriesDescription); printrec(command);
  sprintf(command, "X dimension (columns): %6d      Width (mm): %10.5f \n",imageDim[0], voxelSize[0]); printrec(command);
  sprintf(command, "Y dimension (rows): %9d     Height (mm): %10.5f \n",  imageDim[1], voxelSize[1]); printrec(command);
  sprintf(command, "Z dimension (slices): %7d  Thickness (mm): %10.5f \n",imageDim[2], voxelSize[2]); printrec(command);
  sprintf(command, "Number of Volumes: %10d\nDimensions: %17d\n",imageDim[3], imageDim[3]); printrec(command);

  if(useFolder||outputFileName||xorder||yorder||zorder){  /** Add XYZ and sequence name to REC file **/
      if(xorder||yorder||zorder)
        sprintf(command, "Images Reordered - Note REL Image Position 0020 0037 or REL Image Number 0020 0013\n"); printrec(command);

      sprintf(command, "\t\t\t  X  \t\t   Y  \t\t   Z\t\t SEQ \t\t Img#\n"); printrec(command);

      while (p != NULL) {
         if(useFolder){
            sprintf(command,"%s/%s",base,p->fileName);
            rename(p->fileName, command);
         }
         sprintf(command, "FileName %d = %s  \t%12.6f \t%12.6f \t%12.6f \t %s \t %d\n",
            sliceCounter, p->fileName, p->positionX, p->positionY, p->positionZ,
	    p->acqSequenceName, p->instanceNumber ); printrec(command);

         p = LST_Next(l);
         sliceCounter++;
      }

  }else{  /** Otherwise list only the filenames in the REC file **/

      while (p != NULL) {

         sprintf(command, "FileName %d = %s\n",
              sliceCounter, p->fileName); printrec(command);

         p = LST_Next(l);
         sliceCounter++;
      }
  }

  /*if (p->pixelRepresentation == 1){
  ** Future Development: tags for PET images will be added here **
  */

  endrec();
}

/******************************************************************************/
int checkOrientation(LST_HEAD** l)
{

  PARAMS 	*p, *p1;
  float		normalX, normalY, normalZ;
  int 		orient;
  char          maxDirection;

  if (LST_Count(l) < 3) {
      printf("\nVolume Count < 3 ");
      orient= 5;
      return orient;
  }

  p = LST_Head(l);
  LST_Position(l, p);

  p1 = p;
  p = LST_Next(l);
  
  
  
    normalX=p->imageOrientation[1]*p->imageOrientation[5] - p->imageOrientation[2]*p->imageOrientation[4];
    normalY=p->imageOrientation[2]*p->imageOrientation[3] - p->imageOrientation[0]*p->imageOrientation[5];
    normalZ=p->imageOrientation[0]*p->imageOrientation[4] - p->imageOrientation[1]*p->imageOrientation[3];
    
    

  if (fabs(normalX)>fabs(normalY) && fabs(normalX) >fabs(normalZ))
    maxDirection = 'X';
  else if (fabs(normalY)>fabs(normalX) && fabs(normalY) >fabs(normalZ))
    maxDirection='Y';
  else if (fabs(normalZ)>fabs(normalX) && fabs(normalZ) >fabs(normalY))
    maxDirection ='Z';
  
  if(maxDirection=='X'){
           /* printf("    Sagittal  "); */
           orient = 4;
  }else if(maxDirection=='Y'){
           /* printf("    Coronal   "); */
           orient = 3;
  }else if(maxDirection=='Z'){
          /* printf("    Transverse  "); */
           orient = 2;
  }else{
       /* printf("    Default to Transverse "); */
       orient = 2;
  }
  if(verbose)
     printf("maxDirection=%c Orientation=%d\n",maxDirection, orient);
   
  return orient;
}

/******************************************************************************/
static LST_HEAD*
fillSelectionList(DCM_TAG tag, LST_HEAD** l) {
  PARAMS 	  *p;
  DCM_OBJECT* 	  obj;
  char 		  *s;
  SELECTION_ITEM* item;
  LST_HEAD* 	  rtnList;
  CONDITION 	  cond;

  if (verbose)
      printf("*IN fillSelectionList Ready for LST_Create \n");

  rtnList = LST_Create();

  if (verbose)
      printf("*IN fillSelectionList Ready for LST_Head\n");

  p = LST_Head(l);
  (void)LST_Position(l, p);
  while (p != NULL) {
    unsigned long options = DCM_ORDERLITTLEENDIAN;
    if (p->part10Flag){
    	options &= ~DCM_FILEFORMATMASK;
        options |= DCM_PART10FILE;
	options |= DCM_ACCEPTVRMISMATCH;
    }

    if (verbose){
      printf("fillSel options = %lx\n",options);
      printf(" %s\n", p->fileName);
    }
    cond = DCM_OpenFile(p->fileName, options, &obj);
    if (cond != DCM_NORMAL) {
      printf("Unable to fillSelectionList using: %s\n", p->fileName);
      COND_DumpConditions();
      exit(6);
    }
    if (verbose)
      printf("*IN fillSelectionList Finished opening file \n");

    s = DCM_GetString(&obj, tag);

    if (verbose)
      printf("*IN fillSelectionList Finished DCM_GetString \n");

    DCM_CloseObject(&obj);
    if (verbose)printf("* \n");

    item = LST_Head(&rtnList);
    LST_Position(&rtnList, item);

    while (item != NULL) {
      if (verbose)printf("*item is not NULL\n");
      if (strcmp(s, item->text) == 0) {
	break;
      }
      item = LST_Next(&rtnList);
    }

    if (item == NULL) {

      item = malloc(sizeof(*item));
      memset(item, 0, sizeof(*item));

      if (verbose)printf("*4 \n");
      strcpy(item->text, s);
      if (verbose)printf("*5 \n");

      LST_Enqueue(&rtnList, item);
      /* free(item); WILL Fail */
    }

    CTN_FREE(s);
    p = LST_Next(l);
  }

  if (verbose)printf("*Ready To Return from fillSelectionList\n");

  return rtnList;
}

/******************************************************************************/
static LST_HEAD*
selectImages(LST_HEAD** l, DCM_TAG tag, const char* text)

/* The purpose of this function is to identify dicom images in the list "l" using */
/* the DCM_TAG "tag" with it's "text". These images are placed in a new list */
/* "localCopy" which is returned. */

{
  PARAMS* 	pOrig;
  PARAMS* 	pCopy;
  char* 	s;
  CTNBOOLEAN	removeFlag = FALSE;
  CONDITION 	cond;
  DCM_OBJECT* 	obj;
  LST_HEAD* 	localCopy;

  localCopy = LST_Create();

  pOrig = LST_Head(l);
  LST_Position(l, pOrig);

  while (pOrig != NULL) {
    unsigned long options = DCM_ORDERLITTLEENDIAN;
    if (pOrig->part10Flag){
    	    options &= ~DCM_FILEFORMATMASK;
	    options |= DCM_PART10FILE;
	    options |= DCM_ACCEPTVRMISMATCH;
    }
    if (verbose)printf("selectImages options = %lx\n",options);

    cond = DCM_OpenFile(pOrig->fileName, options, &obj);
    if (cond != DCM_NORMAL) {
      COND_DumpConditions();
      printf("\nMust Exit\n");
      exit(6);
    }

    if(verbose)printf("Opened File %s\n",pOrig->fileName);
  
    s = DCM_GetString(&obj, tag);		/** place tag value of image in "s" **/
    DCM_CloseObject(&obj);

    if (strcmp(s, text) == 0) {			/** If tag values match **/
        pCopy = malloc(sizeof(*pCopy));
       *pCopy = *pOrig;
       LST_Enqueue(&localCopy, pCopy);	 	/** Add to localCopy list **/

       if(useSeqFileName || useFolder){ 	/** If the -u or -f option are selected **/
          pOrig = LST_Remove(l, LST_K_AFTER); 	/** shorten the "l" list. **/
          removeFlag = TRUE;			/** Set flag **/
       }
    }
    CTN_FREE(s);

    if((useSeqFileName || useFolder) && removeFlag){	/** Flag set AND -u or -f **/
        pOrig = LST_Current(l);				/** so point to Current **/
        removeFlag = FALSE;				/** and Reset flag **/
    }else{
        pOrig = LST_Next(l);				/** otherwise, point to Next **/
    }

  }

  free(pOrig);

  return localCopy;
}

/******************************************************************************/
/*
 * Construct a base name for files.  Starts from the sequence name,
 * then constructs a name from 1-5 letters (adding a leading 0 if
 * necessary) plus an underscore and the series number.  The resulting
 * base filename is not checked for validity or uniqueness.
 */
void setBase(LST_HEAD **l, char *out) {
  PARAMS *p;
  char *s;

  p = LST_Head(l);
  assert(p == LST_Position(l, p)); /* test LST package consistency */

  assert(out);		  /* caller must provide valid string space */

  s = p->acqSequenceName;
  assert(s);			/* sequence name must be non-NULL */
  while (*s == '*') s++;	/* skip past any leading asterisks */

  sprintf(out, "%1.5s_%d", s, p->seriesNumber);
  for (s = out; *s; s++) {
    if (' ' == *s) {
      *s = '0';
    }
  }
}


/******************************************************************************/

static void
processOneView(LST_HEAD** l, DCM_TAG tag, const char *text, char* base, int app,
	       CTNBOOLEAN sortOnInstanceNumber, int counter, CTNBOOLEAN xDirect,
	       CTNBOOLEAN yDirect, CTNBOOLEAN zDirect, char* spacingFlag,
	       CTNBOOLEAN computeSliceSpacing, CTNBOOLEAN useRescale,
		   int argc, char **argv, char control)

{
  LST_HEAD*	localCopy;
  PARAMS 	*p;
  char		sequenceName[1024];
  char		newbase[1024];
  float         minVal, maxVal;


  		if(verbose)printf("Ready to selectImages in processOneView\n");
  		/* Select images matching the tag and text and create localCopy */
  localCopy = selectImages(l, tag, text);

  		if(verbose)printf("Ready to Sort By Instance Number in processOneView\n");
  orderOnInstance(&localCopy, sortOnInstanceNumber); /* Sort by REL image number. */

  		if(verbose)printf("Ready to Create newbase in processOneView\n");

  if (useSeqFileName) {
      memset(sequenceName, '\0', sizeof(sequenceName));
      setBase(&localCopy, sequenceName);
      		if(verbose)printf("\nbase=%s  sequenceName=%s  app=%d\n",base,sequenceName,app);
      if(app == 0) {
         sprintf(newbase, "%s_%s", base, sequenceName);
      } else {
         sprintf(newbase, "%s_%s_%d", base, sequenceName, app);
      }
  }else{
      		if(verbose)printf("base=%s  app=%d\n", base, app);
      if(app == 0) {
          sprintf(newbase, "%s", base);
      } else {
          sprintf(newbase, "%s_%d", base, app);
      }
  }
  

  setDimensions(newbase, &localCopy, spacingFlag);	/* Initialize the hdr struc (-t) */
  orientation = checkOrientation(&localCopy);			/* Find orientation of sorted data */

 

  if (verbose) {
    printf("processOneView:: orientation = %d", orientation);
    printf("processOneView:: imageDims = %d %d %d %d\n",imageDim[0], imageDim[1], imageDim[2], imageDim[3]);
    printf("processOneView:: voxelDims = %d %d %d\n", voxelSize[0], voxelSize[1], voxelSize[2]);
  }
  if (orientation < 5){			/* orientation must be 2, 3, or 4 */
      		if(verbose)printf("Ready to setSliceVol in processOneView\n");
      setSliceVol(&localCopy);		/* Find number of volumes of sorted data */
      if(verbose){
	printf("Ready to order slices in processOneView\n");
      }		
      orderXYZSlices(&localCopy, xDirect, yDirect, zDirect);	/* Sort by REL image position */


      if (computeSliceSpacing)
          setSliceThickness(&localCopy);		/* Force computation of sliceThickness -q */

      writeImg(newbase, &localCopy, control, &minVal, &maxVal, useRescale);
      /* Now write the ifh and rec files */
      if(verbose)printf("Ready to ifh and rec files in processOneView\n");
      writeIFH(newbase, &localCopy, control);
      writeAnalyzeHeader(newbase,  &minVal, &maxVal);		
      writeREC(newbase, &localCopy, text, argc, argv, zDirect, xDirect, yDirect, control);

  } else {			/* The number of images in a volume is less than 3, */
       outputFileName = TRUE;	/* write a rec file with the file names */
       writeREC(newbase, &localCopy, text, argc, argv, zDirect, xDirect, yDirect, control);
  }
  printf("\n");

  while ((p = LST_Dequeue(&localCopy)) != NULL)
     free(p);
  LST_Destroy(&localCopy);

}

/******************************************************************************/

static void
processMultipleViews(LST_HEAD** l, DCM_TAG tag, CTNBOOLEAN sortOnInstanceNumber,
		     CTNBOOLEAN xDirection, CTNBOOLEAN yDirection,
		     CTNBOOLEAN zDirection, CTNBOOLEAN computeSliceSpacing,
			 CTNBOOLEAN useRescale,
		     char* base, char* spacingFlag, int argc, char **argv, char control)
{
  LST_HEAD* 	  vList = 0;
  SELECTION_ITEM* item;
  SELECTION_ITEM  *p;
  int 		  append = 0;
  char	 	  newBase[1024];
  int		  fileCounter = 1;

  if(verbose)printf("Ready to fillSelectionList in processMultipleViews\n");

  vList = fillSelectionList(tag, l);
  if (vList == 0)
    exit(1);

  if(verbose)printf("Done With fillSelectionList\n");

  item = LST_Head(&vList);
  LST_Position(&vList, item);
  while (item != NULL) {

    if(verbose)printf("Ready For processOneView \n");

    processOneView(l, tag, item->text, base, append,
		   sortOnInstanceNumber, fileCounter,
	           xDirection, yDirection, zDirection,
	           spacingFlag, computeSliceSpacing,
			   useRescale, argc, argv, control);

    item = LST_Next(&vList);
    fileCounter++;
    append++;
  }

  while ((p = LST_Dequeue(&vList)) != NULL)
       free(p);
  LST_Destroy(&vList);

}

/******************************************************************************/
/* main
**
** Purpose:
**	The purpose of this program is to convert dicom images into analyze
**	3D volumes.
**
** Parameter Dictionary:
**	None
**
** Return Values:
**	None
**
** Notes:
**     This program will now produce analyze images in floating point (4dfp)
**     format. Interfile ifh file, and record file are also written.
**
** Algorithm:
**	The original program is altered to produce analyze images in floating
**	point 4 dimensional format. Series are divided by the division tag.
**	Presently only floating point pixel values are produced.
**	Future implimentation will allow integer pixels (-n).
**	The -f option also increases the speed of the program by removing
**	files from the file list as they are written into analyze images.
**
*******************************************************************************/
int
MAIN(int argc, char **argv)
{
    CTNBOOLEAN xDirection = 		FALSE;
    CTNBOOLEAN yDirection = 		FALSE;
    CTNBOOLEAN zDirection = 		FALSE;
    CTNBOOLEAN computeSliceSpacing = 	FALSE;
    CTNBOOLEAN computeSliceSpace = 	TRUE;
    CTNBOOLEAN sortOnInstanceNumber = 	TRUE;
	CTNBOOLEAN useRescaleFields = 		FALSE;

    char 	**argvf=		argv;
    int		argcf=			argc;
    char 	*base = 		"analyze";
    char 	*sliceThickness = 	NULL;
    LST_HEAD* 	l;
    PARAMS 	*p;

    char* 	selectionText = 	0;
    int 	group = 0, element = 	0;
    DCM_TAG 	selectionTag = 		0;
    DCM_TAG 	rangeTag = 		0;
    char* 	lowerRange = 		0;
    char* 	upperRange = 		0;
    char* 	spacingFlag = 		"0";
    DCM_TAG 	divisionTag = DCM_MAKETAG(0x0020,0x000e); /* Default is Series Instance UID */
    char        control =               'O';


    while (--argc > 0 && (*++argv)[0] == '-') {

        if(verbose)printf("argc = %d argv = %s\n",argc, *argv);

		switch ((*argv)[1]) {

		case 'b': /* Output base filename follows the -b. The default is - analyze -  */
			argc--; argv++;
			base = *argv;
			break;
		case 'c': /* Slice Spacing: Compute using Image Position (0020 0032). (override all other choices) */
			computeSliceSpacing = TRUE;
			break;
		case 'd': /* Separate volumes by the header tag given by the group and element number. */
			if (argc < 2)
			usageerror();
			argc--; argv++;
			sscanf(*argv, "%x", &group);
			argc--; argv++;
			sscanf(*argv, "%x", &element);
			divisionTag = DCM_MAKETAG(group,element);
			break;
		case 'f': /* Directories will be created, and dicom files will be moved. */
			useFolder = TRUE;
			break;
		case 'g': /* Add image name, XYZ relative position, and number to rec file */
			outputFileName = TRUE;
			break;
		case 'n': /* Default is 4dfp Analyze orientation */
			/*floatFlag = FALSE;*/ /* Switch to integer pixels 4dint Analyze orientation */
			printf("Integer Pixels not available\n");
			break;
		case 'q': /* Slice Spacing: Do not Compute. ** default */
			computeSliceSpace = FALSE;
			break;
		case 'r':
			useRescaleFields = TRUE;
			break;
		case 't': /* Slice Spacing: Use number or 'S' or 'T' [-t <flt> S T ] */
				  /* 'S' : Use Slice Spacing 0018 0088 (the default if it exists) */
				  /* 'T' : Use Slice Thickness 0018 0050 if it exists*/

			argc--; argv++; sliceThickness = *argv; spacingFlag = "I";

			if(strcmp(sliceThickness,"T")==0){spacingFlag = "T"; *sliceThickness = '\0';}
			if(strcmp(sliceThickness,"S")==0){spacingFlag = "S"; *sliceThickness = '\0';}

			break;
		case 'u': /* Output files named using sequence (0018 0024). Rel series num appended. */
			useSeqFileName = TRUE;
			break;
		case 'v': /* verbose  */
			verbose = TRUE;
			break;
		case 'X': /* Sagittal: image positions will be ordered low to high */
			printf ("Sagittal- force low to high X order image positions\n");
			xDirection = TRUE;
			break;
		case 'Y': /* Coronal:	image positions will be high to low */
			printf ("Coronal- force high to low Y order image positions\n");
			yDirection = TRUE;
			break;
		case 'Z': /* Transverse: image positions will be high to low */
			printf ("Transverse- force high to low Z order image positions\n");
			zDirection = TRUE;
			break;
		case '@':
			argc--; argv++;
				control = *argv[0];
				break;
		default:
			break;
		}
    }

/****************************/
/* Initialize               */
/****************************/

    if (!(imageDim = (int *) malloc(4*sizeof(int)))) errm(program);
    if (!(voxelSize = (float *)malloc(3*sizeof(float)))) errm(program);

    THR_Init();
#if 0
    DCM_Debug(verbose);
#endif
    if (argc < 1)
	usageerror();

    l = LST_Create();

    if (l == NULL) {
      fprintf(stderr, "Could not create a new list\n");
      exit(8);
    }

    if(verbose)printf("Ready to readImageParams\n");

    readImageParams(argc, argv, selectionTag, selectionText,
		    rangeTag, lowerRange, upperRange,
		    sliceThickness, &l);

    

   if(verbose)printf("divisionTag = %d\n",divisionTag);

    processMultipleViews(&l, divisionTag, sortOnInstanceNumber,
			   xDirection, yDirection, zDirection,
			   computeSliceSpacing, useRescaleFields, base,
			   spacingFlag, argcf, argvf, control);

    while ((p = LST_Dequeue(&l)) != NULL)
      free(p);
    LST_Destroy(&l);

    THR_Shutdown();
    exit(0);
}
