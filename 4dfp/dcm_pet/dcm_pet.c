/*
$Header: /data/petsun4/data1/src_solaris/dcm_pet/RCS/dcm_pet.c,v 1.26 2014/05/05 16:11:25 jon Exp jon $ 
*/
/*
          Copyright (C) 2012 Washington University

*/
/* Copyright marker.  Copyright will be inserted above.  Do not remove */
/*
**		   Mallinckrodt Institute of Radiology
**		Washington University School of Medicine
**
** Authors:		SSM, Jon Christensen
** Intent:		PET DICOM images are converted to 4dfp Analyze volume(s)
**			using the Analyze 4dfp coordinant system.
**			 (Analyze 7.5 img, files hdr, with ifh and rec)
**
** Revision:		$Revision: 1.26 $
** Status:		$State: Exp $
**

$Log: dcm_pet.c,v $
Revision 1.26  2014/05/05 16:11:25  jon
-n option to calculate midpoint and start from DICOM frame duration

 * Revision 1.25  2013/08/28  15:02:39  jon
 * Meta Transfer Syntax added to rec file
 *
 * Revision 1.24  2013/06/05  21:25:44  jon
 * -a option will suppress listing more than the first 3 command line items in the rec file
 *
Revision 1.23  2013/06/05 15:46:53  jon
header and log for RCS added

*/
/******************************************************************************/
static char rcsid[] = "$Id: dcm_pet.c,v 1.26 2014/05/05 16:11:25 jon Exp jon $";
static char program[] = "dcm_pet";

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/fcntl.h>
#include <math.h>

#include "Getifh.h"
#include "endianio.h"
#include "rec.h"

#include "dicom.h"
#include "ctnthread.h"
#include "condition.h"
#include "lst.h"
#include "utility.h"
#include "dicom_objects.h"
#include "dcmprivate.h"
#include "dicom_sq.h"
#include "ifh.h"
#include "dcm_pet.h"

static	int verbose = 0;
static	const int BUF_SIZE = 1024;
static	CTNBOOLEAN useFileName = FALSE;
static	CTNBOOLEAN useFolder = FALSE;
static	CTNBOOLEAN outputFileName = FALSE;
static	CTNBOOLEAN midpoint_Frame_Ref = FALSE;
static	CTNBOOLEAN midpoint_Frame_Dur = FALSE;
static	CTNBOOLEAN suppress_listing_files = FALSE;

int	orientation;
int 	*imageDim;
float	*voxelSize;

#define STATIC static
#define MAIN main

/******************************************************************************/
typedef struct {
  void* reserved[2];
  char text[1024];
} SELECTION_ITEM;

/******************************************************************************/
static void
usageerror()
{

printf ("%s\n", rcsid);
printf ("Convert PET DICOM images, of same subject, to Analyze 4dfp.\n");
printf ("Usage:\tdcm_pet [-b base] [-d gggg eeee] [-t <flt> or S] [-v] <file(s)>\n");

printf (" e.g.,\tdcm_pet *\n");
printf ("   or,\tdcm_pet -b PIBCT101 *IMA\n");
printf ("   or,\tdcm_pet -d 0054 1300 petdir\n");

printf ("Options:\n");
printf("\tPET scan details are formated in the rec file. \n");
printf("\tPET images are scaled by rescale slope (0028 1053). \n");

printf("\t[-a] Suppress listing files in rec file. \n");

printf("\t[-b text] Output name will use text with numbering. \n");
printf("\tDefault output name will use analyze and number. \n");

printf("\t[-d gggg eeee] Divide PET volumes with DICOM tag. (Default is 0020 0011) \n");
printf("\ti.e. PET frames may be separated into 4dfp volumes with -d 0054 1300. \n");
printf("\t     Divide volumes by ID series time with -d 0008 0031. \n");

printf("\t[-m] Use 0054 1300 frame ref time (seconds) for rec file Midpoint and Start. \n");
printf("\t[-n] Use 0018 1242 frame duration (seconds) for rec file Midpoint and Start. \n");
printf("\t(Default will calculate Midpoint using injection time 0018 1072, 0008 0032, 0018 1242.) \n");

printf("\n\t           Slice Spacing: Default Will Use Slice Thickness 0018 0050.\n");
printf("\t[-t S]     Slice Spacing: Use Slice Spacing 0018 0088 [-t S]\n");
printf("\t[-t <flt>] Slice Spacing: Use input value. [-t <flt>]\n");

printf("\t[-g]       Debugging Option: Add file information to rec file.\n");
printf("\t[-v]       Debugging Option: verbose mode\n");

printf("\t[-@ <b|l>]\toutput big or little endian (default CPU endian)\n");

/* Disabled: */

/*printf("\t[-c]	Force computation of slice thickness. \n");			*/
/*printf("\t[-f]	Dir for DICOM in pwd created for each Analyze image. \n");	*/
/*printf("\t[-i]         					\n");			*/
/*printf("\t[-r gggg eeee lower upper]       Select images by tag and range\n");	*/
/*printf("\t[-s gggg eeee text]       Select images by tag and text\n");		*/
/*printf("\t[-u]	Output name will use series description and series number. \n");*/
/*printf("\t[-X]	Sagittal:	image positions will be ordered low to high\n");*/
/*printf("\t[-Y]	Coronal:	image positions will be high to low\n");	*/
/*printf("\t[-Z]	Transverse:	image positions will be high to low\n");	*/
/*printf("\t[-2]	Debugging Option: extreme verbosity\n");			*/

exit (1);

}

/******************************************************************************
**  FUNCTION: 
**	DumpElements
**  AUTHOR: Jon Christensen
**  FUNCTIONAL DESCRIPTION:
**	The purpose of this function is to parse sequenced header items of PET 
**	images. Adapted from the DCM_DumpElements debugging function in dcm.c.
**	The LST_HEAD d.sq from the DCM_ELEMENT sequence tag is sequenced.
**  FORMAL PARAMETERS:
**	The caller object must be a DCM object. "vm" is assumed to always be 1.
**	The PARAMS structure "p" is updated. Group and element 
**	tags are identified, matched and p is updated.
**  RETURN VALUE:
**	Recursive calls allow for nested sequences to be identified.
*/

CONDITION
DumpElements(DCM_OBJECT ** callerObject, long vm, PARAMS* p)
{
    PRV_GROUP_ITEM * groupItem;
    PRV_ELEMENT_ITEM * elementItem;
    PRIVATE_OBJECT ** object;
    CONDITION cond;
    DCM_SEQUENCE_ITEM * sq;
    char  scratch[128];
    char  temp[128];
    int   stringLength;
   
    object = (PRIVATE_OBJECT **) callerObject;

    groupItem = LST_Head(&(*object)->groupList);
    if (groupItem != NULL)
	(void) LST_Position(&(*object)->groupList, groupItem);

    while (groupItem != NULL) {
	elementItem = LST_Head(&groupItem->elementList);
	if (elementItem != NULL)
	    (void) LST_Position(&groupItem->elementList, elementItem);
	    
	while (elementItem != NULL) {

		sprintf(temp,"%04x %04x ", DCM_TAG_GROUP(elementItem->element.tag),
					   DCM_TAG_ELEMENT(elementItem->element.tag));
								   
		if(verbose == 2)printf("DumpElements: %s ", temp);
		if(verbose == 2)(void) printf("//%31s//", elementItem->element.description);

		if (elementItem->element.d.ot == NULL)
		   (void) printf("Warning: Data on disk\n");
	    else {
		
		switch (elementItem->element.representation) {
		case DCM_AE:
		case DCM_AS:
		case DCM_CS:
		case DCM_DA:
		case DCM_DT:
		    stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
		    strncpy(scratch, elementItem->element.d.string, stringLength);
		    if(strcmp(temp, "0008 0100 ") == 0)
		       strncpy(p->rad_dose, elementItem->element.d.string, stringLength);
		    scratch[stringLength] = '\0';
		
		    break;
		case DCM_DS:
		case DCM_IS:
		case DCM_LO:
		case DCM_LT:
		case DCM_PN:
		case DCM_SH:
		case DCM_UT:
		    stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
		    strncpy(scratch, elementItem->element.d.string, stringLength);
		    if(strcmp(temp, "0018 0031 ") == 0)
		       strncpy(p->rad, elementItem->element.d.string, stringLength);
		    if(strcmp(temp, "0018 1071 ") == 0)
		       strncpy(p->rad_volume, elementItem->element.d.string, stringLength);
		    if(strcmp(temp, "0018 1074 ") == 0)
		       strncpy(p->rad_dose, elementItem->element.d.string, stringLength);
		    if(strcmp(temp, "0018 1075 ") == 0)
		       strncpy(p->half_life, elementItem->element.d.string, stringLength);
		    if(strcmp(temp, "0018 1076 ") == 0)
		       strncpy(p->rad_fraction, elementItem->element.d.string, stringLength);
		    if(strcmp(temp, "0008 0100 ") == 0)
		       strncpy(p->id_code_value, elementItem->element.d.string, stringLength);
		    if(strcmp(temp, "0008 0102 ") == 0)
		       strncpy(p->id_code_scheme, elementItem->element.d.string, stringLength);
		    if(strcmp(temp, "0008 0104 ") == 0)
		       strncpy(p->id_code_meaning, elementItem->element.d.string, stringLength);
		    
		    scratch[stringLength] = '\0';

		    break;
		case DCM_SL:
#ifdef MACOS
		    if(verbose == 2)(void) printf("%8lx %ld\n", *elementItem->element.d.sl,
				  *elementItem->element.d.sl);
#else
		    if(verbose == 2)(void) printf("%8x %d\n", *elementItem->element.d.sl,
				  *elementItem->element.d.sl);
#endif
		    break;
		case DCM_SS:
		    if(verbose == 2)(void) printf("%4x %d\n", *elementItem->element.d.ss,
				  *elementItem->element.d.ss);
		    break;
		case DCM_SQ:
		    if(verbose == 2)printf("\n");
		    sq = LST_Head(&elementItem->element.d.sq);
		    if (sq != NULL)
			    (void) LST_Position(&elementItem->element.d.sq, sq);
		    while (sq != NULL) {
			    (void) DumpElements(&sq->object, vm, p);
			    sq = LST_Next(&elementItem->element.d.sq);
		    }
		    if (verbose == 2)printf("DCM Dump Sequence Complete\n");
		    break;
		case DCM_ST:
		    stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
		    strncpy(scratch, elementItem->element.d.string, stringLength);
		    
		    if(verbose == 2)(void) printf("%s\n", scratch);
		    scratch[stringLength] = '\0';
		
		    break;
		case DCM_TM:
		case DCM_UI:
		    stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
		    strncpy(scratch, elementItem->element.d.string, stringLength);
		    
		    if(strcmp(temp, "0018 1072 ") == 0)
		       strncpy(p->inject_time, elementItem->element.d.string, stringLength);
		    if(strcmp(temp, "0018 1073 ") == 0)
		       strncpy(p->rad_stop, elementItem->element.d.string, stringLength);

		    scratch[stringLength] = '\0';
		    
		    break;
		case DCM_AT:
		case DCM_UL:
#ifdef MACOS
		    if(verbose == 2)(void) printf("%8lx %ld\n", *elementItem->element.d.ul,
				  *elementItem->element.d.ul);
#else
		    if(verbose == 2)(void) printf("%8x %d\n", *elementItem->element.d.ul,
				  *elementItem->element.d.ul);
#endif
		    break;
		case DCM_US:
		    if(verbose == 2)(void) printf("%4x %d\n", *elementItem->element.d.us,
				  *elementItem->element.d.us);
		    break;
		case DCM_OB:
		case DCM_UN:
		    break;

		case DCM_OT:
		case DCM_OW:
		
		default:
		    (void) printf("Warning: Unimplemented VR\n");
		    break;
		}
	    }
	    elementItem = LST_Next(&groupItem->elementList);
	}
	groupItem = LST_Next(&(*object)->groupList);
    }
    return DCM_NORMAL;
}

/******************************************************************************
**  FUNCTION: 
**	parseParams
**  FUNCTIONAL DESCRIPTION:
**	The purpose of this function is to parse the DCM_OBJECT "obj" and load 
**	values into the PARAMS structure "p". DCM_ELEMENT structures are initialized,
**	pointers to the p struc are set and tag values are then parsed into 
**	the DCM_ELEMENT with DCM_ParseObject. The sequenced header tags are parsed 
**	by calls to DumpElements.
*/
static CONDITION parseParams(DCM_OBJECT** obj, const char* sliceThickness,
			PARAMS* p)
{
  char	imagePosition[DICOM_DS_LENGTH*3 + 10];
  char	imageOrientation[DICOM_DS_LENGTH*6 + 10];
  char	petSeriesType[DICOM_DS_LENGTH*2 + 10];
  char	series_num_txt[DICOM_IS_LENGTH+1];
  char	frametxt[DICOM_IS_LENGTH+1] = "0";
  char	frdurtxt[DICOM_IS_LENGTH+1] = "0";
  char	frame_reftxt[DICOM_DS_LENGTH+1] = "0";
  char	imagetxt[DICOM_IS_LENGTH+1] = "0";
  int	PET,MR,CT,SR;

  /* tag  representation "descrip" mult length d->NULL */
  DCM_ELEMENT e[] = {
    { DCM_IMGROWS, DCM_US, "", 1, sizeof(p->rows), NULL},				/* 0028 0010 */
    { DCM_IMGCOLUMNS, DCM_US, "", 1, sizeof(p->columns), NULL},				/* 0028 0011 */
    { DCM_IMGBITSALLOCATED, DCM_US, "", 1, sizeof(p->bitsAllocated), NULL},		/* 0028 0100 */
    { DCM_IMGBITSSTORED, DCM_US, "", 1, sizeof(p->bitsStored), NULL},			/* 0028 0101 */
    { DCM_IMGHIGHBIT, DCM_US, "", 1, sizeof(p->highBit), NULL},				/* 0028 0102 */
    { DCM_IMGPIXELREPRESENTATION, DCM_US, "", 1, sizeof(p->pixelRepresentation), NULL},	/* 0028 0103 */
    { DCM_RELSERIESINSTANCEUID, DCM_UI, "", 1, sizeof(p->seriesUID), NULL},		/* 0020 000e */
    { DCM_RELIMAGENUMBER, DCM_IS, "", 1, sizeof(imagetxt), NULL}     			/* 0020 0013 */
  };
  DCM_ELEMENT e0[] = {
    { DCM_RELSERIESNUMBER, DCM_IS, "", 1, sizeof(series_num_txt), NULL},           	/* 0020 0011 */
    { DCM_IMGPHOTOMETRICINTERP, DCM_CS, "", 1, sizeof(p->photometricInterpretation), NULL}/* 0028 0004 */
  };
  DCM_ELEMENT e1[] = {
    { DCM_IMGPIXELSPACING, DCM_DS, "", 1, sizeof(p->pixelSpacing), NULL},           /* 0028 0030 */
    { DCM_ACQSLICETHICKNESS, DCM_DS, "", 1, sizeof(p->sliceThickness), NULL },      /* 0018 0050 */
    { DCM_RELIMAGEPOSITIONPATIENT, DCM_DS, "", 1, sizeof(imagePosition), NULL},     /* 0020 0032 */
    { DCM_RELIMAGEORIENTATIONPATIENT, DCM_DS, "", 1, sizeof(imageOrientation), NULL}/* 0020 0037 */
  };
  DCM_ELEMENT e2[] = {
    { DCM_IMGPLANARCONFIGURATION, DCM_US, "", 1, sizeof(p->planarConfiguration), NULL}};/*0028 0006*/
  DCM_ELEMENT e3[] = {
    { DCM_IDSTUDYDATE, DCM_DS, "", 1, sizeof(p->studyDate), NULL}};		/* 0008 0020 */
  DCM_ELEMENT e4[] = {
    { DCM_IDSERIESDESCR, DCM_DS, "", 1, sizeof(p->seriesDescription), NULL}};	/* 0008 103e */
  DCM_ELEMENT e5[] = {
    { DCM_ACQPROTOCOLNAME, DCM_DS, "", 1, sizeof(p->protocolName), NULL }};	/* 0018 1030 */
  DCM_ELEMENT e6[] = {
    { DCM_PATNAME, DCM_DS, "", 1, sizeof(p->patientName), NULL}};		/* 0010 0010 */
  DCM_ELEMENT e7[] = {
    { DCM_PATID, DCM_DS, "", 1, sizeof(p->patientID), NULL}};			/* 0010 0020 */
  DCM_ELEMENT e8[] = { 
    { DCM_IMGLARGESTIMAGEPIXELVALUE, DCM_DS, "", 1, sizeof(p->imgLargestPix), NULL }, /* 0028 0107 */
    { DCM_IMGSMALLESTIMAGEPIXELVALUE, DCM_DS, "", 1, sizeof(p->imgSmallestPix), NULL} /* 0028 0106 */
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
     { DCM_IDMODALITY, DCM_DS, "", 1, sizeof(p->idModality), NULL },		/* 0008 0060 */
     { DCM_METATRANSFERSYNTAX, DCM_UI, "", 1, sizeof(p->metaTransferUID), NULL} /* 0002 0010 */
   };
  DCM_ELEMENT e14[] = {
    { DCM_IDINSTITUTIONNAME, DCM_DS, "", 1, sizeof(p->idInstitution), NULL}};	/* 0008 0080 */
  
    /* PET */
  DCM_ELEMENT e15[] = {
    {DCM_NMISERIESTYPE, DCM_CS, "", 1, sizeof(petSeriesType), NULL}};		/* 0054 1000 */
  DCM_ELEMENT e16[] = {
    {DCM_IMGCORRECTEDIMAGE, DCM_CS, "", 1, sizeof(p->image_correction), NULL}};	/* 0028 0051 */
  DCM_ELEMENT e17[] = {
    {DCM_NMIIMAGEINDEX, DCM_US, "", 1, sizeof(p->nmi_image_index), NULL}};	/* 0054 1330 */
  DCM_ELEMENT e18[] = {
    {DCM_NMICOUNTSSOURCE, DCM_CS, "", 1, sizeof(p->counts_source), NULL}};	/* 0054 1002 */
  DCM_ELEMENT e19[] = {
    {DCM_NMINUMBEROFSLICES, DCM_US, "", 1, sizeof(p->nmi_slice_max), NULL}};	/* 0054 0081 */
  DCM_ELEMENT e20[] = {
    {DCM_IMGNUMBEROFFRAMES, DCM_IS, "", 1, sizeof(frametxt), NULL}};		/* 0028 0008 */
  DCM_ELEMENT e21[] = {
    {DCM_NMIATTENUATIONCORRECTIONMETHOD, DCM_LO, "", 1, sizeof(p->atten_correction), NULL}};/* 0054 1101 */
  DCM_ELEMENT e22[] = {
    {DCM_NMISCATTERCORRECTIONMETHOD, DCM_LO, "", 1, sizeof(p->scatt_correction), NULL}};/* 0054 1105 */
  DCM_ELEMENT e23[] = {
    {DCM_NMIDECAYCORRECTION, DCM_CS, "", 1, sizeof(p->decay_correction), NULL}};/* 0054 1102 */
  /* DCM_ELEMENT e24[] = {*/
  /*   {DCM_A , DCM_DS, "", 1, sizeof(p-> ), NULL}};*/	/*   */
  /*DCM_ELEMENT e25[] = {*/
  /*   {DCM_A , DCM_DS, "", 1, sizeof(p-> ), NULL}};*/	/*   */
  DCM_ELEMENT e27[] = {
    {DCM_IDSERIESTIME, DCM_TM, "", 1, sizeof(p->series_time), NULL}};		/* 0008 0031 */
  DCM_ELEMENT e28[] = {
    {DCM_IDACQUISITIONTIME, DCM_TM, "", 1, sizeof(p->acq_start_time), NULL}};	/* 0008 0032 */
  DCM_ELEMENT e29[] = {
    {DCM_NMIFRAMEREFERENCETIME, DCM_DS, "", 1, sizeof(frame_reftxt), NULL}};	/* 0054 1300 */
  DCM_ELEMENT e30[] = {
    {DCM_ACQACTUALFRAMEDURATION, DCM_IS, "", 1, sizeof(frdurtxt), NULL}};	/* 0018 1242 */
  DCM_ELEMENT e31[] = {
    {DCM_IMGRESCALEINTERCEPT, DCM_AS, "", 1, sizeof(p->rescale_intercept), NULL}};/* 0028 1052 */
  DCM_ELEMENT e32[] = {
    {DCM_IMGRESCALESLOPE, DCM_DS, "", 1, sizeof(p->rescale_slope), NULL}};	/* 0028 1053 */
  DCM_ELEMENT e33[] = {
    {DCM_IMGLOSSYIMAGECOMPRESSION, DCM_AS, "", 1, sizeof(p->lossy_comp), NULL}}; /* 0028 2110 */
  DCM_ELEMENT e34[] = {
    {DCM_ACQRADIOPHARMVOLUME, DCM_DS, "", 1, sizeof(p->rad_volume), NULL}};	/* 0018 1071 */
  DCM_ELEMENT e35[] = {
    {DCM_NMIUNITS, DCM_CS, "", 1, sizeof(p->nmi_units), NULL}};			/* 0054 1001 */
    /* DCM_ELEMENT e36 */
    /* DCM_ELEMENT e37 */
  DCM_ELEMENT e38[] = {
    {DCM_NMIRECONSTRUCTIONMETHOD, DCM_LO, "", 1, sizeof(p->reconstruct_method), NULL}};/* 0054 1103 */
  DCM_ELEMENT e39[] = {
    {DCM_NMIRANDOMSCORRECTIONMETHOD, DCM_CS, "", 1, sizeof(p->randoms_method), NULL}};/* 0054 1100 */
  DCM_ELEMENT e40[] = {
    {DCM_ACQRADIOPHARMACEUTICAL, DCM_LO, "", 1, sizeof(p->rad), NULL}};		/* 0018 0031 */
  DCM_ELEMENT e42[] = {
    {DCM_NMIDECAYFACTOR, DCM_LO, "", 1, sizeof(p->decay_factor), NULL}};	/* 0054 1321 */
  DCM_ELEMENT e43[] = {
    {DCM_NMINUMBEROFTIMESLICES, DCM_US, "", 1, sizeof(p->nmi_number_of_time_slices), NULL}}; /* 0054 0101 */
     
  CONDITION cond;
  
  e[0].d.us = &p->rows;				/*  DCM_IMGROWS		*/
  e[1].d.us = &p->columns;			/*  DCM_IMGCOLUMNS	*/
  e[2].d.us = &p->bitsAllocated;		/*  DCM_IMGBITSALLOCATED*/
  e[3].d.us = &p->bitsStored;			/*  DCM_IMGBITSSTORED	*/
  e[4].d.us = &p->highBit;			/*  DCM_IMGHIGHBIT	*/
  e[5].d.us = &p->pixelRepresentation;		/*  DCM_IMGPIXELREPRESENTATION */
  e[6].d.string = p->seriesUID;			/*  DCM_RELSERIESINSTANCEUID */
  e[7].d.string = imagetxt;			/*  DCM_RELIMAGENUMBER	*/
  /**************************************************************************/
  e0[0].d.string = series_num_txt;		/*  DCM_RELSERIESNUMBER	*/
  e0[1].d.string = p->photometricInterpretation;/*  DCM_IMGPHOTOMETRICINTERP */
  /**************************************************************************/
  e1[0].d.string = p->pixelSpacing;	/*  DCM_IMGPIXELSPACING		*/
  e1[1].d.string = p->sliceThickness;	/*  DCM_ACQSLICETHICKNESS	*/
  e1[2].d.string = imagePosition;	/*  DCM_RELIMAGEPOSITIONPATIENT	*/
  e1[3].d.string = imageOrientation;	/*  DCM_RELIMAGEORIENTATIONPATIENT */
  /**************************************************************************/
  e2[0].d.us = 	   &p->planarConfiguration; /*  DCM_IMGPLANARCONFIGURATION*/
  e3[0].d.string =  p->studyDate;	    /*  DCM_IDSTUDYDATE	*/
  e4[0].d.string =  p->seriesDescription;   /*  DCM_IDSERIESDESCR	*/
  e5[0].d.string =  p->protocolName;	    /*  DCM_ACQPROTOCOLNAME	*/
  e6[0].d.string =  p->patientName;	    /*  DCM_PATNAME		*/
  e7[0].d.string =  p->patientID;	    /*  DCM_PATID		*/
  /**************************************************************************/
  e8[0].d.us =  &p->imgLargestPix;	/*  DCM_IMGLARGESTIMAGEPIXELVALUE  */
  e8[1].d.us =  &p->imgSmallestPix;	/*  DCM_IMGSMALLESTIMAGEPIXELVALUE */
  /**************************************************************************/
  e9[0].d.string  = p->sliceSpacing;     	/*  DCM_ACQSLICESPACING	*/
  e10[0].d.string = p->acqSequenceName;		/*  DCM_ACQSEQUENCENAME	*/
  e11[0].d.string = p->idManufacturer;		/*  DCM_IDMANUFACTURER	*/
  e12[0].d.string = p->idModel;			/*  DCM_IDMANUFACTURERMODEL */
  /**************************************************************************/
  e13[0].d.string = p->idModality;		/*  DCM_IDMODALITY	    */
  e13[1].d.string = p->metaTransferUID;         /*  DCM_METATRANSFERSYNTAX  */
  /**************************************************************************/
  e14[0].d.string = p->idInstitution;		/*  DCM_IDINSTITUTIONNAME   */
  /**************************************************************************/

  /* PET */
  e15[0].d.string = petSeriesType;	/*  DCM_NMISERIESTYPE     */
  e16[0].d.string = p->image_correction;/*  DCM_IMGCORRECTEDIMAGE */
  e17[0].d.us =    &p->nmi_image_index;	/*  DCM_NMIIMAGEINDEX     */
  e18[0].d.string = p->counts_source;	/*  DCM_NMICOUNTSSOURCE   */
  e19[0].d.us =    &p->nmi_slice_max;	/*  DCM_NMINUMBEROFSLICES */
  e20[0].d.string = frametxt;		/*  DCM_IMGNUMBEROFFRAMES */
  e21[0].d.string = p->atten_correction;/*  DCM_NMIATTENUATIONCORRECTIONMETHOD */
  e22[0].d.string = p->scatt_correction;/*  DCM_NMISCATTERCORRECTIONMETHOD  */
  e23[0].d.string = p->decay_correction;/*  DCM_NMIDECAYCORRECTION */
  /* e24[0] ;*/	/*  DCM_  */
  /* e25[0] ;*/	/*  DCM_  */
  /* e26[0] ;*/	/*  DCM_  */
  e27[0].d.string = p->series_time;		/*  DCM_IDSERIESTIME  */
  e28[0].d.string = p->acq_start_time;		/*  DCM_IDACQUISITIONTIME */
  e29[0].d.string = frame_reftxt;		/*  DCM_NMIFRAMEREFERENCETIME  */
  e30[0].d.string = frdurtxt;			/*  DCM_ACQACTUALFRAMEDURATION */
  e31[0].d.string = p->rescale_intercept;	/* DCM_IMGRESCALEINTERCEPT    */
  e32[0].d.string = p->rescale_slope;		/*  DCM_IMGRESCALESLOPE	*/
  e33[0].d.string = p->lossy_comp;		/*  DCM_IMGLOSSYIMAGECOMPRESSION */
  e34[0].d.string = p->rad_volume;		/*  DCM_ACQRADIOPHARMVOLUME */
  e35[0].d.string = p->nmi_units;		/*  DCM_NMIUNITS */
  /* e36[0] */
  /* e37[0] */
  e38[0].d.string = p->reconstruct_method;	/*DCM_NMIRECONSTRUCTIONMETHOD  */
  e39[0].d.string = p->randoms_method;		/*  DCM_NMIRANDOMSCORRECTIONMETHOD */
  e40[0].d.string = p->rad;			/*  DCM_ACQRADIOPHARMACEUTICAL */
  e42[0].d.string = p->decay_factor;		/*  DCM_NMIDECAYFACTOR */
  e43[0].d.us =	   &p->nmi_number_of_time_slices;/*  DCM_NMINUMBEROFSLICES */
  
  /**************************************************************************/
  
  cond = DCM_ParseObject(obj, e13, (int)DIM_OF(e13), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  
  /* if((strcmp(p->idModality,"MR") == 0))printf("\nModality MR\n"); */
  /* if((strcmp(p->idModality,"CT") == 0))printf("\nModality CT\n"); */
  if((strcmp(p->idModality,"SR") == 0))printf("\nModality SR\n");
  
  /**************************************************************************/
  cond = DCM_ParseObject(obj, e, (int)DIM_OF(e), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
    printf ("parseParams: Problem Parsing e: rows, columns, bits, seriesUID, pix rep, image number\n");
    (void)COND_PopCondition(TRUE);
    COND_DumpConditions();
    return 1;
  }
  cond = DCM_ParseObject(obj, e0, (int)DIM_OF(e0), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
     printf ("Problem Parsing e0: Series Num, Photo Interp in parseParams\n");
     (void)COND_PopCondition(TRUE);
  }
  
  cond = DCM_ParseObject(obj, e1, (int)DIM_OF(e1), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {
     printf ("Problem Parsing e1: in parseParams\n");
     printf ("Check DCM_IMGPIXELSPACING 0028 0030\n");		/* p->pixelSpacing */
     printf ("Check DCM_ACQSLICETHICKNESS 0018 0050\n");	/* p->sliceThickness */
     printf ("Check DCM_RELIMAGEPOSITIONPATIENT 0020 0032\n");	/* imagePosition */
     printf ("Check DCM_RELIMAGEORIENTATIONPATIENT 0020 0037\n"); /* imageOrientation */
     /*(void)COND_PopCondition(TRUE);*/
     COND_DumpConditions();
     return 1;
  }
  {
    char *c;
    c = strtok(imagePosition, "\\");	p->positionX = atof(c);
    c = strtok(NULL, "\\");		p->positionY = atof(c);
    c = strtok(NULL, "\\");		p->positionZ = atof(c);
  }
  {
    char *o;
    o = strtok(imageOrientation, "\\");	p->orientrowX = strtod(o, NULL);
    o = strtok(NULL, "\\");		p->orientrowY = strtod(o, NULL);
    o = strtok(NULL, "\\");		p->orientrowZ = strtod(o, NULL);
    o = strtok(NULL, "\\");		p->orientcolX = strtod(o, NULL);
    o = strtok(NULL, "\\");		p->orientcolY = strtod(o, NULL);
    o = strtok(NULL, "\\");		p->orientcolZ = strtod(o, NULL);
  }
  if (sliceThickness != NULL && *sliceThickness != '\0'  ) {
	  /* -t option used set sliceThickness */
     strcpy(p->sliceThickness, sliceThickness);
  }
  
  cond = DCM_ParseObject(obj, e2, (int)DIM_OF(e2), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE); 
     p->planarConfiguration = 0;
  }
  cond = DCM_ParseObject(obj, e3, (int)DIM_OF(e3), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e4, (int)DIM_OF(e4), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e5, (int)DIM_OF(e5), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e6, (int)DIM_OF(e6), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e7, (int)DIM_OF(e7), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  
  cond = DCM_ParseObject(obj, e8, (int)DIM_OF(e8), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  
  cond = DCM_ParseObject(obj, e9, (int)DIM_OF(e9), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e10, (int)DIM_OF(e10), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e11, (int)DIM_OF(e11), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e12, (int)DIM_OF(e12), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
                           /* e13 modality is parsed first to check for SR */
  cond = DCM_ParseObject(obj, e14, (int)DIM_OF(e14), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e15, (int)DIM_OF(e15), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e16, (int)DIM_OF(e16), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e17, (int)DIM_OF(e17), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e18, (int)DIM_OF(e18), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e19, (int)DIM_OF(e19), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e20, (int)DIM_OF(e20), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
 /*****************************************************************/
  if(petSeriesType[0] == 'W')strcpy(p->nmi_series_type, "WHOLE BODY");
  if(petSeriesType[0] == 'S')strcpy(p->nmi_series_type, "STATIC");
  if(petSeriesType[0] == 'D')strcpy(p->nmi_series_type, "DYNAMIC");
  if(petSeriesType[0] == 'G')strcpy(p->nmi_series_type, "GATED");
  /*****************************************************************/
  cond = DCM_ParseObject(obj, e21, (int)DIM_OF(e21), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e22, (int)DIM_OF(e22), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e23, (int)DIM_OF(e23), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  
  /*cond = DCM_ParseObject(obj, e24, (int)DIM_OF(e24), NULL, 0, NULL);	*/
  /*if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}		*/
  /*cond = DCM_ParseObject(obj, e25, (int)DIM_OF(e25), NULL, 0, NULL);	*/
  /*if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}		*/


  cond = DCM_ParseObject(obj, e27, (int)DIM_OF(e27), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e28, (int)DIM_OF(e28), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e29, (int)DIM_OF(e29), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e30, (int)DIM_OF(e30), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e31, (int)DIM_OF(e31), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e32, (int)DIM_OF(e32), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e33, (int)DIM_OF(e33), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e34, (int)DIM_OF(e34), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e35, (int)DIM_OF(e35), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  
  /*cond = DCM_ParseObject(obj, e36, (int)DIM_OF(e36), NULL, 0, NULL);	*/
  /*if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}		*/
  /*cond = DCM_ParseObject(obj, e37, (int)DIM_OF(e37), NULL, 0, NULL);	*/
  /*if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}		*/
  
  cond = DCM_ParseObject(obj, e38, (int)DIM_OF(e38), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e39, (int)DIM_OF(e39), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e40, (int)DIM_OF(e40), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e42, (int)DIM_OF(e42), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  cond = DCM_ParseObject(obj, e43, (int)DIM_OF(e43), NULL, 0, NULL);
  if (cond != DCM_NORMAL) {(void)COND_PopCondition(TRUE);}
  
 /*****************************************************************/
 
  if((strcmp(p->idModality,"PT") == 0)){
    p->frame_ref =	   atof(frame_reftxt);	/* NMI Frame Ref Time	     0054 1300 */
    p->frame_duration =	   atof(frdurtxt);	/* ACQ Actual Frame Duration 0018 1242 */
    p->img_number_frames = atoi(frametxt);	/* IMG number of frames	     0028 0008 */
  }
  p->seriesNumber =	   atoi(series_num_txt);	/* REL series number 0020 0011 */
  p->relativeImgNumber =   atoi(imagetxt);		/* REL image number  0020 0013 */
  
  /* Set instanceNumber. For PET use [nmi_image_index 0054 1330 ] */
  if((strcmp(p->idModality,"PT") == 0))p->instanceNumber = p->nmi_image_index;
  
  /* MRI will sort by instanceNumber that is REL image number 0020 0013 */
  if((strcmp(p->idModality,"MR") == 0)){
  	p->instanceNumber = atoi(imagetxt);
  	p->img_number_frames = 1;
  }
  
  /*****************************************************************************/
  if((strcmp(p->idModality,"PT") == 0)){
 
	DCM_ELEMENT e26;
	DCM_SEQUENCE_ITEM *sq;
	e26.tag = DCM_NMIRADIOPHARMINFOSEQ; /* 0054 0016 SQ parsed by DumpElements*/
       
	cond = DCM_GetSequenceList(obj, e26.tag, &e26.d.sq);
	if (cond != DCM_NORMAL)
               printf("cond != DCM_NORMAL GetSequenceList(obj, e26.tag, &e26.d.sq)\n");

	sq = LST_Head(&e26.d.sq);
	if (sq != NULL)
	   (void) LST_Position(&e26.d.sq, sq);

	while (sq != NULL) {
	   (void) DumpElements(&sq->object, 1, p);
	   sq = LST_Next(&e26.d.sq);
	}
	/*****************************************************************************/
	if(verbose == 2 && strcmp(p->idModality,"PT") == 0){
	    printf("Returned From Parsing Header Sequences:\n");
            printf("DCM_IDCODEVALUE: Rad Code %s\n",		  p->id_code_value );
            printf("DCM_IDCODEMEANING: Rad Code Meaning %s\n",	  p->id_code_meaning );
            printf("DCM_ACQRADIOPHARMSTARTTIME: inject_time %s\n",p->inject_time );
            printf("DCM_ACQRADIONUCLIDETOTALDOSE: rad_dose %s\n", p->rad_dose );
            printf("DCM_ACQRADIOPHARMVOLUME: rad_volume %s\n", 	  p->rad_volume );
            printf("DCM_ACQRADIOPHARMACEUTICAL: rad %s\n", 	  p->rad );
	}
  }
  
  return 0;
}

/******************************************************************************/
#define BINARY_CHECK(a, b, str, f) if ((a) != (b)) { \
   fprintf(stderr, " Error: Values do not match (%d %d) for %s in %s\n", (a), (b), (str), (f)); \
   exit(3); \
}

/******************************************************************************/
#define STRING_CHECK(a, b, str, f) if (strcmp((a), (b)) != 0) { \
  fprintf(stderr, " Error: Values do not match (%s %s) for %s in %s\n", (a), (b), (str), (f)); \
  exit(3); \
}

/******************************************************************************/
static CONDITION testSelectionText(DCM_OBJECT** obj,
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

/******************************************************************************/
static CONDITION testRangeText(DCM_OBJECT** obj,
			       DCM_TAG tag, const char* lowerRange,
			       const char* upperRange)
{
  char buf[1024];
  CONDITION cond;
  DCM_ELEMENT e;
  long lowerVal;
  long upperVal;
  long testVal;
  DCM_ELEMENT e1;
  char buf1[1024];

  memset(&e, 0, sizeof(e));
  e.tag = tag;
  e.d.string = buf;
  e.length = sizeof(buf) - 1;

  memset(&e1, 0, sizeof(e1));
  e1.tag = DCM_RELIMAGENUMBER;
  e1.d.string = buf1;
  e1.length = sizeof(buf1) - 1;


  cond = DCM_ParseObject(obj, &e, 1, 0, 0, NULL);
  if (cond != DCM_NORMAL) {
    COND_DumpConditions();
    exit(1);
  }

  cond = DCM_ParseObject(obj, &e1, 1, 0, 0, NULL);
  if (cond != DCM_NORMAL) {
    COND_DumpConditions();
    exit(1);
  }

  testVal = atol(buf);
  lowerVal = atol(lowerRange);
  upperVal = atol(upperRange);

  if (testVal < lowerVal)
    return 1;

  if (testVal > upperVal)
    return 1;

  return 0; /* testVal is within lower and upper values */
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
    
    if (verbose == 2)printf("readParamsForOneImage options = %x\n",options);
    
    cond = DCM_OpenFile(fileName, options, &obj);
    
    if (cond != DCM_NORMAL) {
      printf("Error readParamsForOneImage DCM_OpenFile %s\n",fileName);
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
  
  if (rangeTag != 0) {
    cond = testRangeText(&obj, rangeTag, lowerRange, upperRange);
    if (cond != 0) {
      free(p);
      DCM_CloseObject(&obj);
    }
  }
  
  cond = parseParams(&obj, sliceThickness, p);
  
  if (cond != 0) {
        printf ("Error in readParamsForOneImage call to parseParams %s\n",fileName);
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
    if (verbose == 2)
      printf("readImageParams: %s\n", *argv);
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
    
      if (verbose == 2){
          printf("readImageParams: DCM_ORDERLITTLEENDIAN=%d \n",DCM_ORDERLITTLEENDIAN);
          printf("readImageParams: FILEFORMATMASK=%x PART10FILE=%x ACCEPTVRMISMATCH=%x \n",\
                  DCM_FILEFORMATMASK,DCM_PART10FILE,DCM_ACCEPTVRMISMATCH);
      }
      options = DCM_ORDERLITTLEENDIAN;			/* 0x02 */
      options &= ~DCM_FILEFORMATMASK;			/* 0x80 */
         if (verbose == 2)printf("readImageParams: options = %x\n",options);	/* 0 */
      
      options |= DCM_PART10FILE;			/* 0x80 */
         if (verbose == 2)printf("readImageParams: options = %x\n",options);	/* 80 */
      
      options |= DCM_ACCEPTVRMISMATCH;			/* 0x4000 */
         if (verbose == 2)printf("readImageParams: options = %x\n",options);	/* 4080 */
         if (verbose == 2)printf("readImageParams: NOW DCM_OpenFile *argv, options, &obj\n",options);
      
      cond = DCM_OpenFile(*argv, options, &obj);
      if (cond != DCM_NORMAL) {
        printf("readImageParams: Problem reading DCM_PART10FILE\n");
        
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
    } else if (rangeTag != 0) {
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
static void 
setDimensions(char* base, LST_HEAD** l, char* spacingFlag)
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
	     if (p1.rows == 0)p1 = *p;

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

	     if (verbose == 2)printf("setDimensions:: X %f %f\tY %f %f\tZ %f %f\n", \
	     p1.positionX, p->positionX, p1.positionY, p->positionY, p1.positionZ, p->positionZ );

	     p = LST_Next(l);
	}
	p = LST_Head(l);

	/* dsr->dime.dim[0] number of dimensions for 4dfp = 4 */
	imageDim[0] = p1.columns;		/* xdim (columns) dsr->dime.dim[1] */
	imageDim[1] = p1.rows;			/* ydim (rows) dsr->dime.dim[2] */
	imageDim[2] = LST_Count(l)/volCount;	/* zdim [nodes in LST/# in seq] */
	  					/* imageDim[2] dsr->dime.dim[3] set in setSliceVol */
	imageDim[3] = volCount;			/* volumes: dsr->dime.dim[4] see setSliceVol */

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
		 				/* zpix pixel sizes (mm)voxelSize[2] */
	  					/* -t sets sliceThickness number in parseParams */

	if (strcmp(spacingFlag, "I") == 0) {		/* Use -t <flt> for pixdim[3] */
	     voxelSize[2] = atof(p1.sliceThickness);
	     if (verbose==2)printf("setDimensions:I Input sliceThickness is %f\n",atof(p1.sliceThickness));

	}else{					/* Check for selected tag or use default */

	     if (strcmp(spacingFlag, "T") == 0) {	/* Use -t T (0018 0050) */
	           voxelSize[2] = atof(p1.sliceThickness);
	           if (verbose==2)printf("setDimensions:T Flag: Use p1.sliceThickness for pixdim[3]\n");

	     }else if (p1.sliceSpacing[0] != '\0' || strcmp(spacingFlag, "S") == 0) {
	           voxelSize[2] = atof(p1.sliceSpacing);
	           if (verbose==2)printf("setDimensions:S Flag: SliceSpacing = %f\n",atof(p1.sliceSpacing));

	     }else if (p1.sliceThickness[0] != '\0') {
	          voxelSize[2] = atof(p1.sliceThickness);
	          if (verbose==2)printf("setDimensions:Default: SliceSpacing NULL and -t not used \n");

	     }else{
	          spacingFlag = "C"; /* p1.sliceThickness Is NULL so set Flag to compute */
	          if (verbose==2)printf("setDimensions:Compute: sliceSpacing NULL & sliceThickness NULL & -t Not Used \n");
	     }
	}

	if (verbose){
	      printf("\nsetDimensions: imageDim[0] xdim %d\n",imageDim[0] );
	      printf("setDimensions: imageDim[1] ydim %d\n",imageDim[1] );
	      printf("setDimensions: pixdim[3] %f\n", voxelSize[2]);
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
**	The purpose of this function is to orient an image stack that has
**	been sorted by instance number. The number of volumes and
**	orientation must be known. This function sorts again using
**	the REL image position to sort the image into the Analyze
**	Coordinate System. Change the default ordering with xyzDirection.      
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
  PARAMS *p1;
  LST_HEAD *lp;
  
  int	i,j,img_in_frame;
  
  int	orient =	orientation;
  int	imgCount=	imageDim[2];	/* zdim - images in one frame */
  int	VolCount =	imageDim[3];	/* vols - number of frames */
  
  p = LST_Head(l);
  (void)LST_Position(l, p);
 
  p = LST_Next(l);
   if (VolCount == 1){
                           /********************************************************/
      if (orient == 2){	   /** Transverse default is Z Image Position low to high **/
         if (zDirection) { /** with first image the most ventral **/
             printf("\nTransverse Z-axis ordered high to low\n");
             LST_Sort(l, sizeof(*p), zOrderhl);
         } else {
             LST_Sort(l, sizeof(*p), zOrderlh);
         }
      }                    /********************************************************/
      if (orient == 3){	   /** Coronal default is Y Image Position low to high    **/
         if (yDirection) { /** with first image most posterior **/
             printf("\nCoronal Y-axis ordered high to low\n");
             LST_Sort(l, sizeof(*p), yOrderhl);
         } else {
             LST_Sort(l, sizeof(*p), yOrderlh);
         }
      }                    /********************************************************/
      if (orient == 4){	   /** Sagittal default is X Image Position high to low   **/
         if (xDirection) { /** with first image on the right **/
             printf("\nReordering Sagittal X-axis low to high\n");
             LST_Sort(l, sizeof(*p), xOrderlh);
         } else {
             LST_Sort(l, sizeof(*p), xOrderhl);             
         }
      }
   }
                         /**************************************/
   if (VolCount > 1){	 /** Check orientation of each volume **/

      p = LST_Head(l);
      (void)LST_Position(l, p);
      
      lp = LST_Create();
      p1 = LST_Head(&lp);
           LST_Position(&lp,p1);
      
      if(verbose)printf("orderXYZSlices: ");
      
      for(i=0;i<VolCount;i++){
     
           for(j=0;j<imgCount;j++){
              p1=LST_Pop(l);
                 LST_Push(&lp,p1);
           }
           
           img_in_frame = LST_Count(&lp);
           
           if(verbose==2)printf("orderXYZSlices: Frame=%d images in frame=%d\n",i+1,img_in_frame);
           if(verbose)printf("."); 
           
           p1 = LST_Head(&lp);
                LST_Position(&lp,p1);
           
           if (orient == 2) {  /** Transverse **/
              if (zDirection) {
                 LST_Sort(&lp, sizeof(*p1), zOrderhl);
              } else {
                 LST_Sort(&lp, sizeof(*p1), zOrderlh);
              }
           }
           if (orient == 3) {  /** Coronal **/
              if (yDirection) {
                 LST_Sort(&lp, sizeof(*p1), yOrderhl);
              } else {
                 LST_Sort(&lp, sizeof(*p1), yOrderlh);
              }
           }
           if (orient == 4) {  /** Sagittal **/
              if (xDirection) {
                 LST_Sort(&lp, sizeof(*p1), xOrderlh);
              } else {
                 LST_Sort(&lp, sizeof(*p1), xOrderhl);
              }
           }
           
           for(j=0;j<imgCount;j++){
              p1=LST_Dequeue(&lp);
                 LST_Enqueue(l,p1);
           }                      
      }
   }
  
}

/*****************************************************************************/
static void writeImg(const char* base, LST_HEAD** l, char control)
{ 
  CONDITION	cond;
  CTNBOOLEAN	PET=0;
  DCM_OBJECT*	obj;
  PARAMS	*p;
  FILE		*outfp;
  char		fileName[1024],extension[6]="4dfp";
  int		i, test=0, pixelCount, slicetot, volcount;
  void		*fpixels = NULL;	/* pixel pointer */
  float		*img;
  S16		*imgi;
  double	rescale=1.0;
  
  sprintf(extension,"%s","4dfp");
  
  sprintf(fileName, "%s.%s.%s", base, extension, "img");
  if (!(outfp = fopen (fileName, "wb")))  errw (program, fileName);
   printf("Writing %s\n",fileName);

/*  dsr->hist.orient == 2 orientation == 4 Sagittal
    dsr->hist.orient == 1 orientation == 3 Coronal
    dsr->hist.orient == 0 orientation == 2 Transverse
*/

  p = LST_Head(l);
  (void)LST_Position(l, p);
  while (p != NULL) {
  
     unsigned long options = DCM_ORDERLITTLEENDIAN;
     if (p->part10Flag){
    	options &= ~DCM_FILEFORMATMASK;
	options |= DCM_PART10FILE;
	options |= DCM_ACCEPTVRMISMATCH;
     }
     
     if (verbose == 2){
       printf("writeImg: options = %x\n",options);
       printf("writeImg:  %s\n", p->fileName);
     }
     
     cond = DCM_OpenFile(p->fileName, options, &obj);
     if (cond != DCM_NORMAL) {
       COND_DumpConditions();
       exit(6);
     }
     
     fpixels = extractPixels(&obj, p);
     
     if(p->img_number_frames > 0){
		pixelCount = p->rows * p->columns * p->img_number_frames;
     }else{
		pixelCount = p->rows * p->columns;
     }

     if (!(img = (float *) calloc (pixelCount, sizeof (float)))){
        perror(fileName);
        exit(5);
     }
     
     if (verbose == 2){
        printf("writeImg: Returned From extractPixels\n");
     }
     
     if((strcmp(p->idModality,"PT") == 0)){
        PET = 1;
        rescale = atof(p->rescale_slope);
        if(rescale == 0) rescale = 1;
     }else{
        PET = 0;
     }
     
     float* p1;
     i=0;
     p1 = fpixels;
     while(i < pixelCount) {
	     if(PET)img[i] = *p1 * rescale;
	     i++;	  
	     p1++;
     }

     test++;
     if (verbose == 2)printf("writeImg: test=%d rescale=%f write p->fileName= %s\n",test,rescale,p->fileName);
     
     if (ewrite (img,pixelCount,control,outfp)) errw(program,fileName);

     free(fpixels);
     free(img);

     DCM_CloseObject(&obj);

     p = LST_Next(l);
     
  } /* end of while */
  
  fclose (outfp);
}
/******************************************************************************/
/*  FUNCTION: setFrameMinMax
**	
**  FUNCTIONAL DESCRIPTION:
**	The pixel values are all checked and the max and min for each frame is 
**	determined. 
**  FORMAL PARAMETERS:
**	Two float arrays of minimum and maximum pixel values, image list,
**	sliceNum is the number of images in a single frame.
**  RETURN VALUE:
**	The float arrays are passed by reference.
*/

static void setFrameMinMax(float* Min, float* Max, int sliceNum, LST_HEAD** l)
{
  CONDITION	cond;
  CTNBOOLEAN 	PET=0;
  void *	pixels = NULL;
  PARAMS *	p;
  DCM_OBJECT*	obj;
  int		pixelCount, actual, count, imgCountInFrame, imgCount, imgInFrameMax, frameCount;
  float		rescale=1.0;
  
  float minVal =  65537;
  float maxVal = -65538;
    
  imgCount = imgCountInFrame = frameCount = 1;
  
  imgInFrameMax = sliceNum;

  p = LST_Head(l);
  (void)LST_Position(l, p);
  
  while (p != NULL) {
    
    unsigned long options = DCM_ORDERLITTLEENDIAN;
    if (p->part10Flag){
       options &= ~DCM_FILEFORMATMASK;
       options |= DCM_PART10FILE;
       options |= DCM_ACCEPTVRMISMATCH;
    }

    if (verbose == 2){
      printf("setFrameMinMax: options = %x\n",options);
      printf(" INDEX %d %s\n",p->nmi_image_index, p->fileName);
    }
    
    cond = DCM_OpenFile(p->fileName, options, &obj);
    if (cond != DCM_NORMAL) {
      printf("setFrameMinMax: DCM_OpenFile Failed \n");
      COND_DumpConditions();
      exit(6);
    }
    
    pixels = extractPixels(&obj, p);
    pixelCount = p->rows * p->columns; /* one image at a time*/
        
    if (verbose == 2){
      printf("setFrameMinMax: Returned From extractPixels\n");
    }
    
    float* p1;
    p1 = pixels;
    count = pixelCount;
       
    while(count-- > 0) {
       
	  if (*p1 > maxVal) maxVal = *p1;
	  if (*p1 < minVal) minVal = *p1;
	  
	  p1++;
    }
		
    DCM_CloseObject(&obj);
		
    if(verbose == 2)printf("setFrameMinMax:Frame=%d NMI_INDEX=%d REL_IMG=%d imgCount=%d INSTANCE=%d\n",\
       frameCount,p->nmi_image_index,p->relativeImgNumber,imgCount,p->instanceNumber);
    
    if(imgCountInFrame == imgInFrameMax){
            
            if((strcmp(p->idModality,"PT") == 0)){
               
                rescale = atof(p->rescale_slope);
                if (rescale == 0) rescale = 1;
                
                minVal = minVal * rescale;
                maxVal = maxVal * rescale;
            }
            
            Min[frameCount] = minVal;
            Max[frameCount] = maxVal;
	    
	    if(verbose)printf("setFrameMinMax:Frame=%d Rescaled by %.2f Max[%d]=%.2f Min[%d]=%.2f\n",\
                       frameCount, rescale, frameCount, Max[frameCount], frameCount, Min[frameCount]);

	    minVal = 65537;
	    maxVal = -65538;
	    
	    imgCountInFrame = 0;
	    frameCount++;
    }
    imgCount++;
    imgCountInFrame++;
    
    p = LST_Next(l);
  }

}

/******************************************************************************/
/*  FUNCTION: setSliceVol
**	
**  FUNCTIONAL DESCRIPTION:
**	Number of volumes is determined. 
**	Correct orientation must be known. The volcount is 
**	incrimented if there is a change in image position 
**	AND is the same as the position of the first slice. 
**	
**  FORMAL PARAMETERS:
**	Image list
**  RETURN VALUE:
**	Input structure passed by reference. 
*/

static void setSliceVol(LST_HEAD** l)
{
  
  PARAMS 	*p, *p1;  
  int 		slicenum = 	1;
  int 		volcount = 	1;
  int 		fmriAdjust;
  float 	Xpos, Ypos, Zpos;
  float 	deltaX=0, deltaY=0, deltaZ=0;

  fmriAdjust = LST_Count(l); /* Check for single image volume */
  
  if(verbose == 2)printf("setSliceVol: fmriAdjust= %d\n",fmriAdjust); 
   
  if(fmriAdjust == 1){
      p = LST_Head(l);
     if(p->img_number_frames > 0) imageDim[3] = p->img_number_frames;
     
     if(verbose)printf("setSliceVol: Number of Frames imageDim[3]= %d\n",imageDim[3]); 
     
     return;
  }
  
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
  if (volcount == fmriAdjust) volcount = 1; /* Every image has same location */
  
  imageDim[3] = volcount;
  imageDim[2] = slicenum/volcount;
  
  if (verbose){
  	printf("\nsetSliceVol: Total Images= %d\n",slicenum);
  	printf("setSliceVol: imageDim[2] Images In Frame= %d\n",slicenum/volcount); 
	printf("setSliceVol: imageDim[3] Frames= %d\n",volcount);
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
    
    if(verbose)printf("setSliceThickness: meanSpacing = %f sd = %f\n",meanSpacing,sdSpacing);
    
    voxelSize[2] = meanSpacing;
  } else {
    fprintf(stderr, "Not enough slices to compute slice spacing \n");
    exit(1);
  }
}

/******************************************************************************/
static void writeIFH(const char *base, LST_HEAD** l, char control)
{
  char	fileName[BUF_SIZE];
  int	chars_written;
  char	extension[6]= "4dfp";

  chars_written = snprintf(fileName, BUF_SIZE, "%s.%s.img",
                           base, extension);
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
void writeREC(const char *base, LST_HEAD** l, const char* text, 
              float* imgMin, float* imgMax, int argc, char **argv, 
              CTNBOOLEAN zorder, CTNBOOLEAN xorder, CTNBOOLEAN yorder, 
              char control)
{
  PARAMS 	*p;
  char 		fileName[1024];
  char 		imgName[1024];
  char 		command[1024], *dest;
  int 		imgInFrameCounter, frame, imgCounter = 1;
  int 		k, imgInSeries, len;
  int 		PET, MR, CT, SR;
  float		decay=1, midpoint;

  float		frame_start_time, frame_start_time_check=0, frame_half;
  float		ftemp, frame_start_sec, inject_time_sec, running_tot=0, tot=0;
  int		itemp, hour_sec, hour_test, min_sec, frame_start_flag=0, PMAM_flag=0;

  char 		extension[6]= "4dfp";
  
  float *minVal;
  float *maxVal;

  if (useFolder) {
     mkdir(base, 0777);
  }

    minVal = calloc(imageDim[3]+1, sizeof(float));
    if (minVal == NULL) {
       fprintf(stderr, "writeREC: Could not calloc minVal\n");
       exit(9);
    }

    maxVal = calloc(imageDim[3]+1, sizeof(float));
    if (maxVal == NULL) {
       fprintf(stderr, "writeREC: Could not calloc maxVal\n");
       exit(9);
    }
  
    setFrameMinMax(minVal, maxVal, imageDim[2], l);
  
    *imgMin =  65535;
    *imgMax = -65536;
  
    for(frame=1;frame<=imageDim[3];frame++){
        if(minVal[frame] < *imgMin) *imgMin = minVal[frame];
        if(maxVal[frame] > *imgMax) *imgMax = maxVal[frame];
    }
  
    if(verbose == 2){
       for(frame=1;frame<=imageDim[3];frame++){
           printf("writeREC: frame=%d minVal[frame]=%.2f maxVal[frame]=%.2f \n",
           frame,minVal[frame],maxVal[frame]);
       }
    }
  
  sprintf(fileName, "%s.%s.img", base, extension);
  if(suppress_listing_files) argc = 4;
  printf("Writing %s.rec\n",fileName);

  startrecle (fileName, argc, argv, rcsid, control);

  p = LST_Head(l);
  (void)LST_Position(l, p);
  
  if(outputFileName || verbose){
  	sprintf (command, "Patient_Name:	\t%s\n",p->patientName);printrec(command);
  	sprintf (command, "Patient_ID:	\t%s\n",p->patientID);		printrec(command);
  	sprintf (command, "Institution:	\t%s\n",p->idInstitution);	printrec(command);
  	sprintf (command, "ID Study Date:	\t%s\n",p->studyDate);  printrec(command);
  	sprintf (command, "Volume Division Tag:\t%s\n",text);		printrec(command);
	sprintf (command, "ACQ Sequence Name:\t%s\n",p->acqSequenceName);printrec(command);
	sprintf (command, "ACQ Protocol Name:\t%s\n",p->protocolName);	printrec(command);
	sprintf (command, "ID Sequence Description\t%s\n\n",p->seriesDescription); printrec(command);
  }
  
  if(orientation == 2)      printrec("Orientation:             \tTransverse\n");
  else if(orientation == 3) printrec("Orientation:             \tCoronal\n");
  else if(orientation == 4) printrec("Orientation:             \tSagittal\n");
  else if(orientation > 4){
    sprintf (command,"Orientation:             \t = %d\n",orientation);
    printrec(command);
  }

  sprintf (command, "X dimension (columns): %6d      Width (mm): %10.5f \n",
	  imageDim[0], voxelSize[0]);			printrec(command);
  sprintf (command, "   Y dimension (rows): %6d     Height (mm): %10.5f \n",
	  imageDim[1], voxelSize[1]);			printrec(command);
  sprintf (command, " Z dimension (slices): %6d  Thickness (mm): %10.5f \n",
	  imageDim[2], voxelSize[2]);			printrec(command);
  sprintf (command, "    Number of Volumes: %6d\n ",
	  imageDim[3]);					printrec(command);
  sprintf (command, " Min Value: %.2f\n", *imgMin);	printrec(command);
  sprintf (command, " Max Value: %.2f\n", *imgMax);	printrec(command);
  
  /* Convert injection time into seconds */
  inject_time_sec = atof(p->inject_time);
	
  if (inject_time_sec == 0){ /* correct GE standard violation */
     inject_time_sec = (atof(p->acq_start_time))-60; /* one minute */
	 if(verbose)printf("writeREC: 0018 1072 Injection_Time set to acq_start_time %.2f\n",inject_time_sec);
  }

  ftemp = (inject_time_sec/10000); itemp = ftemp;
  hour_sec = (itemp * 60) * 60;
  hour_test = hour_sec;

  ftemp = ftemp - itemp; ftemp = ftemp * 100; itemp = ftemp;
  min_sec = itemp * 60;

  ftemp = ftemp - itemp; ftemp = ftemp * 100;
  inject_time_sec = hour_sec + min_sec + ftemp;

  running_tot = inject_time_sec;
  
  if(verbose)printf("writeREC: 0018 1072 Injection_Time converted to seconds. %.2f\n",inject_time_sec);
  
  if(outputFileName || verbose){
     sprintf (command, "0008 0020 Scan Date:\t%s\n",p->studyDate);  printrec(command);
     sprintf (command, "0008 0031 Scan Time:\t%s\n",p->series_time);printrec(command);
     sprintf (command, "0008 0060 idModality:\t%s\n",p->idModality);printrec(command);
  }
  /* Check the transfer syntax. Find the last occurance of the number 1, and check for big or little endian syntax */
  if(verbose){	
	len = strlen(p->metaTransferUID); dest=strrchr(p->metaTransferUID, '1');
	printf("0002 0010 metaTransferUID:%s len=%d %s \n",p->metaTransferUID,len,dest);
	if(strcmp(dest,"1.2.2")==0)printf ("BIG ENDIAN Transfer\n");
	if(strcmp(dest,"1.2.1")==0)printf ("LITTLE ENDIAN Transfer\n");
  }
  if(strcmp(p->idModality,"PT")==0){
  
    if(strcmp(p->decay_correction,"NONE")==0){
        decay = 1.0;
    }else{
        decay = atof(p->decay_factor);
    }
    
    if(outputFileName || verbose){
       sprintf (command, "0054 1000 PET_Series_Type:\t%s\n",	p->nmi_series_type);	printrec(command);
       sprintf (command, "0054 0101 Number_of_Time_Slices:\t%d\n",p->nmi_number_of_time_slices);printrec(command);
       sprintf (command, "0054 0081 Slices_in_Series:\t%d\n",	p->nmi_slice_max);	printrec(command);
       sprintf (command, "0054 1002 Source_of_Counts:\t%s\n",	p->counts_source);	printrec(command);
       sprintf (command, "0028 0051 Corrections_Applied:\t%s\n",p->image_correction);	printrec(command);
       sprintf (command, "0054 1101 Attenuation_Method:\t%s\n",	p->atten_correction);	printrec(command);
       sprintf (command, "0054 1100 Randoms_Method:\t%s\n",	p->randoms_method);	printrec(command);
       sprintf (command, "0054 1103 Reconstruct_Method:\t%s\n",	p->reconstruct_method);	printrec(command);
       sprintf (command, "0054 1001 NMI_Units:\t\t%s\n",	p->nmi_units);		printrec(command);
       sprintf (command, "0054 1105 Scatter_Correction:\t%s\n",	p->scatt_correction);	printrec(command);
       sprintf (command, "0054 1102 Decay_Correction:\t%s\n",	p->decay_correction);	printrec(command);
       sprintf (command, "0008 0032 ACQ_Start:\t%s\n",		p->acq_start_time);	printrec(command);
       sprintf (command, "0054 1300 Ref_Time:\t%.3f\n",		p->frame_ref);		printrec(command);
       sprintf (command, "0008 0100 ID_Code_Value:\t%s\n",	p->id_code_value);	printrec(command);
       sprintf (command, "0008 0104 ID_Code_Meaning:\t%s\n",	p->id_code_meaning);	printrec(command);
       sprintf (command, "0018 0031 Radiopharmaceutical:\t%s\n",p->rad);		printrec(command);
       sprintf (command, "0018 1071 Rad_Volume:\t%s\n",		p->rad_volume);		printrec(command);
       sprintf (command, "0018 1072 Injection_Time:\t%s\n",	p->inject_time);	printrec(command);
       sprintf (command, "Converted Injection_Time:\t%.2f\n",	inject_time_sec);	printrec(command);
       sprintf (command, "0018 1074 Rad_Dose_In_Bq:\t%s\n",	p->rad_dose);		printrec(command);
       sprintf (command, "0018 1075 Half_Life_In_Seconds\t%s\n",p->half_life);		printrec(command);
       sprintf (command, "0018 1076 Rad_Positron_Frac:\t%s\n",	p->rad_fraction);	printrec(command);
       sprintf (command, "0002 0010 Meta Transfer Syntax\t%s\n",p->metaTransferUID);	printrec(command);
       sprintf (command, "0028 1052 Rescale_Intercept:\t%s\n\n",p->rescale_intercept);	printrec(command);
    }
    
    /*imgInSeries = p->nmi_slice_max;*/	/* Siemens Dynamic Series is one frame. Phillips Dynamic Series is all frames. */
      imgInSeries = imageDim[2];
    
    imgInFrameCounter = 1;	/* Count images in a single frame */
    imgCounter = 1;		/* Count images */
    frame = 1;			/* Count frames */
    
    /* Print Column Header */
    sprintf (command, "Frame     \tLength(msec)\tMidpoint(sec)\tStart(msec)\t Frame_Min\t Frame_Max\t Decay_Fac\tRescale \n");printrec(command);
    /*sprintf (command, "     \t Length(msec) \tMidpoint(sec)   \tStart(msec)  \tFrame_Min  \tFrame_Max \tDecay_Fac \tRescale \n");printrec(command);*/

    while (p != NULL) {
       if(imageDim[3] > 1){			/* Multivolume */

		 /* if(verbose == 2)printf("writeREC:imgInFrameCounter= %d imgInSeries= %d imgCounter=%d\n", */
		 /* imgInFrameCounter,imgInSeries,imgCounter); */
           
		if(imgInFrameCounter == imgInSeries){ /* Add row of dynamic frame data to rec file. */
                    
		    if(strcmp(p->decay_correction,"NONE")==0){
             		decay = 1.0;
		    }else{
               		decay = atof(p->decay_factor);
		    }
		    /* Convert p->acq_start_time 0008 0032 (24 hour format) into seconds */
		    frame_start_time = atof(p->acq_start_time);

	           ftemp = frame_start_time/10000;
	           itemp = ftemp;
	           hour_sec = (itemp * 60) * 60;

	           /* scans at midnight? */
	           if (hour_sec < hour_test || PMAM_flag == 1){
		      PMAM_flag = 1;
		      hour_sec = 86400 + hour_sec;
		      hour_test = hour_sec;
	           }
		
	           ftemp = ftemp - itemp;
	           ftemp = ftemp * 100;
	           itemp = ftemp;
	           min_sec = itemp * 60;
	     	
	           ftemp = ftemp - itemp;
	           ftemp = ftemp * 100;
	           frame_start_sec = hour_sec + min_sec + ftemp;
	           ftemp = frame_start_sec;
		
	           /* Frame start time is time in seconds from the injection time */
	           frame_start_sec = frame_start_sec - inject_time_sec;
	           midpoint = frame_start_sec + ((p->frame_duration/1000) / 2);

		   if(midpoint_Frame_Ref){ /* use -m option to force use of frame ref and frame duration DICOM tags */
			midpoint = p->frame_ref / 1000; /* Use frame_ref rather than injection time */
			running_tot = running_tot + p->frame_duration;
			frame_start_sec = running_tot / 1000;
		   }

		   if(midpoint_Frame_Dur){ /* use -n option to force use of DICOM frame duration to calculate midpoint and frame start */
			tot = tot + p->frame_duration;	
			frame_start_sec = tot / 1000;

			running_tot = running_tot + frame_half + ((p->frame_duration / 1000)/2);
			frame_half = (p->frame_duration / 1000)/2;

			if(frame == 1){
				frame_start_sec = 0; tot = 0; 
				running_tot = (p->frame_duration / 1000)/2;
				frame_half = running_tot;
			}
			
			midpoint = running_tot;
		   }

		   if(frame_start_time_check == frame_start_sec)frame_start_flag++;
		   frame_start_time_check = frame_start_sec;
			   
		   if(verbose == 2 && frame == 1)printf("writeREC %d Ref=%.3f Dur=%.0f Mid=%.2f frame_start_sec=%.2f frame_start_time=%f p->acq_start_time=%s\n", 
                                frame, p->frame_ref, p->frame_duration, midpoint, frame_start_sec, ftemp, p->acq_start_time);
		   if(verbose == 2)printf("writeREC %d Ref=%.3f %.0f %.2f frame_start_sec=%.2f frame_start_time=%f p->acq_start_time=%s\n", 
                                frame, p->frame_ref, p->frame_duration, midpoint, frame_start_sec, ftemp, p->acq_start_time);
	     
		   sprintf (command, "Frame_%d \t%10.0f \t%10.2f \t%10.0f \t%10.2f \t%10.2f \t%10.6f \t%10s\n",
		      frame,p->frame_duration, midpoint, frame_start_sec*1000, minVal[frame], maxVal[frame], decay, p->rescale_slope);
		   printrec(command);
             
		   imgInFrameCounter = 0;
		   frame++;
		     
		} /* imgInFrameCounter == imgInSeries */
        }else{ /* Add single frame information of the whole body or static scan to rec file. */
             if(imgCounter == 1){
		        sprintf (command, "Frame_%d \t%10.0f \t%10.2f \t%10.0f \t%10.2f \t%10.2f \t%10.6f \t%10s\n",
		        frame,p->frame_duration, midpoint, frame_start_sec*1000, minVal[frame], maxVal[frame], decay, p->rescale_slope);
		        printrec(command);
             }
        }
        imgInFrameCounter++;
        imgCounter++;
        
        p = LST_Next(l);
        
    } /* while (p != NULL) */
  } /* if(PT) */
  
  p = LST_Head(l);
  (void)LST_Position(l, p);
  
  imgCounter = 1;
  
  if(useFolder||outputFileName||xorder||yorder||zorder){  /** Add XYZ and sequence name to REC file **/
      
      if(xorder||yorder||zorder){
         sprintf(command, "Images Reordered - Note REL Image Position 0020 0037 or REL Image Number 0020 0013\n");
         printrec(command);
      }
      if(strcmp(p->idModality,"PT")==0)sprintf(command, "\t\t\t  X  \t\t   Y  \t\t   Z\t\t Img# \t\t Factor\n");
      if(strcmp(p->idModality,"MR")==0)sprintf(command, "\t\t\t  X  \t\t   Y  \t\t   Z\t\t Img#\n");
         printrec(command);

      while (p != NULL) {
         if(useFolder){ 
            sprintf(command,"%s/%s",base,p->fileName);
            rename(p->fileName, command);
         }
	    if(strcmp(p->idModality,"PT")==0){
               sprintf(command, "FileName %d = %s  \t%12.6f \t%12.6f \t%12.6f \t %s \t %d \t %s \t %s\n", 
			imgCounter, p->fileName, p->positionX, p->positionY, p->positionZ, 
			p->acqSequenceName, p->instanceNumber, p->decay_factor, p->rescale_slope );
			printrec(command);
         }
	    if(strcmp(p->idModality,"MR")==0){
		sprintf(command, "FileName %d = %s  \t%12.6f \t%12.6f \t%12.6f \t %s \t %d\n", 
		imgCounter, p->fileName, p->positionX, p->positionY, p->positionZ, 
	        p->acqSequenceName, p->instanceNumber); printrec(command);
         }
         p = LST_Next(l);
         imgCounter++; 
      } 
           
  }else{
  
      if(strcmp(p->idModality,"MR")==0){
	     imgCounter = 1;
		  
	     while (p != NULL) {
		    sprintf(command, "FileName %d = %s\n",imgCounter, p->fileName); printrec(command);
		    p = LST_Next(l);
		    imgCounter++;
	     }
      }	  
  
  }
  
  endrec();
	
  if(frame_start_flag > 0)printf("writeREC: Warning: Acquisition Start Time is the same for %d frames \n", frame_start_flag);
  if(verbose)printf("writeREC: Done writing rec file. \n");
  /*free(minVal);*/
  /*free(maxVal);*/

}

/******************************************************************************/
char getMajorAxisFromDirectionCosine(float x,float y,float z)
{
	char	axis = '\0';
	float	obliquityThresholdCosineValue = 0.8;
	float	absX, absY, absZ;
	char	orientationX, orientationY ,orientationZ;
		
	if(verbose == 2)printf("Directional Cosines: x=%.10f y=%.10f z=%.10f \n",x,y,z); 
		
	orientationX = (x < 0) ? 'R' : 'L';
	orientationY = (y < 0) ? 'A' : 'P';
	orientationZ = (z < 0) ? 'F' : 'H';
		
	absX = fabs(x);
	absY = fabs(y);
	absZ = fabs(z);

	if(verbose)printf("getMajorAxisFromDirectionCosine: absX =%.10f absY =%.10f absZ =%.10f \n",
	absX, absY, absZ);

	/* The tests here really don't need to check the other dimensions,	*/
	/* just the threshold, since the sum of the squares should be == 1.0	*/
	/* but just in case ...							*/
		
	if (absX>obliquityThresholdCosineValue && absX>absY && absX>absZ) {
		axis=orientationX;
	}
	else if (absY>obliquityThresholdCosineValue && absY>absX && absY>absZ) {
		axis=orientationY;
	}
	else if (absZ>obliquityThresholdCosineValue && absZ>absX && absZ>absY) {
		axis=orientationZ;
	}
	
	return axis;
}

/******************************************************************************/
/*  FUNCTION: checkRowColOrientation
**	
**  FUNCTIONAL DESCRIPTION:
**	The directional cosines found in RelImageOrientationPatient 0020 0037
**	are used to determine orientation of the subject in the 3D volume.
**	The change in the directional cosines accross the volume are checked.
**	If there is a change, an error message is printed. 
**	"getMajorAxisFromDirectionCosine" is called to determine the row and col
**	orientation for each axis. 
**  FORMAL PARAMETERS:
**	the pointer to the volume of images
**  RETURN VALUE:
**	An integer representing the Analyze Coordinant System orientation.
*/

int checkRowColOrientation(LST_HEAD** l)
{

  PARAMS 	*p, *p1;
  int 		listcount,orient = 0;
  float		rowX,rowY,rowZ,colX,colY,colZ;
  float		deltaXrow, deltaYrow, deltaZrow, deltaXrowsum, deltaYrowsum, deltaZrowsum;
  float		deltaXcol, deltaYcol, deltaZcol, deltaXcolsum, deltaYcolsum, deltaZcolsum;
  char		colAxis, rowAxis = '\0';
  /* value to return: */
  int		oblique =	1;
  int		transverse =	2;
  int		coronal =	3;
  int		sagittal =	4;

  p = LST_Head(l);
  LST_Position(l, p);
  
  p1 = p;
  p = LST_Next(l);
  
  deltaXrowsum=deltaYrowsum=deltaZrowsum=0;
  deltaXcolsum=deltaYcolsum=deltaZcolsum=0;
  
  while (p != NULL) {
	deltaXrow = labs(p1->orientrowX) - labs(p->orientrowX);
	deltaYrow = labs(p1->orientrowY) - labs(p->orientrowY);
	deltaZrow = labs(p1->orientrowZ) - labs(p->orientrowZ);
    
	deltaXrowsum += labs(deltaXrow);
	deltaYrowsum += labs(deltaYrow);
	deltaZrowsum += labs(deltaZrow);
    
  	deltaXcol = labs(p1->orientcolX) - labs(p->orientcolX);
  	deltaYcol = labs(p1->orientcolY) - labs(p->orientcolY);
  	deltaZcol = labs(p1->orientcolZ) - labs(p->orientcolZ);
    
  	deltaXcolsum += labs(deltaXcol);
 	deltaYcolsum += labs(deltaYcol);
  	deltaZcolsum += labs(deltaZcol);
    
  	p1 = p;
  	p = LST_Next(l);
  }
  
  if( deltaXrowsum != 0 || deltaYrowsum != 0 || deltaZrowsum != 0 ){
      printf("\n*error Xrow=%f Yrow=%f Zrow=%f  ",deltaXrowsum, deltaYrowsum, deltaZrowsum);
  }
  if( deltaXcolsum != 0 || deltaYcolsum != 0 || deltaZcolsum != 0 ){
      printf("\n*error Xcol=%f Ycol=%f Zcol=%f\n",deltaXcolsum, deltaYcolsum, deltaZcolsum);
  }
  
  rowX = p1->orientrowX; rowY = p1->orientrowY; rowZ = p1->orientrowZ;
  colX = p1->orientcolX; colY = p1->orientcolY; colZ = p1->orientcolZ;
    
  if(verbose == 2)printf("\nrowX =%f rowY =%f rowZ =%f",rowX, rowY, rowZ);
  if(verbose == 2)printf("\ncolX =%f colY =%f colZ =%f",colX, colY, colZ);
    
  rowAxis = getMajorAxisFromDirectionCosine(rowX,rowY,rowZ);
  colAxis = getMajorAxisFromDirectionCosine(colX,colY,colZ);
    
    if (rowAxis != '\0' && colAxis != '\0') {
		if      ((rowAxis == 'R' || rowAxis == 'L') && (colAxis == 'A' ||\
		          colAxis == 'P')) orient=transverse;

		else if ((colAxis == 'R' || colAxis == 'L') && (rowAxis == 'A' ||\
		          rowAxis == 'P')) orient=transverse;
		
		else if ((rowAxis == 'R' || rowAxis == 'L') && (colAxis == 'H' ||\
		          colAxis == 'F')) orient=coronal;

		else if ((colAxis == 'R' || colAxis == 'L') && (rowAxis == 'H' ||\
		          rowAxis == 'F')) orient=coronal;
		
		else if ((rowAxis == 'A' || rowAxis == 'P') && (colAxis == 'H' ||\
		          colAxis == 'F')) orient=sagittal;

		else if ((colAxis == 'A' || colAxis == 'P') && (rowAxis == 'H' ||\
		          rowAxis == 'F')) orient=sagittal;
    } else {
		orient=oblique;
    }
    
  if(orient == 4)printf("  Sag  ");
  if(orient == 3)printf("  Cor  ");
  if(orient == 2)printf("  Tra  ");
  if(orient == 1)printf("  Obl  ");
  if(orient == 0)printf("  NA  ");
  
  return orient;
}

/******************************************************************************/
LST_HEAD* fillSelectionList(DCM_TAG tag, LST_HEAD** l)
{
  PARAMS 	  *p;
  DCM_OBJECT* 	  obj;
  char 		  *s;
  SELECTION_ITEM* item;
  LST_HEAD* 	  rtnList;
  CONDITION 	  cond;
  
  if (verbose == 2)
      printf("*IN fillSelectionList: Ready for LST_Create \n");
      
  rtnList = LST_Create();

  if (verbose)
      printf("fillSelectionList: Ready for LST_Head\n");
      
  p = LST_Head(l);
  (void)LST_Position(l, p);
  while (p != NULL) {
    unsigned long options = DCM_ORDERLITTLEENDIAN;
    if (p->part10Flag){
    	options &= ~DCM_FILEFORMATMASK;
        options |= DCM_PART10FILE;
	options |= DCM_ACCEPTVRMISMATCH;
    }

    if (verbose == 2){
      printf("*IN fillSelectionList: fillSel options = %x\n",options);
      printf(" %s\n", p->fileName);
    }
    cond = DCM_OpenFile(p->fileName, options, &obj);
    if (cond != DCM_NORMAL) {
      printf("Unable to fillSelectionList using: %s\n", p->fileName);
      COND_DumpConditions();
      exit(6);
    }
    if (verbose == 2)
      printf("*IN fillSelectionList: Finished opening file \n");
    
    s = DCM_GetString(&obj, tag);
    
    if (verbose == 2)
      printf("*IN fillSelectionList: Finished DCM_GetString \n");
    
    DCM_CloseObject(&obj);
    if (verbose == 2)printf(" *");
    
    item = LST_Head(&rtnList);
    LST_Position(&rtnList, item);
    
    while (item != NULL) {
      if (verbose == 2)printf("*IN fillSelectionList: *item is not NULL\n");
      if (strcmp(s, item->text) == 0) {
	break;
      }
      item = LST_Next(&rtnList);
    }

    if (item == NULL) {
      
      item = malloc(sizeof(*item));
      memset(item, 0, sizeof(*item));
      
      if (verbose == 2)printf("*4 \n");
      strcpy(item->text, s);
      if (verbose == 2)printf("*5 \n");
      
      LST_Enqueue(&rtnList, item);
      /* free(item); WILL Fail */
    }
    
    CTN_FREE(s);
    p = LST_Next(l);
  }
  
  if (verbose == 2)printf("fillSelectionList: Ready To Return rtnList\n");

  return rtnList;
}

/*****************************************************************************
**  FUNCTION: selectImages
**	
**  FUNCTIONAL DESCRIPTION:
**	The purpose of this function is to identify dicom images in the 
**	list "l" using the DCM_TAG "tag" with it's "text". These images 
**	are placed in a new list "localCopy" which is returned. 
*/
LST_HEAD* selectImages(LST_HEAD** l, DCM_TAG tag, const char* text)
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
    if (verbose == 2)printf("selectImages: options = %x\n",options);

    cond = DCM_OpenFile(pOrig->fileName, options, &obj);
    if (cond != DCM_NORMAL) {
      COND_DumpConditions();
      printf("\nMust Exit\n");
      exit(6);
    }
    
    if(verbose == 2)printf("Opened File %s\n",pOrig->fileName);
    
    s = DCM_GetString(&obj, tag);		/** place tag value of image in "s" **/
    DCM_CloseObject(&obj);
    
    if (strcmp(s, text) == 0) {			/** If tag values match **/
        pCopy = malloc(sizeof(*pCopy));
       *pCopy = *pOrig;
       LST_Enqueue(&localCopy, pCopy);	 	/** Add to localCopy list **/
        					
       if(useFileName || useFolder){ 	/** If the -u, or -f option are selected **/
          pOrig = LST_Remove(l, LST_K_AFTER); 	/** shorten the "l" list. **/
          removeFlag = TRUE;			/** Set flag **/
       }
    }
    CTN_FREE(s);
    
    if((useFileName || useFolder) && removeFlag){	/** Flag set AND -u, or -f **/
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
/*  FUNCTION: setBase
**	
**  FUNCTIONAL DESCRIPTION:
**	The p->seriesDescription of the first image is used.
**	The "*" and blank space characters are replaced with a "_".
**	Parens will be changed to hyphens.
**	The series number is appended to create the final output name.
*/
char setBase(LST_HEAD** l, char s[1024])
{

  PARAMS	*p;
  char		*string;
  char		namingTemp[1024];
  int		i=0;

  p = LST_Head(l);
  (void)LST_Position(l, p);
  
  memset(namingTemp, '\0', sizeof(namingTemp));	
  strcpy(namingTemp, p->seriesDescription);
  
  string = &namingTemp[0];
  
  while(*string != '\0') {
  
     if(*string == '*')*string = '^';
     if(*string == ' ')*string = '_';
     if(*string == '(')*string = '_';
     if(*string == ')')*string = '_';
     if(*string == '-')*string = '_';
     
     s[i] = *string;
     string++;
     i++;
  }
  /* Append Rel Series Number to string */
  sprintf(s, "%s_%d",s, p->seriesNumber);
  
}

/******************************************************************************/
/*  FUNCTION: processOneView
**	
**  FUNCTIONAL DESCRIPTION:
**	This function assembles the files for one volume of image data, ifh header, 
**	Analyze hdr data, and rec file information. The list l of DICOM images
**	is sorted, file naming is set, image orientation is determined, and ifh,  
**	hdr, img, and rec files are written. Image lists are managed.
*/

void processOneView(LST_HEAD** l, DCM_TAG tag, const char *text, char* base, int app,
	       CTNBOOLEAN sortOnInstanceNumber, int counter, CTNBOOLEAN xDirect,
	       CTNBOOLEAN yDirect, CTNBOOLEAN zDirect, char* spacingFlag,
	       CTNBOOLEAN computeSliceSpacing, int argc, char **argv, char control)
{
  LST_HEAD*	localCopy;
  PARAMS 	*p;
  char		outfileName[1024];
  char		newbase[1024];
  float         minVal, maxVal;

  		if(verbose)printf("processOneView: Ready to selectImages\n");
  		/* Select images matching the tag and text and create localCopy */
  localCopy = selectImages(l, tag, text); 
  
  		if(verbose)printf("processOneView: Ready to Sort By Instance Number\n");
  orderOnInstance(&localCopy, sortOnInstanceNumber);	/* instance number is set in parseParms */

  		if(verbose)printf("processOneView: Create newbase\n");
  
  if (useFileName) {
      memset(outfileName, '\0', sizeof(outfileName));
      setBase(&localCopy, outfileName);
      		if(verbose)printf("\nprocessOneView: base=%s  outfileName=%s  app=%d\n",base,outfileName,app);
      if(app == 0) {
         sprintf(newbase, "%s_%s", base, outfileName);
      } else {
         sprintf(newbase, "%s_%s_%d", base, outfileName, app);
      }
  }else{
      		if(verbose)printf("processOneView: base=%s  app=%d\n", base, app);
      if(app == 0) {
          sprintf(newbase, "%s", base);
      } else {
          sprintf(newbase, "%s_%d", base, app);
      }
  }
  printf("File%d= %s   %s ", counter, newbase, text);

  		if(verbose)printf("\nprocessOneView: Ready to setDimensions \n");
  setDimensions(newbase, &localCopy, spacingFlag);	/* Initialize the hdr struc */
  		if(verbose)printf("\nprocessOneView: Ready to check orientation \n");
  orientation = checkRowColOrientation(&localCopy);	/* Check directional cosines 0020 0037 */
  
  if (verbose) {
    printf("processOneView:: orientation = %d", orientation);
    printf("processOneView:: imageDims = %d %d %d %d\n",imageDim[0], imageDim[1], imageDim[2], imageDim[3]);
    printf("processOneView:: voxelDims = %d %d %d\n", voxelSize[0], voxelSize[1], voxelSize[2]);
  }

  if (orientation < 5 && orientation > 1){		/* Orientation must be 2, 3, or 4 */
      			if(verbose)printf("\nprocessOneView: Ready to setSliceVol");
	setSliceVol(&localCopy);			/* find number of volumes of sorted data */

	if (computeSliceSpacing)
           setSliceThickness(&localCopy);		/* Force computation of sliceThickness */
          
      			if(verbose)printf("\nprocessOneView: Ready to orderXYZSlices \n");
	orderXYZSlices(&localCopy, xDirect, yDirect, zDirect);
      			if(verbose)printf("\nprocessOneView: Ready to writeImg ");
	writeImg(newbase, &localCopy, control);
			if(verbose)printf("\nprocessOneView: Ready for writeIFH ");
	writeIFH(newbase, &localCopy, control);
      			if(verbose)printf("\nprocessOneView: Ready for writeREC \n");
        writeREC(newbase, &localCopy, text, &minVal, &maxVal, argc, argv, zDirect, xDirect, yDirect, control);
			if(verbose)printf("\nprocessOneView: Ready for writeAnalyzeHeader ");
	writeAnalyzeHeader(newbase,  &minVal, &maxVal);
  } else {			/* Oblique volume */
       outputFileName = TRUE;	/* write a rec file with the file names */
       writeREC(newbase, &localCopy, text, &minVal, &maxVal, argc, argv, zDirect, xDirect, yDirect, control);
  }
  printf("\n");

  /* while ((p = LST_Dequeue(&localCopy)) != NULL) */
  /*    free(p);                                   */
  /* LST_Destroy(&localCopy);                      */
  
}

/******************************************************************************/
/*  FUNCTION: processMultipleViews
**	
**  FUNCTIONAL DESCRIPTION:
**	This function creates an unsorted list of multiple volumes using the 
**	division tag. The list of images in each volume (item) is then
**      sent to processOneView for creation of the Analyze image.
**	Number to append to filenames is incremented.
**  FORMAL PARAMETERS:
**	LST_HEAD** l and the DCM_TAG tag used to select a volume of images.
**	Flags: sortOnInstanceNumber,computeSliceSpacing, spacingFlag,
**	XYZ Direction is used to alter the default orientation if selected.
**	base is the base name used for creating the output filenames
**	Command line arguments are passed for rec file creation.
*/

void processMultipleViews(LST_HEAD** l, DCM_TAG tag, CTNBOOLEAN sortOnInstanceNumber,
		     CTNBOOLEAN xDirection, CTNBOOLEAN yDirection, CTNBOOLEAN zDirection, 
		     CTNBOOLEAN computeSliceSpacing, char* base, char* spacingFlag, 
		     int argc, char **argv, char control)
{
  LST_HEAD* 	  vList = 0;
  SELECTION_ITEM* item;
  SELECTION_ITEM  *p;
  int 		  append = 0;
  char	 	  newBase[1024];
  int		  fileCounter = 1;
  
  if(verbose)printf("processMultipleViews: Run fillSelectionList\n");
  
  vList = fillSelectionList(tag, l); /* Create an unsorted list of multiple division tag items */
  if (vList == 0)
    exit(1);
    
  if(verbose)printf("processMultipleViews: Done With fillSelectionList\n");
  
  item = LST_Head(&vList);
  LST_Position(&vList, item);
  
  while (item != NULL) {
     if(verbose)printf("processMultipleViews: Ready For processOneView \n");
    
     processOneView(l, tag, item->text, base, append,
		   sortOnInstanceNumber, fileCounter,
	           xDirection, yDirection, zDirection,
	           spacingFlag, computeSliceSpacing, 
                   argc, argv, control);
		   
     fileCounter++;
     append++;
     item = LST_Next(&vList);
  }
  
  if(verbose)printf("processMultipleViews: Done in processMultipleViews \n");
  
  /* while ((p = LST_Dequeue(&vList)) != NULL) */
  /*      free(p);                             */
  /* LST_Destroy(&vList);                      */
  
}

/******************************************************************************/
/* main
**
** Purpose:
**	The purpose of this program is to convert dicom images into analyze
**	4D volumes.
**
** Notes:
**	This program will now produce images in the Analyze 4dfp Coordinant
**	System. Interfile ifh file, and record file are also written.
**     
** Algorithm:
**	DICOM PET conversion to Analyze 4dfp orientation:
**	
**	The PT images are separated into Analyze volumes using the series number 
**	0020 0011, or using the -d option.
**
**	Steps to position images in the PET scan are; 1.Order by the instance number
**      that is set in parseParms 2. Find the orientation using 0020 0037 
**	3. Number of volumes, imageDim[2] (dime.dim[3]), is then found using setSliceVol
**	4. Images are placed in Analyze image orientation using relative image
**      position 0020 0032 in the function orderXYZSlices 5. Pixel values are 
**  	rescaled if needed, and the min and max for each frame is determined.
** 
**  PET scan information is formated in the rec file.
**  
**	These items must remain constant through the series:
**      Series Type	0054 1000 (DYNAMIC, WHOLE BODY, STATIC)
**	IMG Corrected	0028 0051 (type of reconstruction)
**	NMI Source	0054 1002 (EMISSION or TRANSMISSION)
**	DYNAMIC Series Type will output a multivolume Analyze image.
**	
*******************************************************************************/

int main(int argc, char **argv)
{
    CTNBOOLEAN xDirection = 		FALSE;
    CTNBOOLEAN yDirection = 		FALSE;
    CTNBOOLEAN zDirection = 		FALSE;
    CTNBOOLEAN computeSliceSpacing = 	FALSE;
    CTNBOOLEAN sortOnInstanceNumber = 	TRUE;
    
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
    DCM_TAG 	divisionTag = DCM_MAKETAG(0x0020,0x0011);  /* Default for PET is REL Series Number */
    char        control =               'O';


    while (--argc > 0 && (*++argv)[0] == '-') {
	
        if(verbose)printf("argc = %d argv = %s\n",argc, *argv);

	switch (*(argv[0] + 1)) {

	case 'a': /* Suppress listing command line files in rec file.  */
		  /* The first three items on the command line will be recorded  */
	    suppress_listing_files = TRUE;
	    break;
	
	case 'b': /* Output base filename follows the -b. The default base name is - analyze -  */
	    argc--; argv++;
	    base = *argv;
	    break;
	case 'c': /* Slice Spacing: Compute using Image Position (0020 0032). (override all other choices) */
	    /*computeSliceSpacing = TRUE;*/
	    printf("computeSliceSpacing = FALSE\n");
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
	case 'g': /* Add image name, XYZ relative position, and image number number to rec file. */
	    outputFileName = TRUE;
	    break;
	case 'i': /*  */
	    break;
	case 'm': /* Use frame ref time for "Midpoint" and "Start" columns in rec file.  */
                  /* Default will calculate with injection time. */
	    midpoint_Frame_Ref = TRUE; midpoint_Frame_Dur = FALSE;
	    break;
	case 'n': /* Use frame duration for "Midpoint" and "Start" columns in rec file.  */
                  /* Default will calculate with injection time. */
	    midpoint_Frame_Dur = TRUE; midpoint_Frame_Ref = FALSE;
	    break;
	case 'p': /* */
	    /* */
	    break;
	case 'q': /*  */
	    break;

	case 't': /* Slice Spacing: Use number or 'S' or 'T' */
	    argc--; argv++; sliceThickness = *argv; spacingFlag = "I";
	    if(strcmp(sliceThickness,"T")==0){spacingFlag = "T"; *sliceThickness = '\0';}
	    if(strcmp(sliceThickness,"S")==0){spacingFlag = "S"; *sliceThickness = '\0';}
	    break;

	case 'u': /**  **/
	    /* */
	    break;
	case 'v': /* verbose  */
	    verbose = 1;
	    break;
	case '2': /* ridiculous verbosity  */
	    verbose = 2;
	    break;
	    
	case 'w': /*  */
	    break;
	case 'x': /*  */
	    break;
	case 'y': /*  */
	    break;
	case 'z': /*  */
	    break;
	     
	case 'X': 
	    printf ("Sagittal- force low to high X order image positions in each frame\n");
	    xDirection = TRUE;
	    break;
	case 'Y': 
	    printf ("Coronal- force high to low Y order image positions in each frame\n");
	    yDirection = TRUE;
	    break;
	case 'Z': 
	    printf ("Transverse- force high to low Z order image positions in each frame\n");
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
    
    if(verbose)printf("main: Ready to readImageParams\n");
    /* Load argv into "l" if selection criteria are met */
    
    readImageParams(argc, argv, selectionTag, selectionText,
		    rangeTag, lowerRange, upperRange,
		    sliceThickness, &l);
    
    if(verbose)printf("main: divisionTag = %d\n",divisionTag);
    /* Send "l" for processing. Divide analyze volumes by divisionTag */
    
    processMultipleViews(&l, divisionTag, sortOnInstanceNumber,
			   xDirection, yDirection, zDirection,
			   computeSliceSpacing, base, spacingFlag,
			   argcf, argvf, control);
    
    /*while ((p = LST_Dequeue(&l)) != NULL) */
    /*   free(p);                           */
    /*LST_Destroy(&l);                      */
    
    if(verbose)printf("main: Done\n");
    
    THR_Shutdown();
    exit(0);
}
