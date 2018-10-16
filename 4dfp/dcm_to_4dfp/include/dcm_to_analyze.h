/*
          Copyright (C) 2001 Washington University

*/
/* Copyright marker.  Copyright will be inserted above.  Do not remove */
/*
**		     Electronic Radiology Laboratory
**		   Mallinckrodt Institute of Radiology
**		Washington University School of Medicine
**
** Module Name(s):	
** Author, Date:	Stephen M. Moore, 14-Oct-2001
** Intent:		
** Last Update:		$Author: mohanar $, $Date: 2006/08/25 20:55:49 $
** Source File:		$RCSfile: dcm_to_analyze.h,v $
** Revision:		$Revision: 1.1 $
** Status:		$State: Exp $
*/


typedef struct {
  void *reserved[2];
  int part10Flag;
  U16 rows;
  U16 columns;
  U16 bitsAllocated;
  U16 bitsStored;
  U16 highBit;
  U16 pixelRepresentation;
  U16 planarConfiguration;
  int instanceNumber;
  float positionX;
  float positionY;
  float positionZ;
  float imageOrientation[6];
  char seriesUID[DICOM_UI_LENGTH + 1];
  char pixelSpacing[DICOM_DS_LENGTH * 2 + 2];
  char sliceThickness[DICOM_DS_LENGTH + 1];
  char photometricInterpretation[DICOM_CS_LENGTH+1];
  char fileName[1024];
  int seriesNumber;
  
  char studyDate[DICOM_DA_LENGTH + 1];
  char seriesDescription[DICOM_PN_LENGTH + 1];
  char protocolName[DICOM_LO_LENGTH + 1];
  char patientName[DICOM_LO_LENGTH + 1];
  char patientID[DICOM_PN_LENGTH + 1];
  
  char acqSequenceName[DICOM_LO_LENGTH + 1];
  char idManufacturer[DICOM_LO_LENGTH + 1];
  char idModel[DICOM_LO_LENGTH + 1];
  char idModality[32 + 1];
  char idInstitution[DICOM_LO_LENGTH + 1];
  
  U16	imgLargestPix;
  U16	imgSmallestPix;
  
  char  sliceSpacing[DICOM_DS_LENGTH + 1];
  
  float  rescaleIntercept;
  float  rescaleSlope;

} PARAMS;

void* extractPixels(DCM_OBJECT** obj, PARAMS *p);






