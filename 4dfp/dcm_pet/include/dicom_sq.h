/*
          Copyright (C) 1993, 1994, RSNA and Washington University

          The software and supporting documentation for the Radiological
          Society of North America (RSNA) 1993, 1994 Digital Imaging and
          Communications in Medicine (DICOM) Demonstration were developed
          at the
                  Electronic Radiology Laboratory
                  Mallinckrodt Institute of Radiology
                  Washington University School of Medicine
                  510 S. Kingshighway Blvd.
                  St. Louis, MO 63110
          as part of the 1993, 1994 DICOM Central Test Node project for, and
          under contract with, the Radiological Society of North America.

          THIS SOFTWARE IS MADE AVAILABLE, AS IS, AND NEITHER RSNA NOR
          WASHINGTON UNIVERSITY MAKE ANY WARRANTY ABOUT THE SOFTWARE, ITS
          PERFORMANCE, ITS MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR
          USE, FREEDOM FROM ANY COMPUTER DISEASES OR ITS CONFORMITY TO ANY
          SPECIFICATION. THE ENTIRE RISK AS TO QUALITY AND PERFORMANCE OF
          THE SOFTWARE IS WITH THE USER.

          Copyright of the software and supporting documentation is
          jointly owned by RSNA and Washington University, and free access
          is hereby granted as a license to use this software, copy this
          software and prepare derivative works based upon this software.
          However, any distribution of this software source code or
          supporting documentation or derivative works (source code and
          supporting documentation) must include the three paragraphs of
          the copyright notice.
          
*/
/* Copyright marker.  Copyright will be inserted above.  Do not remove */

/*
**				DICOM 93
**		     Electronic Radiology Laboratory
**		   Mallinckrodt Institute of Radiology
**		Washington University School of Medicine
**
** Module Name(s):
** Author, Date:	Stephen M. Moore, 21-Jul-93
** Intent:		This is the include file for the DICOM SQ
**			facility.  This facility provides functions for
**			creating and parsing well known sequences.
** Last Update:		$Author: smm $, $Date: 1999/02/06 16:15:15 $
** Source File:		$RCSfile: dicom_sq.h,v $
** Revision:		$Revision: 1.27 $
** Status:		$State: Exp $
*/

#ifndef DICOM_SQ_IS_IN
#define DICOM_SQ_IS_IN 1

#ifdef  __cplusplus
extern "C" {
#endif

typedef enum {
    SQ_K_REFPATIENTSOPINSTANCEUID,
    SQ_K_REFVISITSOPINSTANCEUID,
    SQ_K_REFSTUDYSOPINSTANCEUID,
    SQ_K_REFRESULTSSOPINSTANCEUID,
    SQ_K_REFINTERPRETATIONSOPINSTANCEUID,
    SQ_K_GREYSCALEIMAGEMODULE,
    SQ_K_COLORIMAGEMODULE,
    SQ_K_REFFILMSESSION,
    SQ_K_REFBASICIMAGEBOX,
    SQ_K_REFPRINTJOB,
    SQ_K_ADMITDIAGNOSISCODE,
    SQ_K_DISCHARGEDIAGNOSISCODE,
    SQ_K_PROCEDURECODE,
    SQ_K_REQPROCEDURECODE,
    SQ_K_REQDIAGNOSISCODE,
    SQ_K_REQINTERPRETATIONAPPROVER,
    SQ_K_REFSTUDYCOMPONENT,
    SQ_K_REFSERIESSEQUENCE,
    SQ_K_REFIMAGESEQUENCE,
    SQ_K_REFOVERLAYSEQUENCE,
    SQ_K_REFCURVESEQUENCE,
    SQ_K_CURVEREFOVERLAYSEQUENCE,
    /* all the rest from new NM IOD definition */
    SQ_K_PATIENTORIENTATIONCODESEQUENCE,
    SQ_K_PATIENTORIENTATIONMODIFIERCODESEQUENCE,
    SQ_K_PATIENTGANTRYRELATIONSHIPCODESEQUENCE,
    SQ_K_ANATOMICREGIONSEQUENCE,
    SQ_K_ANATOMICREGIONMODIFIERSEQUENCE,
    SQ_K_PRIMARYANATOMICSTRUCTURESEQUENCE,
    SQ_K_PRIMARYANATOMICSTRUCTUREMODIFIERSEQUENCE,
    SQ_K_ENERGYWINDOWINFOSEQUENCE,
    SQ_K_ENERGYWINDOWRANGESEQUENCE,
    SQ_K_RADIOPHARMACEUTICALINFOSEQUENCE,
    SQ_K_RADIONUCLIDECODESEQUENCE,
    SQ_K_RADIOPHARMACEUTICALROUTECODESEQUENCE,
    SQ_K_CALIBRATIONDATASEQUENCE,
    SQ_K_RADIOPHARMACEUTICALCODESEQUENCE,
    SQ_K_INTERVENTIONDRUGINFOSEQUENCE,
    SQ_K_INTERVENTIONDRUGCODESEQUENCE,
    SQ_K_DETECTORINFOSEQUENCE,
    SQ_K_VIEWCODESEQUENCE,
    SQ_K_VIEWANGULATIONMODIFIERCODESEQUENCE,
    SQ_K_ROTATIONINFOSEQUENCE,
    SQ_K_GATEDINFOSEQUENCE,
    SQ_K_DATAINFOSEQUENCE,
    SQ_K_TIMESLOTINFOSEQUENCE,
    SQ_K_PHASEINFOSEQUENCE,
/*  Start random addtions again */
    SQ_K_MODALITYLUTSEQUENCE,
    SQ_K_VOILUTSEQUENCE,
    SQ_K_REFERENCEDSOPSEQUENCE,			/* Referenced SOP Sequence */
    SQ_K_FAILEDSOPSEQUENCE,			/* Failed SOP Sequence */

/*  Derived from Waveform, Supplement 30 */
    SQ_K_REFPREVIOUSWAVEFORM,			/* 0008 1148 */
    SQ_K_REFSIMULTANEOUSWAVEFORMS,		/* 0008 114A */
    SQ_K_REFSUBSEQUENTWAVEFORM			/* 0008 114C */
}   SQ_TYPE;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
}   SQ_GENERICSTRUCT;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char referencedSOPClassUID[DICOM_UI_LENGTH + 1];
    char referencedSOPInstanceUID[DICOM_UI_LENGTH + 1];
}   SQ_REFPATIENTSOPINSTANCEUID;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char referencedSOPClassUID[DICOM_UI_LENGTH + 1];
    char referencedSOPInstanceUID[DICOM_UI_LENGTH + 1];
}   SQ_REFVISITSOPINSTANCEUID;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char referencedSOPClassUID[DICOM_UI_LENGTH + 1];
    char referencedSOPInstanceUID[DICOM_UI_LENGTH + 1];
}   SQ_REFSTUDYSOPINSTANCEUID;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char referencedSOPClassUID[DICOM_UI_LENGTH + 1];
    char referencedSOPInstanceUID[DICOM_UI_LENGTH + 1];
}   SQ_REFRESULTSSOPINSTANCEUID;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char referencedSOPClassUID[DICOM_UI_LENGTH + 1];
    char referencedSOPInstanceUID[DICOM_UI_LENGTH + 1];
}   SQ_REFINTERPRETATIONSOPINSTANCEUID;

#define	SQ_K_GREYSCALEPIXELASPECTRATIO	0x1
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    unsigned short samplesPerPixel;
    char photometricInterpretation[DICOM_CS_LENGTH + 1];
    unsigned short rows;
    unsigned short columns;
    char pixelAspectRatio[2 * (DICOM_IS_LENGTH + 1)];
    unsigned short bitsAllocated;
    unsigned short bitsStored;
    unsigned short highBit;
    unsigned short pixelRepresentation;
    unsigned char *pixelData;
}   SQ_GREYSCALEIMAGEMODULE;

#define	SQ_K_COLORPIXELASPECTRATIO	0x1
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    unsigned short samplesPerPixel;
    char photometricInterpretation[DICOM_CS_LENGTH + 1];
    unsigned short planarConfiguration;
    unsigned short rows;
    unsigned short columns;
    char pixelAspectRatio[2 * (DICOM_IS_LENGTH + 1)];
    unsigned short bitsAllocated;
    unsigned short bitsStored;
    unsigned short highBit;
    unsigned short pixelRepresentation;
    unsigned char *pixelData;
}   SQ_COLORIMAGEMODULE;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char refSOPClassUID[DICOM_UI_LENGTH + 1];
    char refSOPInstanceUID[DICOM_UI_LENGTH + 1];
}   SQ_REFPRINT;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_REQPROCEDURECODE;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_PROCEDURECODE;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    char refSOPClassUID[DICOM_UI_LENGTH + 1];
    char refSOPInstanceUID[DICOM_UI_LENGTH + 1];
}   SQ_REFSTUDYCOMPONENT;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    char approvalDates[DICOM_DA_LENGTH + 1];
    char approvalTimes[DICOM_TM_LENGTH + 1];
    char physiciansApproving[DICOM_PN_LENGTH + 1];
}   SQ_REQINTERPRETATIONAPPROVER;

#define	SQ_K_REFSERIES_RETRIEVEAETITLE	(1 << 0)
#define	SQ_K_REFSERIES_STORAGEMEDIAFILESETID	(1 << 1)
#define	SQ_K_REFSERIES_STORAGEMEDIAFILESETUID	(1 << 2)
#define SQ_K_REFSERIES_IMAGELIST	(1 << 3)

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char seriesDate[DICOM_DA_LENGTH + 1];
    char seriesTime[DICOM_TM_LENGTH + 1];
    char seriesInstanceUID[DICOM_UI_LENGTH + 1];
    char retrieveAETitle[DICOM_AE_LENGTH + 1];
    char storageMediaFileSetID[DICOM_SH_LENGTH + 1];
    char storageMediaFileSetUID[DICOM_UI_LENGTH + 1];
    LST_HEAD *imageList;
}   SQ_REFSERIESSEQUENCE;

#define	SQ_K_REFIMAGE_RETRIEVEAETITLE	(1 << 0)
#define	SQ_K_REFIMAGE_STORAGEMEDIAFILESETID	(1 << 1)
#define	SQ_K_REFIMAGE_STORAGEMEDIAFILESETUID	(1 << 2)

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char refSOPClassUID[DICOM_UI_LENGTH + 1];
    char refSOPInstanceUID[DICOM_UI_LENGTH + 1];
    char retrieveAETitle[DICOM_AE_LENGTH + 1];
    char storageMediaFileSetID[DICOM_SH_LENGTH + 1];
    char storageMediaFileSetUID[DICOM_UI_LENGTH + 1];
}   SQ_REFIMAGESEQUENCE;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char refSOPClassUID[DICOM_UI_LENGTH + 1];
    char refSOPInstanceUID[DICOM_UI_LENGTH + 1];
}   SQ_REFOVERLAYSEQUENCE;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char refSOPClassUID[DICOM_UI_LENGTH + 1];
    char refSOPInstanceUID[DICOM_UI_LENGTH + 1];
}   SQ_REFCURVESEQUENCE;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char refSOPClassUID[DICOM_UI_LENGTH + 1];
    char refSOPInstanceUID[DICOM_UI_LENGTH + 1];
    unsigned short refOverlayGroup;
}   SQ_CURVEREFOVERLAYSEQUENCE;

#define	SQ_K_VOILUT_LUTEXPLANATION	(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    struct {
	unsigned short lutEntries;
	unsigned short firstPixelValue;
	unsigned short lutDataBits;
    }   lutDescriptor;
    /* unsigned short lutDescriptor[3]; */
    char lutExplanation[DICOM_LO_LENGTH + 1];
    unsigned short *lutData;
}   SQ_VOILUTSequence;

/* all of the rest from NM IOD */
#define SQ_K_PATIENTORIENTATIONCODESEQ_CODEMEANING	(1 << 0)
#define SQ_K_PATIENTORIENTATIONCODESEQ_MODIFIERCODESEQ 	(1 << 1)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
    LST_HEAD *patientOrientationModifierCodeSequence;
}   SQ_PATIENTORIENTATIONCODESEQUENCE;

#define SQ_K_PATIENTORIENTATIONMODIFIERCODESEQ_CODEMEANING	(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_PATIENTORIENTATIONMODIFIERCODESEQUENCE;

#define SQ_K_PATIENTGANTRYRELATIONSHIPCODESEQ_CODEMEANING	(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_PATIENTGANTRYRELATIONSHIPCODESEQUENCE;

#define SQ_K_ANATOMICREGIONSEQ_CODEMEANING		(1 << 0)
#define SQ_K_ANATOMICREGIONSEQ_MODIFIERSEQ	(1 << 1)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
    LST_HEAD *anatomicRegionModifierSequence;
}   SQ_ANATOMICREGIONSEQUENCE;

#define SQ_K_ANATOMICREGIONMODIFIERSEQ_CODEMEANING		(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_ANATOMICREGIONMODIFIERSEQUENCE;

#define SQ_K_PRIMARYANATOMICSTRUCTURESEQ_CODEMEANING		(1 << 0)
#define SQ_K_PRIMARYANATOMICSTRUCTURESEQ_MODIFIERSEQ 		(1 << 1)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
    LST_HEAD *primaryAnatomicStructureModifierSequence;
}   SQ_PRIMARYANATOMICSTRUCTURESEQUENCE;

#define SQ_K_PRIMARYANATOMICSTRUCTUREMODIFIERSEQ_CODEMEANING	(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_PRIMARYANATOMICSTRUCTUREMODIFIERSEQUENCE;

#define SQ_K_ENERGYWINDOWINFOSEQ_ENERGYWINDOWNAME	(1 << 0)
#define SQ_K_ENERGYWINDOWINFOSEQ_RANGESEQ		(1 << 1)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char energyWindowName[DICOM_SH_LENGTH + 1];
    LST_HEAD *energyWindowRangeSequence;
}   SQ_ENERGYWINDOWINFOSEQUENCE;

#define SQ_K_ENERGYWINDOWRANGESEQ_WINDOWLOWERLIMIT	(1 << 0)
#define SQ_K_ENERGYWINDOWRANGESEQ_WINDOWUPPERLIMIT	(1 << 1)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char energyWindowLowerLimit[DICOM_DS_LENGTH + 1];
    char energyWindowUpperLimit[DICOM_DS_LENGTH + 1];
}   SQ_ENERGYWINDOWRANGESEQUENCE;

#define	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMROUTE	(1 << 0)
#define	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMROUTECODESEQ	(1 << 1)
#define	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMVOLUME	(1 << 2)
#define	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMSTARTTIME	(1 << 3)
#define	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMSTOPTIME	(1 << 4)
#define	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIONUCLIDETOTALDOSE	(1 << 5)
#define	SQ_K_RADIOPHARMACEUTICALINFOSEQ_CALIBRATIONDATASEQ	(1 << 6)
#define	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMACEUTICAL	(1 << 7)
#define	SQ_K_RADIOPHARMACEUTICALINFOSEQ_CODESEQ			(1 << 8)
#define	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMNUCLIDECODESEQ  (1 << 9)

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    LST_HEAD *radioNuclideCodeSequence;
    char radioPharmaceuticalRoute[DICOM_LO_LENGTH + 1];
    LST_HEAD *radioPharmaceuticalRouteCodeSequence;
    char radioPharmaceuticalVolume[DICOM_DS_LENGTH + 1];
    char radioPharmaceuticalStartTime[DICOM_TM_LENGTH + 1];
    char radioPharmaceuticalStopTime[DICOM_TM_LENGTH + 1];
    char radioNuclideTotalDose[DICOM_DS_LENGTH + 1];
    LST_HEAD *calibrationDataSequence;
    char radioPharmaceutical[DICOM_LO_LENGTH + 1];
    LST_HEAD *radioPharmaceuticalCodeSequence;
}   SQ_RADIOPHARMACEUTICALINFOSEQUENCE;

#define	SQ_K_RADIONUCLIDECODESEQ_CODEMEANING	(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_RADIONUCLIDECODESEQUENCE;

#define	SQ_K_RADIOPHARMACEUTICALROUTECODESEQ_CODEMEANING	(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_RADIOPHARMACEUTICALROUTECODESEQUENCE;

#define SQ_K_CALIBRATIONDATASEQ_SYRINGECOUNTS		(1 << 0)
#define SQ_K_CALIBRATIONDATASEQ_RESIDUALSYRINGECOUNTS	(1 << 1)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    unsigned short energyWindowNumber;
    char syringeCounts[DICOM_IS_LENGTH + 1];
    char residualSyringeCounts[DICOM_IS_LENGTH + 1];
}   SQ_CALIBRATIONDATASEQUENCE;

#define	SQ_K_RADIOPHARMACEUTICALCODESEQ_CODEMEANING	(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_RADIOPHARMACEUTICALCODESEQUENCE;

#define SQ_K_INTERVENTIONDRUGINFOSEQ_DRUGNAME		(1 << 0)
#define SQ_K_INTERVENTIONDRUGINFOSEQ_CODESEQ		(1 << 1)
#define SQ_K_INTERVENTIONDRUGINFOSEQ_DRUGSTARTTIME	(1 << 2)
#define SQ_K_INTERVENTIONDRUGINFOSEQ_DRUGSTOPTIME	(1 << 3)
#define SQ_K_INTERVENTIONDRUGINFOSEQ_DRUGDOSE		(1 << 4)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char interventionDrugName[DICOM_LO_LENGTH + 1];
    LST_HEAD *interventionDrugCodeSequence;
    char interventionDrugStartTime[DICOM_TM_LENGTH + 1];
    char interventionDrugStopTime[DICOM_TM_LENGTH + 1];
    char interventionDrugDose[DICOM_DS_LENGTH + 1];
}   SQ_INTERVENTIONDRUGINFOSEQUENCE;

#define	SQ_K_INTERVENTIONDRUGCODESEQ_CODEMEANING	(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_INTERVENTIONDRUGCODESEQUENCE;

#define SQ_K_DETECTORINFOSEQ_COLLIMATORGRIDNAME	(1 << 0)
#define SQ_K_DETECTORINFOSEQ_FIELDOFVIEWSHAPE	(1 << 1)
#define SQ_K_DETECTORINFOSEQ_FIELDOFVIEWDIMENSIONS	(1 << 2)
#define SQ_K_DETECTORINFOSEQ_FOCALDISTANCE	(1 << 3)
#define SQ_K_DETECTORINFOSEQ_XFOCUSCENTER	(1 << 4)
#define SQ_K_DETECTORINFOSEQ_YFOCUSCENTER	(1 << 5)
#define SQ_K_DETECTORINFOSEQ_ZOOMCENTER		(1 << 6)
#define SQ_K_DETECTORINFOSEQ_ZOOMFACTOR		(1 << 7)
#define SQ_K_DETECTORINFOSEQ_CENTERFORROTATIONOFFSET		(1 << 8)
#define SQ_K_DETECTORINFOSEQ_GANTRYDETECTORTILT		(1 << 9)
#define SQ_K_DETECTORINFOSEQ_DISTANCESRCTODETECTOR		(1 << 10)
#define SQ_K_DETECTORINFOSEQ_STARTANGLE		(1 << 11)
#define SQ_K_DETECTORINFOSEQ_RADIALPOSITION		(1 << 12)
#define SQ_K_DETECTORINFOSEQ_VIEWCODESEQUENCE		(1 << 13)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char collimatorGridName[DICOM_SH_LENGTH + 1];
    char collimatorType[DICOM_CS_LENGTH + 1];
    char fieldOfViewShape[DICOM_CS_LENGTH + 1];
    char fieldOfViewDimensions[DICOM_IS_LENGTH + 1];
    char focalDistance[DICOM_IS_LENGTH + 1];
    char xFocusCenter[DICOM_DS_LENGTH + 1];
    char yFocusCenter[DICOM_DS_LENGTH + 1];
    char zoomCenter[DICOM_DS_LENGTH + 1];
    char zoomFactor[DICOM_DS_LENGTH + 1];
    char centerForRotationOffset[DICOM_DS_LENGTH + 1];
    char gantryDetectorTilt[DICOM_DS_LENGTH + 1];
    char distanceSrcToDetector[DICOM_DS_LENGTH + 1];
    char startAngle[DICOM_DS_LENGTH + 1];
    char radialPosition[DICOM_DS_LENGTH + 1];
    char imageOrientationPatient[6*(DICOM_DS_LENGTH + 1)];
    char imagePositionPatient[3*(DICOM_DS_LENGTH + 1)];
    LST_HEAD *viewCodeSequence;
}   SQ_DETECTORINFOSEQUENCE;

#define	SQ_K_VIEWCODESEQ_CODEMEANING	(1 << 0)
#define	SQ_K_VIEWCODESEQ_ANGULATIONMODIFIERCODESEQ		(1 << 1)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
    LST_HEAD *viewAngulationModifierCodeSequence;
}   SQ_VIEWCODESEQUENCE;

#define	SQ_K_VIEWANGULATIONMODIFIERCODESEQ_CODEMEANING	(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char codeValue[DICOM_SH_LENGTH + 1];
    char codingSchemeDesignator[DICOM_SH_LENGTH + 1];
    char codeMeaning[DICOM_LO_LENGTH + 1];
}   SQ_VIEWANGULATIONMODIFIERCODESEQUENCE;

#define SQ_K_ROTATIONINFOSEQ_RADIALPOSITION	(1 << 0)
#define SQ_K_ROTATIONINFOSEQ_DISTANCESRCTODETECTOR	(1 << 1)
#define SQ_K_ROTATIONINFOSEQ_TABLETRAVERSE	(1 << 2)
#define SQ_K_ROTATIONINFOSEQ_TABLEHEIGHT	(1 << 3)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char startAngle[DICOM_DS_LENGTH + 1];
    char angularStep[DICOM_DS_LENGTH + 1];
    char rotationDirection[DICOM_CS_LENGTH + 1];
    char scanArc[DICOM_DS_LENGTH + 1];
    char actualFrameDuration[DICOM_IS_LENGTH + 1];
    char radialPosition[DICOM_DS_LENGTH + 1];
    char distanceSrcToDetector[DICOM_DS_LENGTH + 1];
    unsigned short numberOfFramesInRotation;
    char tableTraverse[DICOM_DS_LENGTH + 1];
    char tableHeight[DICOM_DS_LENGTH + 1];
}   SQ_ROTATIONINFOSEQUENCE;

#define SQ_K_GATEDINFOSEQ_TRIGGERTIME	(1 << 0)
#define SQ_K_GATEDINFOSEQ_FRAMINGTYPE	(1 << 1)
#define SQ_K_GATEDINFOSEQ_DATAINFOSEQ	(1 << 2)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char triggerTime[DICOM_DS_LENGTH + 1];
    char framingType[DICOM_LO_LENGTH + 1];
    LST_HEAD *dataInfoSequence;
}   SQ_GATEDINFOSEQUENCE;

#define SQ_K_DATAINFOSEQ_NOMINALINTERVAL	(1 << 0)
#define SQ_K_DATAINFOSEQ_LOWRRVALUE		(1 << 1)
#define SQ_K_DATAINFOSEQ_HIGHRRVALUE		(1 << 2)
#define SQ_K_DATAINFOSEQ_INTERVALSACQUIRED	(1 << 3)
#define SQ_K_DATAINFOSEQ_INTERVALSREJECTED	(1 << 4)
#define SQ_K_DATAINFOSEQ_TIMESLOTINFOSEQ	(1 << 5)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char frameTime[DICOM_DS_LENGTH + 1];
    char nominalInterval[DICOM_IS_LENGTH + 1];
    char lowRRValue[DICOM_IS_LENGTH + 1];
    char highRRValue[DICOM_IS_LENGTH + 1];
    char intervalsAcquired[DICOM_IS_LENGTH + 1];
    char intervalsRejected[DICOM_IS_LENGTH + 1];
    LST_HEAD *timeSlotInfoSequence;
}   SQ_DATAINFOSEQUENCE;

#define SQ_K_TIMESLOTINFOSEQ_TIMESLOTTIME		(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char timeSlotTime[DICOM_DS_LENGTH + 1];
}   SQ_TIMESLOTINFOSEQUENCE;

#define SQ_K_PHASEINFOSEQ_TRIGGERVECTOR		(1 << 0)
#define SQ_K_PHASEINFOSEQ_NUMBEROFTRIGGERSINPHASE	(1 << 1)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char phaseDelay[DICOM_IS_LENGTH + 1];
    char actualFrameDuration[DICOM_IS_LENGTH + 1];
    char pauseBetweenFrames[DICOM_IS_LENGTH + 1];
    unsigned short numberOfFramesInPhase;
    char triggerVector[DICOM_IS_LENGTH + 1];
    unsigned short numberOfTriggersInPhase;
}   SQ_PHASEINFOSEQUENCE;

#define SQ_K_MODALITYLUT_LUTEXPLANATION	(1 << 0)
typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    struct {
	unsigned short lutEntries;
	unsigned short firstPixelValue;
	unsigned short lutDataBits;
    }   lutDescriptor;
    char lutExplanation[DICOM_LO_LENGTH + 1];
    char modalityLUTType[DICOM_LO_LENGTH + 1];
    unsigned short *lutData;
}   SQ_MODALITYLUTSEQUENCE;

#define SQ_K_REFSOP_STORAGEMEDIAID  (1 << 0)
#define SQ_K_REFSOP_STORAGEMEDIAUID (1 << 1)
typedef struct {
    void* reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char referencedSOPClassUID[DICOM_UI_LENGTH+1];
    char referencedSOPInstanceUID[DICOM_UI_LENGTH+1];
    char storageMediaID[DICOM_SH_LENGTH+1];
    char storageMediaUID[DICOM_UI_LENGTH+1];
} SQ_REFERENCEDSOPSEQUENCE;

typedef struct {
    void* reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char referencedSOPClassUID[DICOM_UI_LENGTH+1];
    char referencedSOPInstanceUID[DICOM_UI_LENGTH+1];
    U16 failureReason;
} SQ_FAILEDSOPSEQUENCE;

typedef struct {
    void *reserved[2];
    SQ_TYPE type;
    long conditionalFields;
    char referencedSOPClassUID[DICOM_UI_LENGTH + 1];
    char referencedSOPInstanceUID[DICOM_UI_LENGTH + 1];
}   SQ_GENERICREFERENCEDITEM;

CONDITION
SQ_BuildSequence(LST_HEAD ** list, SQ_TYPE type, DCM_ELEMENT ** element);
CONDITION
SQ_ParseSequence(DCM_ELEMENT * element, SQ_TYPE * type, LST_HEAD ** list);
CONDITION
SQ_ObjectToSequence(DCM_TAG tag, DCM_OBJECT ** object, DCM_ELEMENT ** element);
CONDITION
SQ_SequenceToObject(LST_HEAD ** list, DCM_OBJECT ** object);
CONDITION
SQ_ReleaseList(LST_HEAD ** list);
CONDITION
SQ_Release(void **sequenceStructure);

char *SQ_Message(CONDITION cond);

#define	SQ_NORMAL			FORM_COND(FAC_SRV, SEV_SUCC, 1)
#define	SQ_MALLOCFAILURE		FORM_COND(FAC_SRV, SEV_ERROR, 2)
#define	SQ_LISTCREATEFAILURE		FORM_COND(FAC_SRV, SEV_ERROR, 3)
#define	SQ_OBJECTCREATEFAILED		FORM_COND(FAC_SRV, SEV_ERROR, 4)
#define	SQ_MODIFICATIONFAILURE		FORM_COND(FAC_SRV, SEV_ERROR, 5)
#define SQ_LISTFAILURE			FORM_COND(FAC_SRV, SEV_ERROR, 6)
#define SQ_NULLLIST			FORM_COND(FAC_SRV, SEV_ERROR, 7)
#define SQ_EMPTYLIST			FORM_COND(FAC_SRV, SEV_ERROR, 8)
#define SQ_OBJECTACCESSFAILED		FORM_COND(FAC_SRV, SEV_ERROR, 9)
#define SQ_INTERNALSEQBUILDFAILED	FORM_COND(FAC_SRV, SEV_ERROR, 10)
#define SQ_INTERNALSEQPARSEFAILED	FORM_COND(FAC_SRV, SEV_ERROR, 11)
#define SQ_BADSEQUENCETYPEELEMENT	FORM_COND(FAC_SRV, SEV_ERROR, 12)

#ifdef  __cplusplus
}
#endif

#endif
