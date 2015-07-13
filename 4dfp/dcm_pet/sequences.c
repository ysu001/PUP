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
** Module Name(s):	SQ_BuildSequence
**			SQ_ParseSequence
**			SQ_ObjectToSequence
**			SQ_SequenceToObject
** Author, Date:	Stephen M. Moore, 21-Jul-93
** Intent:
** Last Update:		$Author: jon $, $Date: 2008/04/03 13:47:15 $
** Source File:		$RCSfile: sequences.c,v $
** Revision:		$Revision: 1.1 $
** Status:		$State: Exp $
*/

static char rcsid[] = "$Revision: 1.1 $ $RCSfile: sequences.c,v $";

#include "include/ctn_os.h"

#if 0
#include <stdio.h>
#include <string.h>
#include <errno.h>
#ifndef	MACOS
#include <stdlib.h>
#endif
#include <stdarg.h>
#include <sys/types.h>
#endif

#include "dicom.h"
#include "condition.h"
#include "lst.h"
#include "dulprotocol.h"
#include "dicom_objects.h"
#include "dicom_sq.h"


static SQ_REFPATIENTSOPINSTANCEUID refPatientInstanceUID;
static DCM_ELEMENT refPatientInstanceUIDElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(refPatientInstanceUID.referencedSOPClassUID),
    &refPatientInstanceUID.referencedSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(refPatientInstanceUID.referencedSOPInstanceUID),
    &refPatientInstanceUID.referencedSOPInstanceUID[0]},
};

static SQ_REFVISITSOPINSTANCEUID refVisitInstanceUID;
static DCM_ELEMENT refVisitInstanceUIDElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(refVisitInstanceUID.referencedSOPClassUID),
    &refVisitInstanceUID.referencedSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(refVisitInstanceUID.referencedSOPInstanceUID),
    &refVisitInstanceUID.referencedSOPInstanceUID[0]},
};

static SQ_REFSTUDYSOPINSTANCEUID refStudyInstanceUID;
static DCM_ELEMENT refStudyInstanceUIDElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(refStudyInstanceUID.referencedSOPClassUID),
    &refStudyInstanceUID.referencedSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(refStudyInstanceUID.referencedSOPInstanceUID),
    &refStudyInstanceUID.referencedSOPInstanceUID[0]},
};

static SQ_REFRESULTSSOPINSTANCEUID refResultsInstanceUID;
static DCM_ELEMENT refResultsInstanceUIDElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(refResultsInstanceUID.referencedSOPClassUID),
    &refResultsInstanceUID.referencedSOPClassUID[0]},
    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(refResultsInstanceUID.referencedSOPInstanceUID),
    &refResultsInstanceUID.referencedSOPInstanceUID[0]},
};

static SQ_REFINTERPRETATIONSOPINSTANCEUID refInterpretationInstanceUID;
static DCM_ELEMENT refInterpretationInstanceUIDElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(refInterpretationInstanceUID.referencedSOPClassUID),
    &refInterpretationInstanceUID.referencedSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(refInterpretationInstanceUID.referencedSOPInstanceUID),
    &refInterpretationInstanceUID.referencedSOPInstanceUID[0]},
};


static SQ_GREYSCALEIMAGEMODULE greyScaleImageModule;
static DCM_ELEMENT greyScaleImageModuleR[] = {
    {DCM_IMGSAMPLESPERPIXEL, DCM_US, "", 1,
	sizeof(greyScaleImageModule.samplesPerPixel),
    (void *) &greyScaleImageModule.samplesPerPixel},

    {DCM_IMGPHOTOMETRICINTERP, DCM_CS, "", 1,
	sizeof(greyScaleImageModule.photometricInterpretation),
    (void *) greyScaleImageModule.photometricInterpretation},
    {DCM_IMGROWS, DCM_US, "", 1, sizeof(greyScaleImageModule.rows),
    (void *) &greyScaleImageModule.rows},

    {DCM_IMGCOLUMNS, DCM_US, "", 1, sizeof(greyScaleImageModule.columns),
    (void *) &greyScaleImageModule.columns},

    {DCM_IMGPIXELASPECTRATIO, DCM_IS, "", 1,
	sizeof(greyScaleImageModule.pixelAspectRatio),
    (void *) greyScaleImageModule.pixelAspectRatio},

    {DCM_IMGBITSALLOCATED, DCM_US, "", 1, sizeof(greyScaleImageModule.bitsAllocated),
    (void *) &greyScaleImageModule.bitsAllocated},

    {DCM_IMGBITSSTORED, DCM_US, "", 1, sizeof(greyScaleImageModule.bitsStored),
    (void *) &greyScaleImageModule.bitsStored},

    {DCM_IMGHIGHBIT, DCM_US, "", 1, sizeof(greyScaleImageModule.highBit),
    (void *) &greyScaleImageModule.highBit},

    {DCM_IMGPIXELREPRESENTATION, DCM_US, "", 1,
	sizeof(greyScaleImageModule.pixelRepresentation),
    (void *) &greyScaleImageModule.pixelRepresentation},
};

static SQ_COLORIMAGEMODULE colorImageModule;
static DCM_ELEMENT ColorImageModuleR[] = {
    {DCM_IMGSAMPLESPERPIXEL, DCM_US, "", 1,
	sizeof(colorImageModule.samplesPerPixel),
    (void *) &colorImageModule.samplesPerPixel},

    {DCM_IMGPHOTOMETRICINTERP, DCM_CS, "", 1,
	sizeof(colorImageModule.photometricInterpretation),
    (void *) colorImageModule.photometricInterpretation},

    {DCM_IMGPLANARCONFIGURATION, DCM_US, "", 1,
	sizeof(colorImageModule.planarConfiguration),
    (void *) &colorImageModule.planarConfiguration},

    {DCM_IMGROWS, DCM_US, "", 1, sizeof(colorImageModule.rows),
    (void *) &colorImageModule.rows},
    {DCM_IMGCOLUMNS, DCM_US, "", 1, sizeof(colorImageModule.columns),
    (void *) &colorImageModule.columns},

    {DCM_IMGPIXELASPECTRATIO, DCM_IS, "", 1,
	sizeof(colorImageModule.pixelAspectRatio),
    (void *) colorImageModule.pixelAspectRatio},

    {DCM_IMGBITSALLOCATED, DCM_US, "", 1, sizeof(colorImageModule.bitsAllocated),
    (void *) &colorImageModule.bitsAllocated},

    {DCM_IMGBITSSTORED, DCM_US, "", 1, sizeof(colorImageModule.bitsStored),
    (void *) &colorImageModule.bitsStored},

    {DCM_IMGHIGHBIT, DCM_US, "", 1, sizeof(colorImageModule.highBit),
    (void *) &colorImageModule.highBit},

    {DCM_IMGPIXELREPRESENTATION, DCM_US, "", 1,
	sizeof(colorImageModule.pixelRepresentation),
    (void *) &colorImageModule.pixelRepresentation},
};


static SQ_REFPRINT refPrint;
static DCM_ELEMENT refPrintElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
    sizeof(refPrint.refSOPClassUID), &refPrint.refSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
    sizeof(refPrint.refSOPInstanceUID), &refPrint.refSOPInstanceUID[0]},
};

static SQ_REQPROCEDURECODE reqProcedureCode;
static DCM_ELEMENT reqProcedureCodeElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
    sizeof(reqProcedureCode.codeValue), &reqProcedureCode.codeValue[0]},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(reqProcedureCode.codingSchemeDesignator),
    &reqProcedureCode.codingSchemeDesignator[0]},

    {DCM_IDCODEMEANING, DCM_LO, "", 1,
    sizeof(reqProcedureCode.codeMeaning), &reqProcedureCode.codeMeaning[0]},
};

static SQ_PROCEDURECODE procedureCode;
static DCM_ELEMENT procedureCodeElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
    sizeof(procedureCode.codeValue), &procedureCode.codeValue[0]},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(procedureCode.codingSchemeDesignator),
    &procedureCode.codingSchemeDesignator[0]},

    {DCM_IDCODEMEANING, DCM_LO, "", 1,
    sizeof(procedureCode.codeMeaning), &procedureCode.codeMeaning[0]},
};

static SQ_REFSTUDYCOMPONENT refStudyComponent;
static DCM_ELEMENT refStudyComponentElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
    sizeof(refStudyComponent.refSOPClassUID), &refStudyComponent.refSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(refStudyComponent.refSOPInstanceUID),
    &refStudyComponent.refSOPInstanceUID[0]},
};

static SQ_REQINTERPRETATIONAPPROVER reqInterpretationApprover;
static DCM_ELEMENT reqInterpretationApproverElements[] = {
    {DCM_RESINTERPAPPROVALDATE, DCM_DA, "", 1,
	sizeof(reqInterpretationApprover.approvalDates),
    &reqInterpretationApprover.approvalDates[0]},

    {DCM_RESINTERPAPPROVALTIME, DCM_TM, "", 1,
	sizeof(reqInterpretationApprover.approvalTimes),
    &reqInterpretationApprover.approvalTimes[0]},

    {DCM_RESPHYSICIANAPPROVINGINTERP, DCM_PN, "", 1,
	sizeof(reqInterpretationApprover.physiciansApproving),
    &reqInterpretationApprover.physiciansApproving[0]},
};

static SQ_REFSERIESSEQUENCE refSeriesSequence;
static DCM_ELEMENT reqRefSeriesSequenceElements[] = {
    {DCM_IDSERIESDATE, DCM_DA, "", 1,
	sizeof(refSeriesSequence.seriesDate),
    &refSeriesSequence.seriesDate[0]},

    {DCM_IDSERIESTIME, DCM_TM, "", 1,
	sizeof(refSeriesSequence.seriesTime),
    &refSeriesSequence.seriesTime[0]},

    {DCM_RELSERIESINSTANCEUID, DCM_UI, "", 1,
	sizeof(refSeriesSequence.seriesInstanceUID),
    &refSeriesSequence.seriesInstanceUID[0]}
};

static DCM_FLAGGED_ELEMENT condRefSeriesSequenceElements[] = {
    {DCM_IDRETRIEVEAETITLE, DCM_AE, "", 1,
	sizeof(refSeriesSequence.retrieveAETitle),
	&refSeriesSequence.retrieveAETitle[0],
    SQ_K_REFSERIES_RETRIEVEAETITLE, &refSeriesSequence.conditionalFields},

    {DCM_MEDIASTORAGEFILESETID, DCM_SH, "", 1,
	sizeof(refSeriesSequence.storageMediaFileSetID),
	&refSeriesSequence.storageMediaFileSetID[0],
    SQ_K_REFSERIES_STORAGEMEDIAFILESETID, &refSeriesSequence.conditionalFields},

    {DCM_MEDIASTORAGEFILESETUID, DCM_UI, "", 1,
	sizeof(refSeriesSequence.storageMediaFileSetUID),
	&refSeriesSequence.storageMediaFileSetUID[0],
    SQ_K_REFSERIES_STORAGEMEDIAFILESETUID, &refSeriesSequence.conditionalFields},

};
static SQ_REFIMAGESEQUENCE refImageSequence;
static DCM_ELEMENT reqRefImageSequenceElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(refImageSequence.refSOPClassUID),
    &refImageSequence.refSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(refImageSequence.refSOPInstanceUID),
    &refImageSequence.refSOPInstanceUID[0]}
};

static DCM_FLAGGED_ELEMENT condRefImageSequenceElements[] = {
    {DCM_IDRETRIEVEAETITLE, DCM_AE, "", 1,
	sizeof(refImageSequence.retrieveAETitle),
	&refImageSequence.retrieveAETitle[0],
    SQ_K_REFIMAGE_RETRIEVEAETITLE, &refImageSequence.conditionalFields},

    {DCM_MEDIASTORAGEFILESETID, DCM_SH, "", 1,
	sizeof(refImageSequence.storageMediaFileSetID),
	&refImageSequence.storageMediaFileSetID[0],
    SQ_K_REFIMAGE_STORAGEMEDIAFILESETID, &refImageSequence.conditionalFields},

    {DCM_MEDIASTORAGEFILESETUID, DCM_UI, "", 1,
	sizeof(refImageSequence.storageMediaFileSetUID),
	&refImageSequence.storageMediaFileSetUID[0],
    SQ_K_REFIMAGE_STORAGEMEDIAFILESETUID, &refImageSequence.conditionalFields}

};

static SQ_REFOVERLAYSEQUENCE refOverlaySequence;
static DCM_ELEMENT reqRefOverlaySequenceElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(refOverlaySequence.refSOPClassUID),
    &refOverlaySequence.refSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(refOverlaySequence.refSOPInstanceUID),
    &refOverlaySequence.refSOPInstanceUID[0]},
};

static SQ_REFCURVESEQUENCE refCurveSequence;
static DCM_ELEMENT reqRefCurveSequenceElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(refCurveSequence.refSOPClassUID),
    &refCurveSequence.refSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(refCurveSequence.refSOPInstanceUID),
    &refCurveSequence.refSOPInstanceUID[0]},
};

static SQ_CURVEREFOVERLAYSEQUENCE curveRefOverlaySequence;
static DCM_ELEMENT reqCurveRefOverlaySequenceElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(curveRefOverlaySequence.refSOPClassUID),
    &curveRefOverlaySequence.refSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(curveRefOverlaySequence.refSOPInstanceUID),
    &curveRefOverlaySequence.refSOPInstanceUID[0]},

    {DCM_CURVEREFOVERLAYGROUP, DCM_US, "", 1,
	sizeof(curveRefOverlaySequence.refOverlayGroup),
    (void *) &curveRefOverlaySequence.refOverlayGroup}
};

static SQ_PATIENTORIENTATIONCODESEQUENCE patientOrientationCodeSequence;
static DCM_ELEMENT reqPatientOrientationCodeSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(patientOrientationCodeSequence.codeValue),
    patientOrientationCodeSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(patientOrientationCodeSequence.codingSchemeDesignator),
    patientOrientationCodeSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condPatientOrientationCodeSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(patientOrientationCodeSequence.codeMeaning),
	patientOrientationCodeSequence.codeMeaning,
	SQ_K_PATIENTORIENTATIONCODESEQ_CODEMEANING,
    &patientOrientationCodeSequence.conditionalFields}
};

static SQ_PATIENTORIENTATIONMODIFIERCODESEQUENCE
    patientOrientationModifierCodeSequence;
static DCM_ELEMENT reqPatientOrientationModifierCodeSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(patientOrientationModifierCodeSequence.codeValue),
    patientOrientationModifierCodeSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(patientOrientationModifierCodeSequence.codingSchemeDesignator),
    patientOrientationModifierCodeSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condPatientOrientationModifierCodeSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(patientOrientationModifierCodeSequence.codeMeaning),
	patientOrientationModifierCodeSequence.codeMeaning,
	SQ_K_PATIENTORIENTATIONMODIFIERCODESEQ_CODEMEANING,
    &patientOrientationModifierCodeSequence.conditionalFields}
};

static SQ_PATIENTGANTRYRELATIONSHIPCODESEQUENCE
    patientGantryRelationshipCodeSequence;
static DCM_ELEMENT reqPatientGantryRelationshipCodeSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(patientGantryRelationshipCodeSequence.codeValue),
    patientGantryRelationshipCodeSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(patientGantryRelationshipCodeSequence.codingSchemeDesignator),
    patientGantryRelationshipCodeSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condPatientGantryRelationshipCodeSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(patientGantryRelationshipCodeSequence.codeMeaning),
	patientGantryRelationshipCodeSequence.codeMeaning,
	SQ_K_PATIENTGANTRYRELATIONSHIPCODESEQ_CODEMEANING,
    &patientGantryRelationshipCodeSequence.conditionalFields}
};

static SQ_ANATOMICREGIONSEQUENCE anatomicRegionSequence;
static DCM_ELEMENT reqAnatomicRegionSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(anatomicRegionSequence.codeValue),
    anatomicRegionSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(anatomicRegionSequence.codingSchemeDesignator),
    anatomicRegionSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condAnatomicRegionSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(anatomicRegionSequence.codeMeaning),
	anatomicRegionSequence.codeMeaning,
	SQ_K_ANATOMICREGIONSEQ_CODEMEANING,
    &anatomicRegionSequence.conditionalFields}
};

static SQ_ANATOMICREGIONMODIFIERSEQUENCE anatomicRegionModifierSequence;
static DCM_ELEMENT reqAnatomicRegionModifierSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(anatomicRegionModifierSequence.codeValue),
    anatomicRegionModifierSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(anatomicRegionModifierSequence.codingSchemeDesignator),
    anatomicRegionModifierSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condAnatomicRegionModifierSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(anatomicRegionModifierSequence.codeMeaning),
	anatomicRegionModifierSequence.codeMeaning,
	SQ_K_ANATOMICREGIONMODIFIERSEQ_CODEMEANING,
    &anatomicRegionModifierSequence.conditionalFields}
};

static SQ_PRIMARYANATOMICSTRUCTURESEQUENCE primaryAnatomicStructureSequence;
static DCM_ELEMENT reqPrimaryAnatomicStructureSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(primaryAnatomicStructureSequence.codeValue),
    primaryAnatomicStructureSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(primaryAnatomicStructureSequence.codingSchemeDesignator),
    primaryAnatomicStructureSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condPrimaryAnatomicStructureSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(primaryAnatomicStructureSequence.codeMeaning),
	primaryAnatomicStructureSequence.codeMeaning,
	SQ_K_PRIMARYANATOMICSTRUCTURESEQ_CODEMEANING,
    &primaryAnatomicStructureSequence.conditionalFields}
};

static SQ_PRIMARYANATOMICSTRUCTUREMODIFIERSEQUENCE
    primaryAnatomicStructureModifierSequence;
static DCM_ELEMENT reqPrimaryAnatomicStructureModifierSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(primaryAnatomicStructureModifierSequence.codeValue),
    primaryAnatomicStructureModifierSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(primaryAnatomicStructureModifierSequence.codingSchemeDesignator),
    primaryAnatomicStructureModifierSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condPrimaryAnatomicStructureModifierSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(primaryAnatomicStructureModifierSequence.codeMeaning),
	primaryAnatomicStructureModifierSequence.codeMeaning,
	SQ_K_PRIMARYANATOMICSTRUCTUREMODIFIERSEQ_CODEMEANING,
    &primaryAnatomicStructureModifierSequence.conditionalFields}
};
static SQ_ENERGYWINDOWINFOSEQUENCE energyWindowInfoSequence;
static DCM_FLAGGED_ELEMENT condEnergyWindowInfoSeqElements[] = {
    {DCM_NMIENERGYWINDOWNAME, DCM_SH, "", 1,
	sizeof(energyWindowInfoSequence.energyWindowName),
	energyWindowInfoSequence.energyWindowName,
	SQ_K_ENERGYWINDOWINFOSEQ_ENERGYWINDOWNAME,
    &energyWindowInfoSequence.conditionalFields}
};

static SQ_ENERGYWINDOWRANGESEQUENCE energyWindowRangeSequence;
static DCM_FLAGGED_ELEMENT condEnergyWindowRangeSeqElements[] = {
    {DCM_NMIENERGYWINDOWLOWERLIMIT, DCM_DS, "", 1,
	sizeof(energyWindowRangeSequence.energyWindowLowerLimit),
	energyWindowRangeSequence.energyWindowLowerLimit,
	SQ_K_ENERGYWINDOWRANGESEQ_WINDOWLOWERLIMIT,
    &energyWindowRangeSequence.conditionalFields},

    {DCM_NMIENERGYWINDOWUPPERLIMIT, DCM_DS, "", 1,
	sizeof(energyWindowRangeSequence.energyWindowUpperLimit),
	energyWindowRangeSequence.energyWindowUpperLimit,
	SQ_K_ENERGYWINDOWRANGESEQ_WINDOWUPPERLIMIT,
    &energyWindowRangeSequence.conditionalFields},
};

static SQ_RADIOPHARMACEUTICALINFOSEQUENCE radioPharmaceuticalInfoSequence;
static DCM_FLAGGED_ELEMENT condRadioPharmaceuticalInfoSeqElements[] = {
    {DCM_ACQRADIOPHARMROUTE, DCM_LO, "", 1,
	sizeof(radioPharmaceuticalInfoSequence.radioPharmaceuticalRoute),
	radioPharmaceuticalInfoSequence.radioPharmaceuticalRoute,
	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMROUTE,
    &radioPharmaceuticalInfoSequence.conditionalFields},

    {DCM_ACQRADIOPHARMVOLUME, DCM_DS, "", 1,
	sizeof(radioPharmaceuticalInfoSequence.radioPharmaceuticalVolume),
	radioPharmaceuticalInfoSequence.radioPharmaceuticalVolume,
	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMVOLUME,
    &radioPharmaceuticalInfoSequence.conditionalFields},

    {DCM_ACQRADIOPHARMSTARTTIME, DCM_TM, "", 1,
	sizeof(radioPharmaceuticalInfoSequence.radioPharmaceuticalStartTime),
	radioPharmaceuticalInfoSequence.radioPharmaceuticalStartTime,
	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMSTARTTIME,
    &radioPharmaceuticalInfoSequence.conditionalFields},

    {DCM_ACQRADIOPHARMSTOPTIME, DCM_TM, "", 1,
	sizeof(radioPharmaceuticalInfoSequence.radioPharmaceuticalStopTime),
	radioPharmaceuticalInfoSequence.radioPharmaceuticalStopTime,
	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMSTOPTIME,
    &radioPharmaceuticalInfoSequence.conditionalFields},

    {DCM_ACQRADIONUCLIDETOTALDOSE, DCM_DS, "", 1,
	sizeof(radioPharmaceuticalInfoSequence.radioNuclideTotalDose),
	radioPharmaceuticalInfoSequence.radioNuclideTotalDose,
	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIONUCLIDETOTALDOSE,
    &radioPharmaceuticalInfoSequence.conditionalFields},

    {DCM_ACQRADIOPHARMACEUTICAL, DCM_LO, "", 1,
	sizeof(radioPharmaceuticalInfoSequence.radioPharmaceutical),
	radioPharmaceuticalInfoSequence.radioPharmaceutical,
	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMACEUTICAL,
    &radioPharmaceuticalInfoSequence.conditionalFields}
};

static SQ_RADIONUCLIDECODESEQUENCE radioNuclideCodeSequence;
static DCM_ELEMENT reqRadioNuclideCodeSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(radioNuclideCodeSequence.codeValue),
    radioNuclideCodeSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(radioNuclideCodeSequence.codingSchemeDesignator),
    radioNuclideCodeSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condRadioNuclideCodeSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(radioNuclideCodeSequence.codeMeaning),
	radioNuclideCodeSequence.codeMeaning,
	SQ_K_RADIONUCLIDECODESEQ_CODEMEANING,
    &radioNuclideCodeSequence.conditionalFields}
};

static SQ_RADIOPHARMACEUTICALROUTECODESEQUENCE
    radioPharmaceuticalRouteCodeSequence;
static DCM_ELEMENT reqRadioPharmaceuticalRouteCodeSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(radioPharmaceuticalRouteCodeSequence.codeValue),
    radioPharmaceuticalRouteCodeSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(radioPharmaceuticalRouteCodeSequence.codingSchemeDesignator),
    radioPharmaceuticalRouteCodeSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condRadioPharmaceuticalRouteCodeSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(radioPharmaceuticalRouteCodeSequence.codeMeaning),
	radioPharmaceuticalRouteCodeSequence.codeMeaning,
	SQ_K_RADIOPHARMACEUTICALROUTECODESEQ_CODEMEANING,
    &radioPharmaceuticalRouteCodeSequence.conditionalFields}
};

static SQ_CALIBRATIONDATASEQUENCE calibrationDataSequence;
static DCM_ELEMENT reqCalibrationDataSeqElements[] = {
    {DCM_NMIENERGYWINDOWNUMBER, DCM_US, "", 1,
	sizeof(calibrationDataSequence.energyWindowNumber),
    (void *) &calibrationDataSequence.energyWindowNumber}
};
static DCM_FLAGGED_ELEMENT condCalibrationDataSeqElements[] = {
    {DCM_ACQSYRINGECOUNTS, DCM_IS, "", 1,
	sizeof(calibrationDataSequence.syringeCounts),
	calibrationDataSequence.syringeCounts,
	SQ_K_CALIBRATIONDATASEQ_SYRINGECOUNTS,
    &calibrationDataSequence.conditionalFields},

    {DCM_NMIRESIDUALSYRINGECOUNTS, DCM_IS, "", 1,
	sizeof(calibrationDataSequence.residualSyringeCounts),
	calibrationDataSequence.residualSyringeCounts,
	SQ_K_CALIBRATIONDATASEQ_RESIDUALSYRINGECOUNTS,
    &calibrationDataSequence.conditionalFields},
};

static SQ_RADIOPHARMACEUTICALCODESEQUENCE
    radioPharmaceuticalCodeSequence;
static DCM_ELEMENT reqRadioPharmaceuticalCodeSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(radioPharmaceuticalCodeSequence.codeValue),
    radioPharmaceuticalCodeSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(radioPharmaceuticalCodeSequence.codingSchemeDesignator),
    radioPharmaceuticalCodeSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condRadioPharmaceuticalCodeSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(radioPharmaceuticalCodeSequence.codeMeaning),
	radioPharmaceuticalCodeSequence.codeMeaning,
	SQ_K_RADIOPHARMACEUTICALCODESEQ_CODEMEANING,
    &radioPharmaceuticalCodeSequence.conditionalFields}
};

static SQ_INTERVENTIONDRUGINFOSEQUENCE interventionDrugInfoSequence;
static DCM_FLAGGED_ELEMENT condInterventionDrugInfoSeqElements[] = {
    {DCM_ACQINTERVENTIONDRUGNAME, DCM_LO, "", 1,
	sizeof(interventionDrugInfoSequence.interventionDrugName),
	interventionDrugInfoSequence.interventionDrugName,
	SQ_K_INTERVENTIONDRUGINFOSEQ_DRUGNAME,
    &interventionDrugInfoSequence.conditionalFields},

    {DCM_ACQINTERVENTIONDRUGSTART, DCM_TM, "", 1,
	sizeof(interventionDrugInfoSequence.interventionDrugStartTime),
	interventionDrugInfoSequence.interventionDrugStartTime,
	SQ_K_INTERVENTIONDRUGINFOSEQ_DRUGSTARTTIME,
    &interventionDrugInfoSequence.conditionalFields},

    {DCM_ACQINTERVENTIONDRUGSTOPTIME, DCM_TM, "", 1,
	sizeof(interventionDrugInfoSequence.interventionDrugStopTime),
	interventionDrugInfoSequence.interventionDrugStopTime,
	SQ_K_INTERVENTIONDRUGINFOSEQ_DRUGSTOPTIME,
    &interventionDrugInfoSequence.conditionalFields},

    {DCM_ACQINTERVENTIONDRUGDOSE, DCM_DS, "", 1,
	sizeof(interventionDrugInfoSequence.interventionDrugDose),
	interventionDrugInfoSequence.interventionDrugDose,
	SQ_K_INTERVENTIONDRUGINFOSEQ_DRUGDOSE,
    &interventionDrugInfoSequence.conditionalFields},
};

static SQ_INTERVENTIONDRUGCODESEQUENCE interventionDrugCodeSequence;
static DCM_ELEMENT reqInterventionDrugCodeSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(interventionDrugCodeSequence.codeValue),
    interventionDrugCodeSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(interventionDrugCodeSequence.codingSchemeDesignator),
    interventionDrugCodeSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condInterventionDrugCodeSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(interventionDrugCodeSequence.codeMeaning),
	interventionDrugCodeSequence.codeMeaning,
	SQ_K_INTERVENTIONDRUGCODESEQ_CODEMEANING,
    &interventionDrugCodeSequence.conditionalFields}
};

static SQ_DETECTORINFOSEQUENCE detectorInfoSequence;
static DCM_ELEMENT reqDetectorInfoSeqElements[] = {
    {DCM_ACQCOLLIMATORTYPE, DCM_CS, "", 1,
	sizeof(detectorInfoSequence.collimatorType),
    detectorInfoSequence.collimatorType},

    {DCM_RELIMAGEORIENTATIONPATIENT, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.imageOrientationPatient),
    detectorInfoSequence.imageOrientationPatient},

    {DCM_RELIMAGEPOSITIONPATIENT, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.imagePositionPatient),
    detectorInfoSequence.imagePositionPatient},
};
static DCM_FLAGGED_ELEMENT condDetectorInfoSeqElements[] = {
    {DCM_ACQCOLLIMATORGRIDNAME, DCM_SH, "", 1,
	sizeof(detectorInfoSequence.collimatorGridName),
	detectorInfoSequence.collimatorGridName,
	SQ_K_DETECTORINFOSEQ_COLLIMATORGRIDNAME,
    &detectorInfoSequence.conditionalFields},

    {DCM_ACQFIELDOFVIEWSHAPE, DCM_CS, "", 1,
	sizeof(detectorInfoSequence.fieldOfViewShape),
	detectorInfoSequence.fieldOfViewShape,
	SQ_K_DETECTORINFOSEQ_FIELDOFVIEWSHAPE,
    &detectorInfoSequence.conditionalFields},

    {DCM_ACQFIELDOFVIEWDIMENSION, DCM_IS, "", 1,
	sizeof(detectorInfoSequence.fieldOfViewDimensions),
	detectorInfoSequence.fieldOfViewDimensions,
	SQ_K_DETECTORINFOSEQ_FIELDOFVIEWDIMENSIONS,
    &detectorInfoSequence.conditionalFields},

    {DCM_ACQFOCALDISTANCE, DCM_IS, "", 1,
	sizeof(detectorInfoSequence.focalDistance),
	detectorInfoSequence.focalDistance,
	SQ_K_DETECTORINFOSEQ_FOCALDISTANCE,
    &detectorInfoSequence.conditionalFields},

    {DCM_ACQXFOCUSCENTER, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.xFocusCenter),
	detectorInfoSequence.xFocusCenter,
	SQ_K_DETECTORINFOSEQ_XFOCUSCENTER,
    &detectorInfoSequence.conditionalFields},

    {DCM_ACQYFOCUSCENTER, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.yFocusCenter),
	detectorInfoSequence.yFocusCenter,
	SQ_K_DETECTORINFOSEQ_YFOCUSCENTER,
    &detectorInfoSequence.conditionalFields},

    {DCM_IMGZOOMCENTER, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.zoomCenter),
	detectorInfoSequence.zoomCenter,
	SQ_K_DETECTORINFOSEQ_ZOOMCENTER,
    &detectorInfoSequence.conditionalFields},

    {DCM_IMGZOOMFACTOR, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.zoomFactor),
	detectorInfoSequence.zoomFactor,
	SQ_K_DETECTORINFOSEQ_ZOOMFACTOR,
    &detectorInfoSequence.conditionalFields},

    {DCM_ACQCENTERROTATIONOFFSET, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.centerForRotationOffset),
	detectorInfoSequence.centerForRotationOffset,
	SQ_K_DETECTORINFOSEQ_CENTERFORROTATIONOFFSET,
    &detectorInfoSequence.conditionalFields},

    {DCM_ACQGANTRYTILT, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.gantryDetectorTilt),
	detectorInfoSequence.gantryDetectorTilt,
	SQ_K_DETECTORINFOSEQ_GANTRYDETECTORTILT,
    &detectorInfoSequence.conditionalFields},

    {DCM_ACQDISTANCESRCTODETECTOR, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.distanceSrcToDetector),
	detectorInfoSequence.distanceSrcToDetector,
	SQ_K_DETECTORINFOSEQ_DISTANCESRCTODETECTOR,
    &detectorInfoSequence.conditionalFields},

    {DCM_NMISTARTANGLE, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.startAngle),
	detectorInfoSequence.startAngle,
	SQ_K_DETECTORINFOSEQ_STARTANGLE,
    &detectorInfoSequence.conditionalFields},

    {DCM_ACQRADIALPOSITION, DCM_DS, "", 1,
	sizeof(detectorInfoSequence.radialPosition),
	detectorInfoSequence.radialPosition,
	SQ_K_DETECTORINFOSEQ_RADIALPOSITION,
    &detectorInfoSequence.conditionalFields},
};

static SQ_VIEWCODESEQUENCE viewCodeSequence;
static DCM_ELEMENT reqViewCodeSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(viewCodeSequence.codeValue),
    viewCodeSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(viewCodeSequence.codingSchemeDesignator),
    viewCodeSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condViewCodeSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(viewCodeSequence.codeMeaning),
	viewCodeSequence.codeMeaning,
	SQ_K_VIEWCODESEQ_CODEMEANING,
    &viewCodeSequence.conditionalFields}
};

static SQ_VIEWCODESEQUENCE viewAngulationModifierCodeSequence;
static DCM_ELEMENT reqViewAngulationModifierCodeSeqElements[] = {
    {DCM_IDCODEVALUE, DCM_SH, "", 1,
	sizeof(viewAngulationModifierCodeSequence.codeValue),
    viewAngulationModifierCodeSequence.codeValue},

    {DCM_IDCODINGSCHEMEDESIGNATOR, DCM_SH, "", 1,
	sizeof(viewAngulationModifierCodeSequence.codingSchemeDesignator),
    viewAngulationModifierCodeSequence.codingSchemeDesignator}
};
static DCM_FLAGGED_ELEMENT condViewAngulationModifierCodeSeqElements[] = {
    {DCM_IDCODEMEANING, DCM_LO, "", 1,
	sizeof(viewAngulationModifierCodeSequence.codeMeaning),
	viewAngulationModifierCodeSequence.codeMeaning,
	SQ_K_VIEWANGULATIONMODIFIERCODESEQ_CODEMEANING,
    &viewAngulationModifierCodeSequence.conditionalFields}
};

static SQ_ROTATIONINFOSEQUENCE rotationInfoSequence;
static DCM_ELEMENT reqRotationInfoSeqElements[] = {
    {DCM_NMISTARTANGLE, DCM_DS, "", 1,
	sizeof(rotationInfoSequence.startAngle),
    rotationInfoSequence.startAngle},

    {DCM_ACQANGULARSTEP, DCM_DS, "", 1,
	sizeof(rotationInfoSequence.angularStep),
    rotationInfoSequence.angularStep},

    {DCM_ACQROTATIONDIRECTION, DCM_CS, "", 1,
	sizeof(rotationInfoSequence.rotationDirection),
    rotationInfoSequence.rotationDirection},

    {DCM_ACQSCANARC, DCM_DS, "", 1,
	sizeof(rotationInfoSequence.scanArc),
    rotationInfoSequence.scanArc},

    {DCM_ACQACTUALFRAMEDURATION, DCM_IS, "", 1,
	sizeof(rotationInfoSequence.actualFrameDuration),
    rotationInfoSequence.actualFrameDuration},

    {DCM_NMINUMBEROFFRAMESINROTATION, DCM_US, "", 1,
	sizeof(rotationInfoSequence.numberOfFramesInRotation),
    (void *) &rotationInfoSequence.numberOfFramesInRotation},
};
static DCM_FLAGGED_ELEMENT condRotationInfoSeqElements[] = {
    {DCM_ACQRADIALPOSITION, DCM_DS, "", 1,
	sizeof(rotationInfoSequence.radialPosition),
	rotationInfoSequence.radialPosition,
	SQ_K_ROTATIONINFOSEQ_RADIALPOSITION,
    &rotationInfoSequence.conditionalFields},

    {DCM_ACQDISTANCESRCTODETECTOR, DCM_DS, "", 1,
	sizeof(rotationInfoSequence.distanceSrcToDetector),
	rotationInfoSequence.distanceSrcToDetector,
	SQ_K_ROTATIONINFOSEQ_DISTANCESRCTODETECTOR,
    &rotationInfoSequence.conditionalFields},

    {DCM_ACQTABLETRAVERSE, DCM_DS, "", 1,
	sizeof(rotationInfoSequence.tableTraverse),
	rotationInfoSequence.tableTraverse,
	SQ_K_ROTATIONINFOSEQ_TABLETRAVERSE,
    &rotationInfoSequence.conditionalFields},

    {DCM_ACQTABLEHEIGHT, DCM_DS, "", 1,
	sizeof(rotationInfoSequence.tableHeight),
	rotationInfoSequence.tableHeight,
	SQ_K_ROTATIONINFOSEQ_TABLEHEIGHT,
    &rotationInfoSequence.conditionalFields},
};

static SQ_GATEDINFOSEQUENCE gatedInfoSequence;
static DCM_FLAGGED_ELEMENT condGatedInfoSeqElements[] = {
    {DCM_ACQTRIGGERTIME, DCM_DS, "", 1,
	sizeof(gatedInfoSequence.triggerTime),
	gatedInfoSequence.triggerTime,
	SQ_K_GATEDINFOSEQ_TRIGGERTIME,
    &gatedInfoSequence.conditionalFields},

    {DCM_ACQFRAMINGTYPE, DCM_LO, "", 1,
	sizeof(gatedInfoSequence.framingType),
	gatedInfoSequence.framingType,
	SQ_K_GATEDINFOSEQ_FRAMINGTYPE,
    &gatedInfoSequence.conditionalFields}
};

static SQ_DATAINFOSEQUENCE dataInfoSequence;
static DCM_ELEMENT reqDataInfoSeqElements[] = {
    {DCM_ACQFRAMETIME, DCM_DS, "", 1,
	sizeof(dataInfoSequence.frameTime),
    dataInfoSequence.frameTime}
};
static DCM_FLAGGED_ELEMENT condDataInfoSeqElements[] = {
    {DCM_ACQNOMINALINTERVAL, DCM_IS, "", 1,
	sizeof(dataInfoSequence.nominalInterval),
	dataInfoSequence.nominalInterval,
	SQ_K_DATAINFOSEQ_NOMINALINTERVAL,
    &dataInfoSequence.conditionalFields},

    {DCM_ACQLOWRRVALUE, DCM_IS, "", 1,
	sizeof(dataInfoSequence.lowRRValue),
	dataInfoSequence.lowRRValue,
	SQ_K_DATAINFOSEQ_LOWRRVALUE,
    &dataInfoSequence.conditionalFields},

    {DCM_ACQHIGHRRVALUE, DCM_IS, "", 1,
	sizeof(dataInfoSequence.highRRValue),
	dataInfoSequence.highRRValue,
	SQ_K_DATAINFOSEQ_HIGHRRVALUE,
    &dataInfoSequence.conditionalFields},

    {DCM_ACQINTERVALSACQUIRED, DCM_IS, "", 1,
	sizeof(dataInfoSequence.intervalsAcquired),
	dataInfoSequence.intervalsAcquired,
	SQ_K_DATAINFOSEQ_INTERVALSACQUIRED,
    &dataInfoSequence.conditionalFields},

    {DCM_ACQINTERVALSREJECTED, DCM_IS, "", 1,
	sizeof(dataInfoSequence.intervalsRejected),
	dataInfoSequence.intervalsRejected,
	SQ_K_DATAINFOSEQ_INTERVALSREJECTED,
    &dataInfoSequence.conditionalFields},
};

static SQ_TIMESLOTINFOSEQUENCE timeSlotInfoSequence;
static DCM_FLAGGED_ELEMENT condTimeSlotInfoSeqElements[] = {
    {DCM_NMITIMESLOTTIME, DCM_DS, "", 1,
	sizeof(timeSlotInfoSequence.timeSlotTime),
	timeSlotInfoSequence.timeSlotTime,
	SQ_K_TIMESLOTINFOSEQ_TIMESLOTTIME,
    &timeSlotInfoSequence.conditionalFields}
};

static SQ_PHASEINFOSEQUENCE phaseInfoSequence;
static DCM_ELEMENT reqPhaseInfoSeqElements[] = {
    {DCM_NMIPHASEDELAY, DCM_IS, "", 1,
	sizeof(phaseInfoSequence.phaseDelay),
    phaseInfoSequence.phaseDelay},

    {DCM_ACQACTUALFRAMEDURATION, DCM_IS, "", 1,
	sizeof(phaseInfoSequence.actualFrameDuration),
    phaseInfoSequence.actualFrameDuration},

    {DCM_NMIPAUSEBETWEENFRAMES, DCM_IS, "", 1,
	sizeof(phaseInfoSequence.pauseBetweenFrames),
    phaseInfoSequence.pauseBetweenFrames},

    {DCM_NMINUMBEROFFRAMESINPHASE, DCM_US, "", 1,
	sizeof(phaseInfoSequence.numberOfFramesInPhase),
    (void *) &phaseInfoSequence.numberOfFramesInPhase},
};
static DCM_FLAGGED_ELEMENT condPhaseInfoSeqElements[] = {
    {DCM_NMITRIGGERVECTOR, DCM_IS, "", 1,
	sizeof(phaseInfoSequence.triggerVector),
	phaseInfoSequence.triggerVector,
	SQ_K_PHASEINFOSEQ_TRIGGERVECTOR,
    &phaseInfoSequence.conditionalFields},

    {DCM_NMINUMBEROFTRIGGERSINPHASE, DCM_US, "", 1,
	sizeof(phaseInfoSequence.numberOfTriggersInPhase),
	(void *) &phaseInfoSequence.numberOfTriggersInPhase,
	SQ_K_PHASEINFOSEQ_NUMBEROFTRIGGERSINPHASE,
    &phaseInfoSequence.conditionalFields}
};

static SQ_MODALITYLUTSEQUENCE modalityLUTSequence;
static DCM_ELEMENT reqModalityLUTSequenceElements[] = {
    {DCM_IMGLUTDESCRIPTOR, DCM_US, "", 3,
	sizeof(modalityLUTSequence.lutDescriptor),
    (void *) &modalityLUTSequence.lutDescriptor},

    {DCM_IMGMODALITYLUTTYPE, DCM_LO, "", 1,
	sizeof(modalityLUTSequence.modalityLUTType),
    (void *) modalityLUTSequence.modalityLUTType}
};
static DCM_FLAGGED_ELEMENT condModalityLUTSequenceElements[] = {
    {DCM_IMGLUTEXPLANATION, DCM_LO, "", 1,
	sizeof(modalityLUTSequence.lutExplanation),
	(void *) modalityLUTSequence.lutExplanation,
	SQ_K_MODALITYLUT_LUTEXPLANATION,
    &modalityLUTSequence.conditionalFields}
};

static SQ_VOILUTSequence voiLUTSequence;
static DCM_ELEMENT reqVOILUTSequenceElements[] = {
    {DCM_IMGLUTDESCRIPTOR, DCM_US, "", 3,
	sizeof(voiLUTSequence.lutDescriptor),
    (void *) &voiLUTSequence.lutDescriptor}
};
static DCM_FLAGGED_ELEMENT condVOILUTSequenceElements[] = {
    {DCM_IMGLUTEXPLANATION, DCM_LO, "", 1,
	sizeof(voiLUTSequence.lutExplanation),
	(void *) voiLUTSequence.lutExplanation,
	SQ_K_VOILUT_LUTEXPLANATION,
    &voiLUTSequence.conditionalFields}
};

static SQ_REFERENCEDSOPSEQUENCE referencedSOPSequence;
static DCM_ELEMENT reqReferencedSOPSequenceElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(referencedSOPSequence.referencedSOPClassUID),
    (void *) referencedSOPSequence.referencedSOPClassUID},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(referencedSOPSequence.referencedSOPInstanceUID),
    (void *) referencedSOPSequence.referencedSOPInstanceUID}
};
static DCM_FLAGGED_ELEMENT condReferencedSOPSequenceElements[] = {
    {DCM_MEDIASTORAGEFILESETID, DCM_SH, "", 1,
	sizeof(referencedSOPSequence.storageMediaID),
	(void *) referencedSOPSequence.storageMediaID,
	SQ_K_REFSOP_STORAGEMEDIAID,
    &referencedSOPSequence.conditionalFields},

    {DCM_MEDIASTORAGEFILESETUID, DCM_SH, "", 1,
	sizeof(referencedSOPSequence.storageMediaUID),
	(void *) referencedSOPSequence.storageMediaUID,
	SQ_K_REFSOP_STORAGEMEDIAUID,
    &referencedSOPSequence.conditionalFields},
};

static SQ_FAILEDSOPSEQUENCE failedSOPSequence;
static DCM_ELEMENT reqFailedSOPSequenceElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(failedSOPSequence.referencedSOPClassUID),
    (void *) failedSOPSequence.referencedSOPClassUID},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(failedSOPSequence.referencedSOPInstanceUID),
    (void *) failedSOPSequence.referencedSOPInstanceUID},

    {DCM_IDFAILUREREASON, DCM_US, "", 1,
	sizeof(failedSOPSequence.failureReason),
    (void *) &failedSOPSequence.failureReason},
};

static SQ_GENERICREFERENCEDITEM genericReferencedItem;
static DCM_ELEMENT reqGenericReferencedItemElements[] = {
    {DCM_IDREFERENCEDSOPCLASSUID, DCM_UI, "", 1,
	sizeof(genericReferencedItem.referencedSOPClassUID),
    &genericReferencedItem.referencedSOPClassUID[0]},

    {DCM_IDREFERENCEDSOPINSTUID, DCM_UI, "", 1,
	sizeof(genericReferencedItem.referencedSOPInstanceUID),
    &genericReferencedItem.referencedSOPInstanceUID[0]},
};

typedef struct {
    SQ_TYPE type;
    DCM_TAG tag;
    DCM_ELEMENT *required;
    int requiredCount;
    DCM_FLAGGED_ELEMENT *conditional;
    int conditionalCount;
    void *structure;
    size_t structureSize;
}   SQ_TABLE;

static SQ_TABLE sequenceTable[] = {
    {SQ_K_REFPATIENTSOPINSTANCEUID, DCM_IDREFERENCEDPATIENTSEQ,
	refPatientInstanceUIDElements, (int) DIM_OF(refPatientInstanceUIDElements),
	NULL, 0,
    &refPatientInstanceUID, sizeof(refPatientInstanceUID)},

    {SQ_K_REFVISITSOPINSTANCEUID, DCM_IDREFERENCEDVISITSEQ,
	refVisitInstanceUIDElements, (int) DIM_OF(refVisitInstanceUIDElements),
	NULL, 0,
    &refVisitInstanceUID, sizeof(refVisitInstanceUID)},

    {SQ_K_REFSTUDYSOPINSTANCEUID, DCM_IDREFERENCEDSTUDYSEQ,
	refStudyInstanceUIDElements, (int) DIM_OF(refStudyInstanceUIDElements),
	NULL, 0,
    &refStudyInstanceUID, sizeof(refStudyInstanceUID)},

    {SQ_K_REFRESULTSSOPINSTANCEUID, DCM_IDREFERENCEDRESULTSSEQ,
	refResultsInstanceUIDElements, (int) DIM_OF(refResultsInstanceUIDElements),
	NULL, 0,
    &refResultsInstanceUID, sizeof(refResultsInstanceUID)},

    {SQ_K_REFINTERPRETATIONSOPINSTANCEUID, DCM_RESREFERENCEDINTERPSEQ,
	refInterpretationInstanceUIDElements, (int) DIM_OF(refInterpretationInstanceUIDElements),
	NULL, 0,
    &refInterpretationInstanceUID, sizeof(refInterpretationInstanceUID)},

/* Added by Anders Israelsson 		*** AI ***		*/

    {SQ_K_GREYSCALEIMAGEMODULE, DCM_BIBPREFORMATGREYSCALEIMAGESEQ,
	greyScaleImageModuleR, (int) DIM_OF(greyScaleImageModuleR),
	NULL, 0,
    &greyScaleImageModule, sizeof(greyScaleImageModule)},

    {SQ_K_COLORIMAGEMODULE, DCM_BIBPREFORMATCOLORIMAGESEQ,
	ColorImageModuleR, (int) DIM_OF(ColorImageModuleR),
	NULL, 0,
    &colorImageModule, sizeof(colorImageModule)},

/* End Added by Anders Israelsson 		*** AI ***		*/

    {SQ_K_REFFILMSESSION, DCM_BFBREFBASICFILMSESSIONSEQ,
	refPrintElements, (int) DIM_OF(refPrintElements),
	NULL, 0,
    &refPrint, sizeof(refPrint)},

    {SQ_K_REFBASICIMAGEBOX, DCM_BFBREFBASICIMAGEBOXSEQ,
	refPrintElements, (int) DIM_OF(refPrintElements),
	NULL, 0,
    &refPrint, sizeof(refPrint)},

    {SQ_K_REFPRINTJOB, DCM_PJREFPRINTJOBSEQ,
	refPrintElements, (int) DIM_OF(refPrintElements),
	NULL, 0,
    &refPrint, sizeof(refPrint)},

    {SQ_K_ADMITDIAGNOSISCODE, DCM_IDADMITDIAGCODESEQUENCE,
	reqProcedureCodeElements, (int) DIM_OF(reqProcedureCodeElements),
	NULL, 0,
    &reqProcedureCode, sizeof(reqProcedureCode)},

    {SQ_K_DISCHARGEDIAGNOSISCODE, DCM_VISDISCHARGEDIAGNOSISCODESEQ,
	reqProcedureCodeElements, (int) DIM_OF(reqProcedureCodeElements),
	NULL, 0,
    &reqProcedureCode, sizeof(reqProcedureCode)},

    {SQ_K_PROCEDURECODE, DCM_IDPROCEDURECODESEQUENCE,
	procedureCodeElements, (int) DIM_OF(procedureCodeElements),
	NULL, 0,
    &procedureCode, sizeof(procedureCode)},

    {SQ_K_REQPROCEDURECODE, DCM_SDYREQUESTEDPROCODESEQ,
	reqProcedureCodeElements, (int) DIM_OF(reqProcedureCodeElements),
	NULL, 0,
    &reqProcedureCode, sizeof(reqProcedureCode)},

    {SQ_K_REFSTUDYCOMPONENT, DCM_IDREFERENCEDSTUDYCOMPONENTSEQ,
	refStudyComponentElements, (int) DIM_OF(refStudyComponentElements),
	NULL, 0,
    &refStudyComponent, sizeof(refStudyComponent)},

    {SQ_K_REQDIAGNOSISCODE, DCM_RESDIAGNOSISCODESEQ,
	reqProcedureCodeElements, (int) DIM_OF(reqProcedureCodeElements),
	NULL, 0,
    &reqProcedureCode, sizeof(reqProcedureCode)},

    {SQ_K_REQINTERPRETATIONAPPROVER, DCM_RESINTERPAPPROVERSEQUENCE,
	reqInterpretationApproverElements,
	(int) DIM_OF(reqInterpretationApproverElements), NULL, 0,
    &reqInterpretationApprover, sizeof(reqInterpretationApprover)},

    {SQ_K_REFSERIESSEQUENCE, DCM_IDREFERENCEDSERIESSEQ,
	reqRefSeriesSequenceElements,
	(int) DIM_OF(reqRefSeriesSequenceElements),
	condRefSeriesSequenceElements,
	(int) DIM_OF(condRefSeriesSequenceElements),
    &refSeriesSequence, sizeof(refSeriesSequence)},

    {SQ_K_REFIMAGESEQUENCE, DCM_IDREFERENCEDIMAGESEQ,
	reqRefImageSequenceElements,
	(int) DIM_OF(reqRefImageSequenceElements),
	condRefImageSequenceElements,
	(int) DIM_OF(condRefImageSequenceElements),
    &refImageSequence, sizeof(refImageSequence)},

    {SQ_K_REFOVERLAYSEQUENCE, DCM_IDREFERENCEDOVERLAYSEQ,
	reqRefOverlaySequenceElements,
	(int) DIM_OF(reqRefOverlaySequenceElements),
    NULL, 0, &refOverlaySequence, sizeof(refOverlaySequence)},

    {SQ_K_REFCURVESEQUENCE, DCM_IDREFERENCEDCURVESEQ,
	reqRefCurveSequenceElements,
	(int) DIM_OF(reqRefCurveSequenceElements),
    NULL, 0, &refCurveSequence, sizeof(refCurveSequence)},

    {SQ_K_CURVEREFOVERLAYSEQUENCE, DCM_CURVEREFOVERLAYSEQUENCE,
	reqCurveRefOverlaySequenceElements,
	(int) DIM_OF(reqCurveRefOverlaySequenceElements), NULL, 0,
    &curveRefOverlaySequence, sizeof(curveRefOverlaySequence)},

    {SQ_K_PATIENTORIENTATIONCODESEQUENCE, DCM_NMIPATIENTORIENTATIONCODESEQ,
	reqPatientOrientationCodeSeqElements,
	(int) DIM_OF(reqPatientOrientationCodeSeqElements),
	condPatientOrientationCodeSeqElements,
	(int) DIM_OF(condPatientOrientationCodeSeqElements),
    &patientOrientationCodeSequence, sizeof(patientOrientationCodeSequence)},

    {SQ_K_PATIENTORIENTATIONMODIFIERCODESEQUENCE,
	DCM_NMIPATIENTORIENTATIONMODIFIERCODESEQ,
	reqPatientOrientationModifierCodeSeqElements,
	(int) DIM_OF(reqPatientOrientationModifierCodeSeqElements),
	condPatientOrientationModifierCodeSeqElements,
	(int) DIM_OF(condPatientOrientationModifierCodeSeqElements),
	&patientOrientationModifierCodeSequence,
    sizeof(patientOrientationModifierCodeSequence)},

    {SQ_K_PATIENTGANTRYRELATIONSHIPCODESEQUENCE,
	DCM_NMIPATIENTGANTRYRELATIONSHIPCODESEQ,
	reqPatientGantryRelationshipCodeSeqElements,
	(int) DIM_OF(reqPatientGantryRelationshipCodeSeqElements),
	condPatientGantryRelationshipCodeSeqElements,
	(int) DIM_OF(condPatientGantryRelationshipCodeSeqElements),
	&patientGantryRelationshipCodeSequence,
    sizeof(patientGantryRelationshipCodeSequence)},

    {SQ_K_ANATOMICREGIONSEQUENCE,
	DCM_IDANATOMICREGIONSEQUENCE,
	reqAnatomicRegionSeqElements,
	(int) DIM_OF(reqAnatomicRegionSeqElements),
	condAnatomicRegionSeqElements,
	(int) DIM_OF(condAnatomicRegionSeqElements),
	&anatomicRegionSequence,
    sizeof(anatomicRegionSequence)},

    {SQ_K_ANATOMICREGIONMODIFIERSEQUENCE,
	DCM_IDANATOMICREGIONMODIFIERSEQ,
	reqAnatomicRegionModifierSeqElements,
	(int) DIM_OF(reqAnatomicRegionModifierSeqElements),
	condAnatomicRegionModifierSeqElements,
	(int) DIM_OF(condAnatomicRegionModifierSeqElements),
	&anatomicRegionModifierSequence,
    sizeof(anatomicRegionModifierSequence)},

    {SQ_K_PRIMARYANATOMICSTRUCTURESEQUENCE,
	DCM_IDPRIMARYANATOMICSTRUCTURESEQ,
	reqPrimaryAnatomicStructureSeqElements,
	(int) DIM_OF(reqPrimaryAnatomicStructureSeqElements),
	condPrimaryAnatomicStructureSeqElements,
	(int) DIM_OF(condPrimaryAnatomicStructureSeqElements),
	&primaryAnatomicStructureSequence,
    sizeof(primaryAnatomicStructureSequence)},

    {SQ_K_PRIMARYANATOMICSTRUCTUREMODIFIERSEQUENCE,
	DCM_IDPRIMARYANATOMICSTRUCTUREMODIFIERSEQ,
	reqPrimaryAnatomicStructureModifierSeqElements,
	(int) DIM_OF(reqPrimaryAnatomicStructureModifierSeqElements),
	condPrimaryAnatomicStructureModifierSeqElements,
	(int) DIM_OF(condPrimaryAnatomicStructureModifierSeqElements),
	&primaryAnatomicStructureModifierSequence,
    sizeof(primaryAnatomicStructureModifierSequence)},

    {SQ_K_ENERGYWINDOWINFOSEQUENCE,
	DCM_NMIENERGYWINDOWINFOSEQ, NULL, 0,
	condEnergyWindowInfoSeqElements,
	(int) DIM_OF(condEnergyWindowInfoSeqElements),
    &energyWindowInfoSequence, sizeof(energyWindowInfoSequence)},

    {SQ_K_ENERGYWINDOWRANGESEQUENCE,
	DCM_NMIENERGYWINDOWRANGESEQ, NULL, 0,
	condEnergyWindowRangeSeqElements,
	(int) DIM_OF(condEnergyWindowRangeSeqElements),
    &energyWindowRangeSequence, sizeof(energyWindowRangeSequence)},

    {SQ_K_RADIOPHARMACEUTICALINFOSEQUENCE,
	DCM_NMIRADIOPHARMINFOSEQ, NULL, 0,
	condRadioPharmaceuticalInfoSeqElements,
	(int) DIM_OF(condRadioPharmaceuticalInfoSeqElements),
	&radioPharmaceuticalInfoSequence,
    sizeof(radioPharmaceuticalInfoSequence)},

    {SQ_K_RADIONUCLIDECODESEQUENCE, DCM_NMIRADIONUCLIDECODESEQUENCE,
	reqRadioNuclideCodeSeqElements,
	(int) DIM_OF(reqRadioNuclideCodeSeqElements),
	condRadioNuclideCodeSeqElements,
	(int) DIM_OF(condRadioNuclideCodeSeqElements),
    &radioNuclideCodeSequence, sizeof(radioNuclideCodeSequence)},

    {SQ_K_RADIOPHARMACEUTICALROUTECODESEQUENCE,
	DCM_NMIRADIOPHARMROUTECODESEQUENCE,
	reqRadioPharmaceuticalRouteCodeSeqElements,
	(int) DIM_OF(reqRadioPharmaceuticalRouteCodeSeqElements),
	condRadioPharmaceuticalRouteCodeSeqElements,
	(int) DIM_OF(condRadioPharmaceuticalRouteCodeSeqElements),
	&radioPharmaceuticalRouteCodeSequence,
    sizeof(radioPharmaceuticalRouteCodeSequence)},

    {SQ_K_CALIBRATIONDATASEQUENCE,
	DCM_NMICALIBRATIONDATASEQUENCE,
	reqCalibrationDataSeqElements,
	(int) DIM_OF(reqCalibrationDataSeqElements),
	condCalibrationDataSeqElements,
	(int) DIM_OF(condCalibrationDataSeqElements),
	&calibrationDataSequence,
    sizeof(calibrationDataSequence)},

    {SQ_K_RADIOPHARMACEUTICALCODESEQUENCE,
	DCM_NMIRADIOPHARMCODESEQUENCE,
	reqRadioPharmaceuticalCodeSeqElements,
	(int) DIM_OF(reqRadioPharmaceuticalCodeSeqElements),
	condRadioPharmaceuticalCodeSeqElements,
	(int) DIM_OF(condRadioPharmaceuticalCodeSeqElements),
	&radioPharmaceuticalCodeSequence,
    sizeof(radioPharmaceuticalCodeSequence)},

    {SQ_K_INTERVENTIONDRUGINFOSEQUENCE, DCM_ACQINTERVENTIONDRUGINFOSEQ,
	NULL, 0,
	condInterventionDrugInfoSeqElements,
	(int) DIM_OF(condInterventionDrugInfoSeqElements),
    &interventionDrugInfoSequence, sizeof(interventionDrugInfoSequence)},

    {SQ_K_INTERVENTIONDRUGCODESEQUENCE,
	DCM_ACQINTERVENTIONDRUGCODESEQ,
	reqInterventionDrugCodeSeqElements,
	(int) DIM_OF(reqInterventionDrugCodeSeqElements),
	condInterventionDrugCodeSeqElements,
	(int) DIM_OF(condInterventionDrugCodeSeqElements),
	&interventionDrugCodeSequence,
    sizeof(interventionDrugCodeSequence)},

    {SQ_K_DETECTORINFOSEQUENCE, DCM_NMIDETECTORINFOSEQUENCE,
	reqDetectorInfoSeqElements,
	(int) DIM_OF(reqDetectorInfoSeqElements),
	condDetectorInfoSeqElements,
	(int) DIM_OF(condDetectorInfoSeqElements),
    &detectorInfoSequence, sizeof(detectorInfoSequence)},

    {SQ_K_VIEWCODESEQUENCE,
	DCM_NMIVIEWCODESEQUENCE,
	reqViewCodeSeqElements,
	(int) DIM_OF(reqViewCodeSeqElements),
	condViewCodeSeqElements,
	(int) DIM_OF(condViewCodeSeqElements),
	&viewCodeSequence,
    sizeof(viewCodeSequence)},

    {SQ_K_VIEWANGULATIONMODIFIERCODESEQUENCE,
	DCM_NMIVIEWANGULATIONMODIFIERCODESEQ,
	reqViewAngulationModifierCodeSeqElements,
	(int) DIM_OF(reqViewAngulationModifierCodeSeqElements),
	condViewAngulationModifierCodeSeqElements,
	(int) DIM_OF(condViewAngulationModifierCodeSeqElements),
	&viewAngulationModifierCodeSequence,
    sizeof(viewAngulationModifierCodeSequence)},

    {SQ_K_ROTATIONINFOSEQUENCE, DCM_NMIROTATIONINFOSEQUENCE,
	reqRotationInfoSeqElements,
	(int) DIM_OF(reqRotationInfoSeqElements),
	condRotationInfoSeqElements,
	(int) DIM_OF(condRotationInfoSeqElements),
    &rotationInfoSequence, sizeof(rotationInfoSequence)},

    {SQ_K_GATEDINFOSEQUENCE, DCM_NMIGATEDINFOSEQUENCE,
	NULL, 0,
	condGatedInfoSeqElements, (int) DIM_OF(condGatedInfoSeqElements),
    &gatedInfoSequence, sizeof(gatedInfoSequence)},

    {SQ_K_DATAINFOSEQUENCE, DCM_NMIDATAINFORMATIONSEQUENCE,
	reqDataInfoSeqElements, (int) DIM_OF(reqDataInfoSeqElements),
	condDataInfoSeqElements, (int) DIM_OF(condDataInfoSeqElements),
    &dataInfoSequence, sizeof(dataInfoSequence)},

    {SQ_K_TIMESLOTINFOSEQUENCE, DCM_NMITIMESLOTINFOSEQUENCE,
	NULL, 0,
	condTimeSlotInfoSeqElements, (int) DIM_OF(condTimeSlotInfoSeqElements),
    &timeSlotInfoSequence, sizeof(timeSlotInfoSequence)},

    {SQ_K_PHASEINFOSEQUENCE, DCM_NMIPHASEINFOSEQUENCE,
	reqPhaseInfoSeqElements,
	(int) DIM_OF(reqPhaseInfoSeqElements),
	condPhaseInfoSeqElements,
	(int) DIM_OF(condPhaseInfoSeqElements),
    &phaseInfoSequence, sizeof(phaseInfoSequence)},

    {SQ_K_MODALITYLUTSEQUENCE, DCM_IMGMODALITYLUTSEQUENCE,
	reqModalityLUTSequenceElements,
	(int) DIM_OF(reqModalityLUTSequenceElements),
	condModalityLUTSequenceElements,
	(int) DIM_OF(condModalityLUTSequenceElements),
    &modalityLUTSequence, sizeof(modalityLUTSequence)},

    {SQ_K_VOILUTSEQUENCE, DCM_IMGVOILUTSEQUENCE,
	reqVOILUTSequenceElements,
	(int) DIM_OF(reqVOILUTSequenceElements),
	condVOILUTSequenceElements,
	(int) DIM_OF(condVOILUTSequenceElements),
    &voiLUTSequence, sizeof(voiLUTSequence)},

    {SQ_K_REFERENCEDSOPSEQUENCE, DCM_IDREFERENCEDSOPSEQUENCE,
	reqReferencedSOPSequenceElements,
	(int) DIM_OF(reqReferencedSOPSequenceElements),
	condReferencedSOPSequenceElements,
	(int) DIM_OF(condReferencedSOPSequenceElements),
    &referencedSOPSequence, sizeof(referencedSOPSequence)},

    {SQ_K_FAILEDSOPSEQUENCE, DCM_IDFAILEDSOPSEQUENCE,
	reqFailedSOPSequenceElements,
	(int) DIM_OF(reqFailedSOPSequenceElements),
	NULL, 0,
    &failedSOPSequence, sizeof(failedSOPSequence)},

    {SQ_K_REFPREVIOUSWAVEFORM, DCM_IDREFERENCEDPREVIOUSWAVEFORM,
	reqGenericReferencedItemElements,
	(int) DIM_OF(reqGenericReferencedItemElements),
	NULL, 0,
    &genericReferencedItem, sizeof(genericReferencedItem)},

    {SQ_K_REFSIMULTANEOUSWAVEFORMS, DCM_IDREFERENCEDSIMULTANEOUSWAVEFORMS,
	reqGenericReferencedItemElements,
	(int) DIM_OF(reqGenericReferencedItemElements),
	NULL, 0,
    &genericReferencedItem, sizeof(genericReferencedItem)},

    {SQ_K_REFSUBSEQUENTWAVEFORM, DCM_IDREFERENCEDSUBSEQUENTWAVEFORM,
	reqGenericReferencedItemElements,
	(int) DIM_OF(reqGenericReferencedItemElements),
	NULL, 0,
    &genericReferencedItem, sizeof(genericReferencedItem)},
};

static CONDITION
    buildSequence(LST_HEAD ** list, SQ_TABLE * sqTable, DCM_ELEMENT ** element);
static CONDITION
    parseSequence(DCM_ELEMENT * element, SQ_TABLE * sqTable, LST_HEAD ** list);

/* various static routines to build those sequences which are part of
   other sequences
*/
static CONDITION addRefImage(SQ_TABLE * sqTable, DCM_OBJECT ** obj);
static CONDITION
addPatientOrientationModifierCodeSequence(SQ_TABLE * sqTable,
					  DCM_OBJECT ** obj);
static CONDITION
    addAnatomicRegionModifierSequence(SQ_TABLE * sqTable, DCM_OBJECT ** obj);
static CONDITION
addPrimaryAnatomicStructureModifierSequence(SQ_TABLE * sqTable,
					    DCM_OBJECT ** obj);
static CONDITION
addRadioNuclideCodeSequence(SQ_TABLE * sqTable,
			    DCM_OBJECT ** obj);
static CONDITION
addRadioPharmRouteCodeSequence(SQ_TABLE * sqTable,
			       DCM_OBJECT ** obj);
static CONDITION
addCalibrationDataSequence(SQ_TABLE * sqTable,
			   DCM_OBJECT ** obj);
static CONDITION
addRadioPharmCodeSequence(SQ_TABLE * sqTable,
			  DCM_OBJECT ** obj);
static CONDITION
addEnergyWindowRangeSequence(SQ_TABLE * sqTable,
			     DCM_OBJECT ** obj);
static CONDITION
addInterventionDrugCodeSequence(SQ_TABLE * sqTable,
				DCM_OBJECT ** obj);
static CONDITION
addViewCodeSequence(SQ_TABLE * sqTable,
		    DCM_OBJECT ** obj);
static CONDITION
addViewAngulationModifierCodeSequence(SQ_TABLE * sqTable,
				      DCM_OBJECT ** obj);
static CONDITION
addDataInfoSequence(SQ_TABLE * sqTable,
		    DCM_OBJECT ** obj);
static CONDITION
addTimeSlotInfoSequence(SQ_TABLE * sqTable,
			DCM_OBJECT ** obj);
static CONDITION
addLUTData(SQ_TABLE * sqTable,
	   DCM_OBJECT ** obj);
static CONDITION
addVOILUTData(SQ_TABLE * sqTable,
	      DCM_OBJECT ** obj);

/* declare all the static routines to parse those sequences which are
   part of other sequences
*/
static CONDITION
parseRefImage(DCM_OBJECT ** obj,
	      SQ_REFSERIESSEQUENCE * series);
static CONDITION
parsePatientOrientationModifierCodeSequence(DCM_OBJECT ** obj,
				   SQ_PATIENTORIENTATIONCODESEQUENCE * seq);
static CONDITION
parseAnatomicRegionModifierSequence(DCM_OBJECT ** obj,
				    SQ_ANATOMICREGIONSEQUENCE * seq);
static CONDITION
parsePrimaryAnatomicStructureModifierSequence(DCM_OBJECT ** obj,
				 SQ_PRIMARYANATOMICSTRUCTURESEQUENCE * seq);
static CONDITION
parseEnergyWindowRangeSequence(DCM_OBJECT ** obj,
			       SQ_ENERGYWINDOWINFOSEQUENCE * seq);
static CONDITION
parseRadioNuclideCodeSequence(DCM_OBJECT ** obj,
			      SQ_RADIOPHARMACEUTICALINFOSEQUENCE * seq);
static CONDITION
parseRadioPharmRouteCodeSequence(DCM_OBJECT ** obj,
				 SQ_RADIOPHARMACEUTICALINFOSEQUENCE * seq);
static CONDITION
parseCalibrationDataSequence(DCM_OBJECT ** obj,
			     SQ_RADIOPHARMACEUTICALINFOSEQUENCE * seq);
static CONDITION
parseRadioPharmCodeSequence(DCM_OBJECT ** obj,
			    SQ_RADIOPHARMACEUTICALINFOSEQUENCE * seq);
static CONDITION
parseInterventionDrugCodeSequence(DCM_OBJECT ** obj,
				  SQ_INTERVENTIONDRUGINFOSEQUENCE * seq);
static CONDITION
parseViewCodeSequence(DCM_OBJECT ** obj,
		      SQ_DETECTORINFOSEQUENCE * seq);
static CONDITION
parseViewAngulationModifierCodeSequence(DCM_OBJECT ** obj,
					SQ_VIEWCODESEQUENCE * seq);
static CONDITION
parseDataInfoSequence(DCM_OBJECT ** obj,
		      SQ_GATEDINFOSEQUENCE * seq);
static CONDITION
parseTimeSlotInfoSequence(DCM_OBJECT ** obj,
			  SQ_DATAINFOSEQUENCE * seq);
static CONDITION
parseLUTData(DCM_OBJECT ** obj,
	     SQ_MODALITYLUTSEQUENCE * seq);
static CONDITION
parseVOILUTData(DCM_OBJECT ** obj,
		SQ_VOILUTSequence * seq);


/* SQ_BuildSequence
**
** Purpose:
**	Build a sequence type DICOM element from the specified list and the
**	given sequence type.
**
** Parameter Dictionary:
**	list		Handle to the list from which the sequence type element
**			is to be built
**	type		The "type" of sequence element to be built
**	element		The sequence type element to be built and returned to
**			the caller
**
** Return Values:
**
**	SQ_LISTCREATEFAILURE
**	SQ_LISTFAILURE
**	SQ_MALLOCFAILURE
**	SQ_MODIFICATIONFAILURE
**	SQ_NORMAL
**	SQ_OBJECTCREATEFAILED
**
** Algorithm:
**	Description of the algorithm (optional) and any other notes.
*/

CONDITION
SQ_BuildSequence(LST_HEAD ** list, SQ_TYPE type, DCM_ELEMENT ** element)
{
    int		index;
    
    for (index = 0; index < (int) DIM_OF(sequenceTable); index++) {
        
	if (type == sequenceTable[index].type)
	    return buildSequence(list, &sequenceTable[index], element);
    }
    return SQ_NORMAL;
}


/* buildSequence
**
** Purpose:
**	Build the sequence type element using the list and the information
**	about the structure of the sequence item from the entry in the
**	sequence table.
**
** Parameter Dictionary:
**	list		Handle to a list from which the sequence element is
**			to be built.
**	sqTable		Pointer to sequence table entry that supplies
**			information about the structure, the required and
**			optional attributes of the sequence item.
**	element		The sequence element to be built and returned to the
**			caller.
**
** Return Values:
**
**	SQ_LISTCREATEFAILURE
**	SQ_LISTFAILURE
**	SQ_MALLOCFAILURE
**	SQ_MODIFICATIONFAILURE
**	SQ_NORMAL
**	SQ_OBJECTCREATEFAILED
**	SQ_OBJECTACCESSFAILED
**	SQ_INTERNALSEQBUILDFAILED
**
** Notes:
**
** Algorithm:
**	Description of the algorithm (optional) and any other notes.
*/
static CONDITION
buildSequence(LST_HEAD ** list, SQ_TABLE * sqTable, DCM_ELEMENT ** element)
{
    CONDITION
    cond;
    DCM_SEQUENCE_ITEM
	* item;
    void
       *node;
    unsigned long
        length;

    *element = malloc(sizeof(**element));
    if (*element == NULL)
	return COND_PushCondition(SQ_MALLOCFAILURE,
			     SQ_Message(SQ_MALLOCFAILURE), sizeof(*element),
				  "buildSequence");

    (*element)->tag = sqTable->tag;
    (*element)->representation = DCM_SQ;
    (*element)->multiplicity = 1;
    (*element)->length = 0;
    (*element)->d.sq = LST_Create();
    if ((*element)->d.sq == NULL) {
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
			 SQ_Message(SQ_LISTCREATEFAILURE), "buildSequence");
    }
    if (*list == NULL) return SQ_NORMAL;

    node = LST_Head(list);
    if (node != NULL)
	(void) LST_Position(list, node);
printf("Position 1 \n");
    while (node != NULL) {

	item = malloc(sizeof(*item));
	if (item == NULL)
	    return COND_PushCondition(SQ_MALLOCFAILURE,
				SQ_Message(SQ_MALLOCFAILURE), sizeof(*item),
				      "buildSequence");
printf("Position 2 \n");
	cond = DCM_CreateObject(&item->object, 0);
	if (cond != DCM_NORMAL)
	    return COND_PushCondition(SQ_OBJECTCREATEFAILED,
			SQ_Message(SQ_OBJECTCREATEFAILED), "buildSequence");
printf("Position 3 \n");
	
	(void) memcpy(sqTable->structure, node, sqTable->structureSize);
printf("Position 4 \n");	

	cond = DCM_ModifyElements(&item->object, sqTable->required,
				  sqTable->requiredCount,
				  sqTable->conditional,
				  sqTable->conditionalCount,
				  NULL);
	if (cond != DCM_NORMAL)
	    return COND_PushCondition(SQ_MODIFICATIONFAILURE,
		       SQ_Message(SQ_MODIFICATIONFAILURE), "buildSequence");
printf("Position 5 \n");
	switch (sqTable->type) {
	case SQ_K_REFSERIESSEQUENCE:
	    cond = addRefImage(sqTable, &item->object);
	    break;

	case SQ_K_PATIENTORIENTATIONCODESEQUENCE:
	    cond = addPatientOrientationModifierCodeSequence(sqTable,
							     &item->object);
	    break;

	case SQ_K_ANATOMICREGIONSEQUENCE:
	    cond = addAnatomicRegionModifierSequence(sqTable, &item->object);
	    break;

	case SQ_K_PRIMARYANATOMICSTRUCTURESEQUENCE:
	    cond = addPrimaryAnatomicStructureModifierSequence(sqTable,
							     &item->object);
	    break;

	case SQ_K_ENERGYWINDOWINFOSEQUENCE:
	    cond = addEnergyWindowRangeSequence(sqTable, &item->object);
	    break;

	case SQ_K_RADIOPHARMACEUTICALINFOSEQUENCE:
printf("Position 6 \n");
	    cond = addRadioNuclideCodeSequence(sqTable, &item->object);
	    if (cond != SQ_NORMAL)
		break;
	    cond = addRadioPharmRouteCodeSequence(sqTable, &item->object);
	    if (cond != SQ_NORMAL)
		break;
	    cond = addCalibrationDataSequence(sqTable, &item->object);
	    if (cond != SQ_NORMAL)
		break;
	    cond = addRadioPharmCodeSequence(sqTable, &item->object);
	    break;

	case SQ_K_INTERVENTIONDRUGINFOSEQUENCE:
	    cond = addInterventionDrugCodeSequence(sqTable, &item->object);
	    break;

	case SQ_K_DETECTORINFOSEQUENCE:
	    cond = addViewCodeSequence(sqTable, &item->object);
	    break;

	case SQ_K_VIEWCODESEQUENCE:
	    cond = addViewAngulationModifierCodeSequence(sqTable,
							 &item->object);
	    break;
	case SQ_K_GATEDINFOSEQUENCE:
	    cond = addDataInfoSequence(sqTable, &item->object);
	    break;

	case SQ_K_DATAINFOSEQUENCE:
	    cond = addTimeSlotInfoSequence(sqTable, &item->object);
	    break;

	case SQ_K_MODALITYLUTSEQUENCE:
	    cond = addLUTData(sqTable, &item->object);
	    break;

	case SQ_K_VOILUTSEQUENCE:
	    cond = addVOILUTData(sqTable, &item->object);
	    break;

	default:
	    cond = SQ_NORMAL;
	    break;

#ifdef SMM
	    return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
		    SQ_Message(SQ_BADSEQUENCETYPEELEMENT), "buildSequence");
#endif
	}
	if (cond != SQ_NORMAL) {
	    return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
		       SQ_Message(SQ_MODIFICATIONFAILURE), "buildSequence");
	}
	cond = LST_Enqueue(&(*element)->d.sq, item);
	if (cond != LST_NORMAL)
	    return COND_PushCondition(SQ_LISTFAILURE,
			       SQ_Message(SQ_LISTFAILURE), "buildSequence");

	cond = DCM_GetObjectSize(&item->object, &length);
	if (cond != DCM_NORMAL)
	    return COND_PushCondition(SQ_OBJECTACCESSFAILED,
			SQ_Message(SQ_OBJECTACCESSFAILED), "buildSequence");
	(*element)->length += 8 + length;
	node = LST_Next(list);
    }
    return SQ_NORMAL;
}

/* SQ_ParseSequence
**
** Purpose:
**	Parse the sequence type element and decipher its type and build a
**	list of the sequence items.
**
** Parameter Dictionary:
**	element		The sequence type DICOM element to be parsed
**	type		The "type" of sequence element that is deciphered and
**			returned to the caller
**	list		Handle to the list that is created and returned to
**			the caller
**
** Return Values:
**
** Notes:
**
** Algorithm:
**	Description of the algorithm (optional) and any other notes.
*/

CONDITION
SQ_ParseSequence(DCM_ELEMENT * element, SQ_TYPE * type, LST_HEAD ** list)
{
    int
        index;


    for (index = 0; index < (int) DIM_OF(sequenceTable); index++) {
	if (element->tag == sequenceTable[index].tag) {
	    *type = sequenceTable[index].type;
	    return parseSequence(element, &sequenceTable[index], list);
	}
    }
    return SQ_NORMAL;

}

/* parseSequence
**
** Purpose:
**	Parse the sequence type DICOM element using information from the
**	sequence table entry and create a list.
**
** Parameter Dictionary:
**	element		DICOM sequence type element to be parsed
**	sqTable		Sequence table entry
**	list		Handle to the list that is built and returned to the
**			caller
**
** Return Values:
**
**	SQ_LISTCREATEFAILURE
**	SQ_LISTFAILURE
**	SQ_MALLOCFAILURE
**	SQ_NORMAL
**
** Notes:
**
** Algorithm:
**	Description of the algorithm (optional) and any other notes.
*/
static CONDITION
parseSequence(DCM_ELEMENT * element, SQ_TABLE * sqTable, LST_HEAD ** list)
{
    CONDITION
    cond;
    DCM_SEQUENCE_ITEM
	* item;
    void
       *structure;

    *list = LST_Create();
    if (*list == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
			 SQ_Message(SQ_LISTCREATEFAILURE), "parseSequence");

    item = LST_Head(&element->d.sq);
    if (item != NULL)
	(void) LST_Position(&element->d.sq, item);

    while (item != NULL) {
	memset(sqTable->structure, 0, sqTable->structureSize);
	cond = DCM_ParseObject(&item->object, sqTable->required,
			       sqTable->requiredCount,
			       sqTable->conditional,
			       sqTable->conditionalCount,
			       NULL);
	if (cond != DCM_NORMAL)
	    return cond;

	structure = malloc(sqTable->structureSize);
	if (structure == NULL)
	    return COND_PushCondition(SQ_MALLOCFAILURE,
				      SQ_Message(SQ_MALLOCFAILURE),
				   sqTable->structureSize, "parseSequence");


	(void) memcpy(structure, sqTable->structure, sqTable->structureSize);

	switch (sqTable->type) {
	case SQ_K_REFSERIESSEQUENCE:
	    cond = parseRefImage(&item->object,
				 (SQ_REFSERIESSEQUENCE *) structure);
	    break;

	case SQ_K_PATIENTORIENTATIONCODESEQUENCE:
	    cond = parsePatientOrientationModifierCodeSequence(&item->object,
			   (SQ_PATIENTORIENTATIONCODESEQUENCE *) structure);
	    break;

	case SQ_K_ANATOMICREGIONSEQUENCE:
	    cond = parseAnatomicRegionModifierSequence(&item->object,
				   (SQ_ANATOMICREGIONSEQUENCE *) structure);
	    break;

	case SQ_K_PRIMARYANATOMICSTRUCTURESEQUENCE:
	    cond = parsePrimaryAnatomicStructureModifierSequence(&item->object,
			 (SQ_PRIMARYANATOMICSTRUCTURESEQUENCE *) structure);
	    break;

	case SQ_K_ENERGYWINDOWINFOSEQUENCE:
	    cond = parseEnergyWindowRangeSequence(&item->object,
				 (SQ_ENERGYWINDOWINFOSEQUENCE *) structure);
	    break;

	case SQ_K_RADIOPHARMACEUTICALINFOSEQUENCE:
	    cond = parseRadioNuclideCodeSequence(&item->object,
			  (SQ_RADIOPHARMACEUTICALINFOSEQUENCE *) structure);
	    if (cond != SQ_NORMAL)
		break;
	    cond = parseRadioPharmRouteCodeSequence(&item->object,
			  (SQ_RADIOPHARMACEUTICALINFOSEQUENCE *) structure);
	    if ( (cond != SQ_NORMAL) && (cond != DCM_ELEMENTNOTFOUND) )
		break;
	    cond = parseCalibrationDataSequence(&item->object,
			  (SQ_RADIOPHARMACEUTICALINFOSEQUENCE *) structure);
	    if ( (cond != SQ_NORMAL) && (cond != DCM_ELEMENTNOTFOUND) )
		break;
	    cond = parseRadioPharmCodeSequence(&item->object,
			  (SQ_RADIOPHARMACEUTICALINFOSEQUENCE *) structure);
	    if ( (cond != SQ_NORMAL) && (cond != DCM_ELEMENTNOTFOUND) )
		cond = cond;
	    else cond = SQ_NORMAL;
	    break;

	case SQ_K_INTERVENTIONDRUGINFOSEQUENCE:
	    cond = parseInterventionDrugCodeSequence(&item->object,
			     (SQ_INTERVENTIONDRUGINFOSEQUENCE *) structure);
	    break;

	case SQ_K_DETECTORINFOSEQUENCE:
	    cond = parseViewCodeSequence(&item->object,
				     (SQ_DETECTORINFOSEQUENCE *) structure);
	    break;

	case SQ_K_VIEWCODESEQUENCE:
	    cond = parseViewAngulationModifierCodeSequence(&item->object,
					 (SQ_VIEWCODESEQUENCE *) structure);
	    break;

	case SQ_K_GATEDINFOSEQUENCE:
	    cond = parseDataInfoSequence(&item->object,
					 (SQ_GATEDINFOSEQUENCE *) structure);
	    break;

	case SQ_K_DATAINFOSEQUENCE:
	    cond = parseTimeSlotInfoSequence(&item->object,
					 (SQ_DATAINFOSEQUENCE *) structure);
	    break;

	case SQ_K_MODALITYLUTSEQUENCE:
	    cond = parseLUTData(&item->object,
				(SQ_MODALITYLUTSEQUENCE *) structure);
	    break;

	case SQ_K_VOILUTSEQUENCE:
	    cond = parseVOILUTData(&item->object,
				   (SQ_VOILUTSequence *) structure);
	    break;

	default:
	    cond = SQ_NORMAL;
/*
	    return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
		    SQ_Message(SQ_BADSEQUENCETYPEELEMENT), "parseSequence");
*/
	}
	if (cond != SQ_NORMAL) {
	    return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
		    SQ_Message(SQ_INTERNALSEQPARSEFAILED), "parseSequence");
	}
	cond = LST_Enqueue(list, structure);
	if (cond != LST_NORMAL)
	    return COND_PushCondition(SQ_LISTFAILURE,
			       SQ_Message(SQ_LISTFAILURE), "parseSequence");

	item = LST_Next(&element->d.sq);
    }
    return SQ_NORMAL;
}

/* SQ_ObjectToSequence
**
** Purpose:
**	Convert a DICOM object to a sequence type element
**
** Parameter Dictionary:
**	tag		Tag number of the DICOM element to be built.
**	object		The DICOM object to be converted
**	element		Handle to DICOM sequence element to be built and
**			returned to caller
**
** Return Values:
**
**	SQ_LISTCREATEFAILURE
**	SQ_LISTFAILURE
**	SQ_MALLOCFAILURE
**	SQ_OBJECTACCESSFAILED
**	SQ_NORMAL
**
** Notes:
**
** Algorithm:
**	Description of the algorithm (optional) and any other notes.
*/
CONDITION
SQ_ObjectToSequence(DCM_TAG tag, DCM_OBJECT ** object, DCM_ELEMENT ** element)
{
    CONDITION
	cond;
    DCM_SEQUENCE_ITEM
	* item;
    unsigned long
        length;

    *element = malloc(sizeof(**element));
    if (*element == NULL)
	return COND_PushCondition(SQ_MALLOCFAILURE,
			     SQ_Message(SQ_MALLOCFAILURE), sizeof(*element),
				  "SQ_ObjectToSequence");

    (*element)->tag = tag;
    (*element)->representation = DCM_SQ;
    (*element)->multiplicity = 1;
    (*element)->length = 0;
    (*element)->d.sq = LST_Create();
    if ((*element)->d.sq == NULL) {
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
		   SQ_Message(SQ_LISTCREATEFAILURE), "SQ_ObjectToSequence");
    } {
	item = malloc(sizeof(*item));
	if (item == NULL)
	    return COND_PushCondition(SQ_MALLOCFAILURE,
				SQ_Message(SQ_MALLOCFAILURE), sizeof(*item),
				      "buildSequence");

	item->object = *object;

	cond = LST_Enqueue(&(*element)->d.sq, item);
	if (cond != LST_NORMAL)
	    return COND_PushCondition(SQ_LISTFAILURE,
			       SQ_Message(SQ_LISTFAILURE), "buildSequence");

	cond = DCM_GetObjectSize(&item->object, &length);
	if (cond != DCM_NORMAL)
	    return COND_PushCondition(SQ_OBJECTACCESSFAILED,
		  SQ_Message(SQ_OBJECTACCESSFAILED), "SQ_ObjectToSequence");
	(*element)->length += 8 + length;
	*object = NULL;
    }
    return SQ_NORMAL;
}

/* SQ_SequenceToObject
**
** Purpose:
**	Convert a sequence into a DICOM object.
**
** Parameter Dictionary:
**	list		Handle to list of sequence items
**	object		DICOM object to be built and returned to caller
**
** Return Values:
**	SQ_NORMAL
**	SQ_NULLLIST
**	SQ_EMPTYLIST
**	SQ_LISTFALURE
**
** Notes:
**
** Algorithm:
**	Description of the algorithm (optional) and any other notes.
*/

CONDITION
SQ_SequenceToObject(LST_HEAD ** list, DCM_OBJECT ** object)
{
    DCM_SEQUENCE_ITEM
	* item;

    if (list == NULL)
	return COND_PushCondition(SQ_NULLLIST, SQ_Message(SQ_NULLLIST),
				  "SQ_SequenceToObject");
    if (*list == NULL)
	return COND_PushCondition(SQ_NULLLIST, SQ_Message(SQ_NULLLIST),
				  "SQ_SequenceToObject");
    if (LST_Count(list) != 1)
	return COND_PushCondition(SQ_EMPTYLIST, SQ_Message(SQ_EMPTYLIST),
				  "SQ_SequenceToObject");

    item = LST_Head(list);
    if (item == NULL)
	return COND_PushCondition(SQ_LISTFAILURE, SQ_Message(SQ_LISTFAILURE),
				  "SQ_SequenceToObject");

    *object = item->object;
    return SQ_NORMAL;
}

/* All build routines for internal sequences */
static CONDITION
addRefImage(SQ_TABLE * sqTable, DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_REFSERIESSEQUENCE *refSeries;

    refSeries = (SQ_REFSERIESSEQUENCE *) sqTable->structure;
    if (refSeries->type != SQ_K_REFSERIESSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addRefImage");

    cond = SQ_BuildSequence(&refSeries->imageList, SQ_K_REFIMAGESEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addRefImage");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}

static CONDITION
addPatientOrientationModifierCodeSequence(SQ_TABLE * sqTable,
					  DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_PATIENTORIENTATIONCODESEQUENCE *patOrientCodeSeq;

    patOrientCodeSeq = (SQ_PATIENTORIENTATIONCODESEQUENCE *) sqTable->structure;
    if (patOrientCodeSeq->type != SQ_K_PATIENTORIENTATIONCODESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
			       "addPatientOrientationModifierCodeSequence");

    cond = SQ_BuildSequence(
		  &patOrientCodeSeq->patientOrientationModifierCodeSequence,
			    SQ_K_PATIENTORIENTATIONMODIFIERCODESEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
			       "addPatientOrientationModifierCodeSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}

static CONDITION
addAnatomicRegionModifierSequence(SQ_TABLE * sqTable, DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_ANATOMICREGIONSEQUENCE *anatRegSeq;

    anatRegSeq = (SQ_ANATOMICREGIONSEQUENCE *) sqTable->structure;
    if (anatRegSeq->type != SQ_K_ANATOMICREGIONSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addAnatomicRegionModifierSequence");

    cond = SQ_BuildSequence(
			    &anatRegSeq->anatomicRegionModifierSequence,
			    SQ_K_ANATOMICREGIONMODIFIERSEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addAnatomicRegionModifierSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}

static CONDITION
addPrimaryAnatomicStructureModifierSequence(SQ_TABLE * sqTable,
					    DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_PRIMARYANATOMICSTRUCTURESEQUENCE *primAnatStructSeq;

    primAnatStructSeq =
	(SQ_PRIMARYANATOMICSTRUCTURESEQUENCE *) sqTable->structure;
    if (primAnatStructSeq->type != SQ_K_PRIMARYANATOMICSTRUCTURESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
			     "addPrimaryAnatomicStructureModifierSequence");

    cond = SQ_BuildSequence(
	       &primAnatStructSeq->primaryAnatomicStructureModifierSequence,
			 SQ_K_PRIMARYANATOMICSTRUCTUREMODIFIERSEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
			     "addPrimaryAnatomicStructureModifierSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}

static CONDITION
addEnergyWindowRangeSequence(SQ_TABLE * sqTable,
			     DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_ENERGYWINDOWINFOSEQUENCE *energyWinInfoSeq;

    energyWinInfoSeq = (SQ_ENERGYWINDOWINFOSEQUENCE *) sqTable->structure;
    if (energyWinInfoSeq->type != SQ_K_ENERGYWINDOWINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addEnergyWindowRangeSequence");

    cond = SQ_BuildSequence(
			    &energyWinInfoSeq->energyWindowRangeSequence,
			    SQ_K_ENERGYWINDOWRANGESEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addEnergyWindowRangeSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}

static CONDITION
addRadioNuclideCodeSequence(SQ_TABLE * sqTable,
			    DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_RADIOPHARMACEUTICALINFOSEQUENCE *radPharmInfoSeq;

    radPharmInfoSeq = (SQ_RADIOPHARMACEUTICALINFOSEQUENCE *) sqTable->structure;
    if (radPharmInfoSeq->type != SQ_K_RADIOPHARMACEUTICALINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addRadioNuclideCodeSequence");

	if (radPharmInfoSeq->radioNuclideCodeSequence == NULL)
	  return SQ_NORMAL;


    cond = SQ_BuildSequence(
			    &radPharmInfoSeq->radioNuclideCodeSequence,
			    SQ_K_RADIONUCLIDECODESEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addRadioNuclideCodeSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}
static CONDITION
addRadioPharmRouteCodeSequence(SQ_TABLE * sqTable,
			       DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_RADIOPHARMACEUTICALINFOSEQUENCE *radPharmInfoSeq;

    radPharmInfoSeq = (SQ_RADIOPHARMACEUTICALINFOSEQUENCE *) sqTable->structure;
    if (radPharmInfoSeq->type != SQ_K_RADIOPHARMACEUTICALINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addRadioPharmRouteCodeSequence");


	if (radPharmInfoSeq->radioPharmaceuticalRouteCodeSequence == NULL)
	  return SQ_NORMAL;


    cond = SQ_BuildSequence(
		     &radPharmInfoSeq->radioPharmaceuticalRouteCodeSequence,
			    SQ_K_RADIOPHARMACEUTICALROUTECODESEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addRadioPharmRouteCodeSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}
static CONDITION
addCalibrationDataSequence(SQ_TABLE * sqTable,
			   DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_RADIOPHARMACEUTICALINFOSEQUENCE *radPharmInfoSeq;

    radPharmInfoSeq = (SQ_RADIOPHARMACEUTICALINFOSEQUENCE *) sqTable->structure;
    if (radPharmInfoSeq->type != SQ_K_RADIOPHARMACEUTICALINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addCalibrationDataSequence");

	if (radPharmInfoSeq->calibrationDataSequence == NULL)
	  return SQ_NORMAL;


    cond = SQ_BuildSequence(
			    &radPharmInfoSeq->calibrationDataSequence,
			    SQ_K_CALIBRATIONDATASEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addCalibrationDataSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}
static CONDITION
addRadioPharmCodeSequence(SQ_TABLE * sqTable,
			  DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_RADIOPHARMACEUTICALINFOSEQUENCE *radPharmInfoSeq;

    radPharmInfoSeq = (SQ_RADIOPHARMACEUTICALINFOSEQUENCE *) sqTable->structure;
    if (radPharmInfoSeq->type != SQ_K_RADIOPHARMACEUTICALINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addRadioPharmCodeSequence");

	if (radPharmInfoSeq->radioPharmaceuticalCodeSequence == NULL)
	  return SQ_NORMAL;


    cond = SQ_BuildSequence(
			  &radPharmInfoSeq->radioPharmaceuticalCodeSequence,
			    SQ_K_RADIOPHARMACEUTICALCODESEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addRadioPharmCodeSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}
static CONDITION
addInterventionDrugCodeSequence(SQ_TABLE * sqTable,
				DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_INTERVENTIONDRUGINFOSEQUENCE *intrDrugInfoSeq;

    intrDrugInfoSeq = (SQ_INTERVENTIONDRUGINFOSEQUENCE *) sqTable->structure;
    if (intrDrugInfoSeq->type != SQ_K_INTERVENTIONDRUGINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addInterventionDrugCodeSequence");

    cond = SQ_BuildSequence(
			    &intrDrugInfoSeq->interventionDrugCodeSequence,
			    SQ_K_INTERVENTIONDRUGCODESEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addInterventionDrugCodeSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}
static CONDITION
addViewCodeSequence(SQ_TABLE * sqTable,
		    DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_DETECTORINFOSEQUENCE *detInfoSeq;

    detInfoSeq = (SQ_DETECTORINFOSEQUENCE *) sqTable->structure;
    if (detInfoSeq->type != SQ_K_DETECTORINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addViewCodeSequence");

    cond = SQ_BuildSequence(
			    &detInfoSeq->viewCodeSequence,
			    SQ_K_VIEWCODESEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addViewCodeSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}
static CONDITION
addViewAngulationModifierCodeSequence(SQ_TABLE * sqTable,
				      DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_VIEWCODESEQUENCE *viewCodeSeq;

    viewCodeSeq = (SQ_VIEWCODESEQUENCE *) sqTable->structure;
    if (viewCodeSeq->type != SQ_K_VIEWCODESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addViewAngulationModifierCodeSequence");

    cond = SQ_BuildSequence(
			    &viewCodeSeq->viewAngulationModifierCodeSequence,
			    SQ_K_VIEWANGULATIONMODIFIERCODESEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addViewAngulationModifierCodeSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}
static CONDITION
addDataInfoSequence(SQ_TABLE * sqTable,
		    DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_GATEDINFOSEQUENCE *gatedInfoSeq;

    gatedInfoSeq = (SQ_GATEDINFOSEQUENCE *) sqTable->structure;
    if (gatedInfoSeq->type != SQ_K_GATEDINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addDataInfoSequence");

    cond = SQ_BuildSequence(
			    &gatedInfoSeq->dataInfoSequence,
			    SQ_K_DATAINFOSEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addDataInfoSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}
static CONDITION
addTimeSlotInfoSequence(SQ_TABLE * sqTable,
			DCM_OBJECT ** obj)
{
    CONDITION cond;
    DCM_ELEMENT *e;
    SQ_DATAINFOSEQUENCE *dataInfoSeq;

    dataInfoSeq = (SQ_DATAINFOSEQUENCE *) sqTable->structure;
    if (dataInfoSeq->type != SQ_K_DATAINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addTimeSlotInfoSequence");

    cond = SQ_BuildSequence(
			    &dataInfoSeq->timeSlotInfoSequence,
			    SQ_K_TIMESLOTINFOSEQUENCE, &e);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQBUILDFAILED,
				  SQ_Message(SQ_INTERNALSEQBUILDFAILED),
				  "addTimeSlotInfoSequence");

    cond = DCM_AddElement(obj, e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}

static CONDITION
addLUTData(SQ_TABLE * sqTable,
	   DCM_OBJECT ** obj)
{
    CONDITION cond;
    SQ_MODALITYLUTSEQUENCE *modalityLUTSeq;
    DCM_ELEMENT e;

    modalityLUTSeq = (SQ_MODALITYLUTSEQUENCE *) sqTable->structure;
    if (modalityLUTSeq->type != SQ_K_MODALITYLUTSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addLUTData");

    memset(&e, 0, sizeof(e));
    e.tag = DCM_IMGLUTDATA;
    e.representation = DCM_US;
    e.multiplicity = 1;
    e.length = modalityLUTSeq->lutDescriptor.lutEntries *
	(modalityLUTSeq->lutDescriptor.lutDataBits / 8);
    e.d.us = modalityLUTSeq->lutData;

    cond = DCM_AddElement(obj, &e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}

static CONDITION
addVOILUTData(SQ_TABLE * sqTable,
	      DCM_OBJECT ** obj)
{
    CONDITION cond;
    SQ_VOILUTSequence *voiLUTSeq;
    DCM_ELEMENT e;

    voiLUTSeq = (SQ_VOILUTSequence *) sqTable->structure;
    if (voiLUTSeq->type != SQ_K_VOILUTSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "addVOILUTData");

    memset(&e, 0, sizeof(e));
    e.tag = DCM_IMGLUTDATA;
    e.representation = DCM_US;
    e.multiplicity = 1;
    e.length = voiLUTSeq->lutDescriptor.lutEntries *
	(voiLUTSeq->lutDescriptor.lutDataBits / 8);
    e.d.us = voiLUTSeq->lutData;

    cond = DCM_AddElement(obj, &e);
    if (cond != DCM_NORMAL)
	return cond;

    return SQ_NORMAL;
}

/* All parse routines for internal sequences */
static CONDITION
parseRefImage(DCM_OBJECT ** obj, SQ_REFSERIESSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->imageList = LST_Create();
    if (seq->imageList == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseRefImage");

    e.tag = DCM_IDREFERENCEDIMAGESEQ;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type, &seq->imageList);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
		    SQ_Message(SQ_INTERNALSEQPARSEFAILED), "parseRefImage");

    if (type != SQ_K_REFIMAGESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
		    SQ_Message(SQ_BADSEQUENCETYPEELEMENT), "parseRefImage");

    seq->conditionalFields |= SQ_K_REFSERIES_IMAGELIST;
    return SQ_NORMAL;
}
static CONDITION
parsePatientOrientationModifierCodeSequence(DCM_OBJECT ** obj,
				    SQ_PATIENTORIENTATIONCODESEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->patientOrientationModifierCodeSequence = LST_Create();
    if (seq->patientOrientationModifierCodeSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
			     "parsePatientOrientationModifierCodeSequence");

    e.tag = DCM_NMIPATIENTORIENTATIONMODIFIERCODESEQ;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type,
			    &seq->patientOrientationModifierCodeSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
			     "parsePatientOrientationModifierCodeSequence");

    if (type != SQ_K_PATIENTORIENTATIONMODIFIERCODESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
			     "parsePatientOrientationModifierCodeSequence");

    seq->conditionalFields |= SQ_K_PATIENTORIENTATIONCODESEQ_MODIFIERCODESEQ;
    return SQ_NORMAL;
}
static CONDITION
parseAnatomicRegionModifierSequence(DCM_OBJECT ** obj,
				    SQ_ANATOMICREGIONSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->anatomicRegionModifierSequence = LST_Create();
    if (seq->anatomicRegionModifierSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseAnatomicRegionModifierSequence");

    e.tag = DCM_IDANATOMICREGIONMODIFIERSEQ;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type, &seq->anatomicRegionModifierSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
				  "parseAnatomicRegionModifierSequence");

    if (type != SQ_K_ANATOMICREGIONMODIFIERSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "parseAnatomicRegionModifierSequence");

    seq->conditionalFields |= SQ_K_ANATOMICREGIONSEQ_MODIFIERSEQ;
    return SQ_NORMAL;
}
static CONDITION
parsePrimaryAnatomicStructureModifierSequence(DCM_OBJECT ** obj,
				  SQ_PRIMARYANATOMICSTRUCTURESEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->primaryAnatomicStructureModifierSequence = LST_Create();
    if (seq->primaryAnatomicStructureModifierSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
			   "parsePrimaryAnatomicStructureModifierSequence");

    e.tag = DCM_IDPRIMARYANATOMICSTRUCTUREMODIFIERSEQ;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type,
			    &seq->primaryAnatomicStructureModifierSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
			   "parsePrimaryAnatomicStructureModifierSequence");

    if (type != SQ_K_PRIMARYANATOMICSTRUCTUREMODIFIERSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
			   "parsePrimaryAnatomicStructureModifierSequence");

    seq->conditionalFields |= SQ_K_PRIMARYANATOMICSTRUCTURESEQ_MODIFIERSEQ;
    return SQ_NORMAL;
}
static CONDITION
parseEnergyWindowRangeSequence(DCM_OBJECT ** obj,
			       SQ_ENERGYWINDOWINFOSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->energyWindowRangeSequence = LST_Create();
    if (seq->energyWindowRangeSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseEnergyWindowRangeSequence");

    e.tag = DCM_NMIENERGYWINDOWRANGESEQ;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type, &seq->energyWindowRangeSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
				  "parseEnergyWindowRangeSequence");

    if (type != SQ_K_ENERGYWINDOWRANGESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "parseEnergyWindowRangeSequence");

    seq->conditionalFields |= SQ_K_ENERGYWINDOWINFOSEQ_RANGESEQ;
    return SQ_NORMAL;
}
static CONDITION
parseRadioNuclideCodeSequence(DCM_OBJECT ** obj,
			      SQ_RADIOPHARMACEUTICALINFOSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->radioNuclideCodeSequence = LST_Create();
    if (seq->radioNuclideCodeSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseRadioNuclideCodeSequence");

    e.tag = DCM_NMIRADIONUCLIDECODESEQUENCE;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type, &seq->radioNuclideCodeSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
				  "parseRadioNuclideCodeSequence");

    if (type != SQ_K_RADIONUCLIDECODESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "parseRadioNuclideCodeSequence");

    seq->conditionalFields |=
	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMNUCLIDECODESEQ;
    return SQ_NORMAL;
}
static CONDITION
parseRadioPharmRouteCodeSequence(DCM_OBJECT ** obj,
				 SQ_RADIOPHARMACEUTICALINFOSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->radioPharmaceuticalRouteCodeSequence = LST_Create();
    if (seq->radioPharmaceuticalRouteCodeSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseRadioPharmRouteCodeSequence");

    e.tag = DCM_NMIRADIOPHARMROUTECODESEQUENCE;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type,
			    &seq->radioPharmaceuticalRouteCodeSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
				  "parseRadioPharmRouteCodeSequence");

    if (type != SQ_K_RADIOPHARMACEUTICALROUTECODESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "parseRadioPharmRouteCodeSequence");

    seq->conditionalFields |=
	SQ_K_RADIOPHARMACEUTICALINFOSEQ_RADIOPHARMROUTECODESEQ;
    return SQ_NORMAL;
}

static CONDITION
parseCalibrationDataSequence(DCM_OBJECT ** obj,
			     SQ_RADIOPHARMACEUTICALINFOSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->calibrationDataSequence = LST_Create();
    if (seq->calibrationDataSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseCalibrationDataSequence");

    e.tag = DCM_NMICALIBRATIONDATASEQUENCE;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type, &seq->calibrationDataSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
				  "parseCalibrationDataSequence");

    if (type != SQ_K_CALIBRATIONDATASEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "parseCalibrationDataSequence");

    seq->conditionalFields |=
	SQ_K_RADIOPHARMACEUTICALINFOSEQ_CALIBRATIONDATASEQ;
    return SQ_NORMAL;
}

static CONDITION
parseRadioPharmCodeSequence(DCM_OBJECT ** obj,
			    SQ_RADIOPHARMACEUTICALINFOSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->radioPharmaceuticalCodeSequence = LST_Create();
    if (seq->radioPharmaceuticalCodeSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseRadioPharmCodeSequence");

    e.tag = DCM_NMIRADIOPHARMCODESEQUENCE;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type, &seq->radioPharmaceuticalCodeSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
				  "parseRadioPharmCodeSequence");

    if (type != SQ_K_RADIOPHARMACEUTICALCODESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "parseRadioPharmCodeSequence");

    seq->conditionalFields |= SQ_K_RADIOPHARMACEUTICALINFOSEQ_CODESEQ;
    return SQ_NORMAL;
}
static CONDITION
parseInterventionDrugCodeSequence(DCM_OBJECT ** obj,
				  SQ_INTERVENTIONDRUGINFOSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->interventionDrugCodeSequence = LST_Create();
    if (seq->interventionDrugCodeSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseInterventionDrugCodeSequence");

    e.tag = DCM_ACQINTERVENTIONDRUGCODESEQ;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type, &seq->interventionDrugCodeSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
				  "parseInterventionDrugCodeSequence");

    if (type != SQ_K_INTERVENTIONDRUGCODESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "parseInterventionDrugCodeSequence");

    seq->conditionalFields |= SQ_K_INTERVENTIONDRUGINFOSEQ_CODESEQ;
    return SQ_NORMAL;
}
static CONDITION
parseViewCodeSequence(DCM_OBJECT ** obj,
		      SQ_DETECTORINFOSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->viewCodeSequence = LST_Create();
    if (seq->viewCodeSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseViewCodeSequence");

    e.tag = DCM_NMIVIEWCODESEQUENCE;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type, &seq->viewCodeSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
	    SQ_Message(SQ_INTERNALSEQPARSEFAILED), "parseViewCodeSequence");

    if (type != SQ_K_VIEWCODESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
	    SQ_Message(SQ_BADSEQUENCETYPEELEMENT), "parseViewCodeSequence");

    seq->conditionalFields |= SQ_K_DETECTORINFOSEQ_VIEWCODESEQUENCE;
    return SQ_NORMAL;
}
static CONDITION
parseViewAngulationModifierCodeSequence(DCM_OBJECT ** obj,
					SQ_VIEWCODESEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->viewAngulationModifierCodeSequence = LST_Create();
    if (seq->viewAngulationModifierCodeSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseViewAngulationModifierCodeSequence");

    e.tag = DCM_NMIVIEWANGULATIONMODIFIERCODESEQ;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type,
			    &seq->viewAngulationModifierCodeSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
				  "parseViewAngulationModifierCodeSequence");

    if (type != SQ_K_VIEWANGULATIONMODIFIERCODESEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "parseViewAngulationModifierCodeSequence");

    seq->conditionalFields |= SQ_K_VIEWCODESEQ_ANGULATIONMODIFIERCODESEQ;
    return SQ_NORMAL;
}
static CONDITION
parseDataInfoSequence(DCM_OBJECT ** obj,
		      SQ_GATEDINFOSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->dataInfoSequence = LST_Create();
    if (seq->dataInfoSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseDataInfoSequence");

    e.tag = DCM_NMIDATAINFORMATIONSEQUENCE;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type, &seq->dataInfoSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
	    SQ_Message(SQ_INTERNALSEQPARSEFAILED), "parseDataInfoSequence");

    if (type != SQ_K_DATAINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
	    SQ_Message(SQ_BADSEQUENCETYPEELEMENT), "parseDataInfoSequence");

    seq->conditionalFields |= SQ_K_GATEDINFOSEQ_DATAINFOSEQ;
    return SQ_NORMAL;
}
static CONDITION
parseTimeSlotInfoSequence(DCM_OBJECT ** obj,
			  SQ_DATAINFOSEQUENCE * seq)
{
    CONDITION cond;
    DCM_ELEMENT e;
    SQ_TYPE type;

    seq->timeSlotInfoSequence = LST_Create();
    if (seq->timeSlotInfoSequence == NULL)
	return COND_PushCondition(SQ_LISTCREATEFAILURE,
				  SQ_Message(SQ_LISTCREATEFAILURE),
				  "parseTimeSlotInfoSequence");

    e.tag = DCM_NMITIMESLOTINFOSEQUENCE;
    cond = DCM_GetSequenceList(obj, e.tag, &e.d.sq);
    if (cond != DCM_NORMAL)
	return cond;

    cond = SQ_ParseSequence(&e, &type, &seq->timeSlotInfoSequence);
    if (cond != SQ_NORMAL)
	return COND_PushCondition(SQ_INTERNALSEQPARSEFAILED,
				  SQ_Message(SQ_INTERNALSEQPARSEFAILED),
				  "parseTimeSlotInfoSequence");

    if (type != SQ_K_TIMESLOTINFOSEQUENCE)
	return COND_PushCondition(SQ_BADSEQUENCETYPEELEMENT,
				  SQ_Message(SQ_BADSEQUENCETYPEELEMENT),
				  "parseTimeSlotInfoSequence");

    seq->conditionalFields |= SQ_K_DATAINFOSEQ_TIMESLOTINFOSEQ;
    return SQ_NORMAL;
}

static CONDITION
parseLUTData(DCM_OBJECT ** obj,
	     SQ_MODALITYLUTSEQUENCE * seq)
{
    DCM_ELEMENT e;
    CONDITION cond;
    void *ctx = NULL;

    e.tag = DCM_IMGLUTDATA;
    cond = DCM_GetElementSize(obj, e.tag, &e.length);
    if (cond != DCM_NORMAL)
	return 0;

    e.d.us = malloc(e.length);
    if (e.d.us == NULL)
	return 0;

    cond = DCM_GetElementValue(obj, &e, &e.length, &ctx);
    if (cond != DCM_NORMAL)
	return 0;

    seq->lutData = e.d.us;
    return SQ_NORMAL;
}

static CONDITION
parseVOILUTData(DCM_OBJECT ** obj,
		SQ_VOILUTSequence * seq)
{
    DCM_ELEMENT e;
    CONDITION cond;
    void *ctx = NULL;

    e.tag = DCM_IMGLUTDATA;
    cond = DCM_GetElementSize(obj, e.tag, &e.length);
    if (cond != DCM_NORMAL)
	return 0;

    e.d.us = malloc(e.length);
    if (e.d.us == NULL)
	return 0;

    cond = DCM_GetElementValue(obj, &e, &e.length, &ctx);
    if (cond != DCM_NORMAL)
	return 0;

    seq->lutData = e.d.us;
    return SQ_NORMAL;
}


/* SQ_ReleaseList
**
** Purpose:
**	Release list of structures allocated by SQ_ParseSequence
**
** Parameter Dictionary:
**	lst	The address of the list pointer holding the structures
**		that were allocated (and are now to be released).
**
** Return Values:
**	SQ_NORMAL
** Notes:
**
*/

CONDITION
SQ_ReleaseList(LST_HEAD ** lst)
{
    LST_NODE *node;

    while ((node = LST_Dequeue(lst)) != NULL) {
	(void) SQ_Release(&node);
    }
    (void) LST_Destroy(lst);
    return SQ_NORMAL;
}

/* SQ_Release
**
** Purpose:
**	Release one structure allocated by SQ_ParseSequence
**
** Parameter Dictionary:
**	node	Address of pointer to structure to be released.
**
** Return Values:
**	SQ_NORMAL
** Notes:
**
*/

CONDITION
SQ_Release(void **node)
{
    SQ_GENERICSTRUCT *s;

    s = (SQ_GENERICSTRUCT *) * node;

    switch (s->type) {
    case SQ_K_GREYSCALEIMAGEMODULE:
	free(((SQ_GREYSCALEIMAGEMODULE *) s)->pixelData);
	break;

    case SQ_K_COLORIMAGEMODULE:
	free(((SQ_COLORIMAGEMODULE *) s)->pixelData);
	break;

    case SQ_K_REFSERIESSEQUENCE:
	break;

    case SQ_K_PATIENTORIENTATIONCODESEQUENCE:
	break;

    case SQ_K_ANATOMICREGIONSEQUENCE:
	break;

    case SQ_K_PRIMARYANATOMICSTRUCTURESEQUENCE:
	break;

    case SQ_K_ENERGYWINDOWINFOSEQUENCE:
	break;

    case SQ_K_RADIOPHARMACEUTICALINFOSEQUENCE:
	break;

    case SQ_K_INTERVENTIONDRUGINFOSEQUENCE:
	break;

    case SQ_K_DETECTORINFOSEQUENCE:
	break;

    case SQ_K_VIEWCODESEQUENCE:
	break;

    case SQ_K_GATEDINFOSEQUENCE:
	break;

    case SQ_K_DATAINFOSEQUENCE:
	break;

    case SQ_K_MODALITYLUTSEQUENCE:
	free(((SQ_MODALITYLUTSEQUENCE *) s)->lutData);
	break;

    case SQ_K_VOILUTSEQUENCE:
	free(((SQ_VOILUTSequence *) s)->lutData);
	break;

    default:
	break;
    }
    free(s);
    *node = NULL;
    return SQ_NORMAL;
}
