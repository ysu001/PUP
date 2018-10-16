/*
 * $Id: 6ae5dac6436ed249f9978540fb895efaee288259 $
 * CUnit unit tests for dcm_to_4dfp.c
 * by Kevin A. Archie, karchie@wustl.edu
 */

#include <stdio.h>
#include <CUnit/Basic.h>

#include "dicom.h"
#include "condition.h"
#include "lst.h"
#include "dicom_objects.h"
#include "ifh.h"
#include "dcm_to_analyze.h"


/* The suite initialization function.
 * Returns zero on success, non-zero otherwise.
 */
int init_suite1(void)
{
  return 0;
}

/* The suite cleanup function.
 * Returns zero on success, non-zero otherwise.
 */
int clean_suite1(void)
{
  return 0;
}

#define BUF_SIZE 512

void test_setBase(void) {
  LST_HEAD *lh;
  PARAMS params;
  char buf[BUF_SIZE];

  strcpy(params.acqSequenceName, "*tse2dl_7");
  params.seriesNumber = 1;

  /* verify it works at all */
  lh = LST_Create();
  LST_Push(&lh, &params);
  setBase(&lh, buf);
  CU_ASSERT_STRING_EQUAL(buf, "tse2d_1");

  /* verify it works without leading *
   * and with nonzero chars already in the output buffer 
   */
  strcpy(params.acqSequenceName, "foobar328");
  params.seriesNumber = 3;
  setBase(&lh, buf);
  CU_ASSERT_STRING_EQUAL(buf, "fooba_3");

  /* how about short names? */
  memset(buf, 0, BUF_SIZE);
  strcpy(params.acqSequenceName, "*l");
  params.seriesNumber = 12;
  setBase(&lh, buf);
  CU_ASSERT_STRING_EQUAL(buf, "l_12");
}

void
test_testSelectionText(void) {
  DCM_OBJECT *obj = 0;
  CONDITION rval;
  DCM_ELEMENT e;
  DCM_TAG tag1 = DCM_IDREFERRINGPHYSICIAN; /* VR type Person's Name */
  DCM_TAG tag2 = DCM_IDOPERATORNAME; /* also VR type PN */
  char buf[1024];
  const char *testString = "Test Value";

  if (DCM_CreateObject(&obj, 0) != DCM_NORMAL || obj == 0)
    CU_FAIL("Unable to create DCM object");

  /* we reuse this element; okay because AddElement makes deep copy */
  memset(&e, 0, sizeof(e));
  e.tag = tag1;
  e.representation = DCM_PN;
  e.d.string = buf;
  e.length = sizeof(buf) - 1;

  if (DCM_AddElement(&obj, &e) != DCM_NORMAL)
    CU_FAIL("Unable to add element 1 to DCM object");
  rval = testSelectionText(&obj, tag1, testString);
  CU_ASSERT_EQUAL(rval, 1);

  e.tag = tag2;
  strcpy(buf, testString);
  if (DCM_AddElement(&obj, &e) != DCM_NORMAL)
    CU_FAIL("Unable to add element 2 to DCM object");
  rval = testSelectionText(&obj, tag1, testString);
  CU_ASSERT_EQUAL(rval, 1);	/* (0000,0001) still has empty text */
  rval = testSelectionText(&obj, tag2, testString);
  CU_ASSERT_EQUAL(rval, 0);	/* (0000,0002) has our text */
}



void
test_checkOrientation(void) {
  PARAMS paramsTra;
  PARAMS paramsCor;
  PARAMS paramsSag;
  PARAMS dummy1;
  PARAMS dummy2;
  float orientTra[] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
  float orientCor[] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
  float orientSag[] = { 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
  const size_t osz = sizeof(orientTra);
  LST_HEAD *lh;
  int rval;

  typedef enum {Tra = 2, Cor = 3, Sag = 4} Orientation;

  if (!(lh = LST_Create()))
    CU_FAIL("could not create list structure");

  // checkOrientation() somewhat mysteriously requires at least three
  // entries in the parameter list, and uses only the second one.
  LST_Push(&lh, &dummy2);
  
  memcpy(&paramsTra.imageOrientation, orientTra, osz);
  LST_Push(&lh, &paramsTra);
  memcpy(&paramsCor.imageOrientation, orientCor, osz);
  LST_Push(&lh, &paramsCor);
  memcpy(&paramsSag.imageOrientation, orientSag, osz);
  LST_Push(&lh, &paramsSag);

  LST_Push(&lh, &dummy1);

  CU_ASSERT_EQUAL(Sag, checkOrientation(&lh));
  LST_Dequeue(&lh);
  CU_ASSERT_EQUAL(Cor, checkOrientation(&lh));
  LST_Dequeue(&lh);
  CU_ASSERT_EQUAL(Tra, checkOrientation(&lh));
  LST_Dequeue(&lh);

  // if there aren't at least three items in the list,
  // checkOrientation() returns an error indication.
  fprintf(stderr, "Ignore this Volume Count warning: ");
  CU_ASSERT_EQUAL(5, checkOrientation(&lh));
}

  
void test_parseParams(void) {
  int i;
  DCM_OBJECT *obj;
  DCM_ELEMENT e;
  const char *thickness = "foo thick";
  PARAMS p, p2;

  // pick numbers we're unlikely to get by accident
  const U16 valRows = 125;
  const U16 valCols = 123;
  const U16 valBitsAllocated = 5;
  const U16 valBitsStored = 3;
  const U16 valHighBit = 7;
  const U16 valPixelRep = 0;
  const char *valSeriesInstanceUID = "1.2.3.4.5.6.7.8";
  const char *valInstanceNumber = "17";
  const char *valPhotometInterp = "MONOCHROME2";
  const char *valSeriesNumber = "19";
  const char *valPixelSpacing = "1.2\\1.2";
  const char *valSliceThickness = "9.8";
  const char *valImagePosition = "1.0\\-1.0\\0.0";
  const char *valImageOrientation = "1\\0\\0\\0\\0\\1";
  const U16 valPlanarConfig = 0;
  const char *valStudyDate = "20061012";
  const char *valStudyDesc = "this was a very nice study";
  const char *valProtocol = "a great protocol";
  const char *valPatientName = "Barr^Foo^X";
  const char *valPatientID = "bf32";
  const U16 valLargestPix = 252;
  const U16 valSmallestPix = 3;
  const char *valSliceSpacing = "8.0";
  const char *valSeqName = "my sequence";
  const char *valManufacturer = "Acme";
  const char *valModelName = "BirdSeed";
  const char *valModality = "MR";
  const char *valInstitution = "WUSM";
  const char *specifiedThickness = "10.2";
  
  
  /*
   * List of tags taken from dcm_to_4dfp
   *
   * Attribute name, multiplicity, and VR individually checked against
   * DICOM standard (some disagreements with the code)
   */
  DCM_ELEMENT ELEMS_TO_PARSE[] = {
    { DCM_IMGROWS, DCM_US, "Rows", 1, sizeof(U16), (void*)&valRows },
    { DCM_IMGCOLUMNS, DCM_US, "Columns", 1, sizeof(U16), (void*)&valCols },
    { DCM_IMGBITSALLOCATED, DCM_US, "Bits Allocated", 1, sizeof(U16),
      (void*)&valBitsAllocated },
    { DCM_IMGBITSSTORED, DCM_US, "Bits Stored", 1, sizeof(U16), 
      (void*)&valBitsStored },
    { DCM_IMGHIGHBIT, DCM_US, "High Bit", 1, sizeof(U16), (void*)&valHighBit },
    { DCM_IMGPIXELREPRESENTATION, DCM_US, "Pixel Representation", 1,
      sizeof(U16), (void*)&valPixelRep },
    { DCM_RELSERIESINSTANCEUID, DCM_UI, "Series Instance UID", 1,
      DICOM_UI_LENGTH, (void*)valSeriesInstanceUID },
    { DCM_RELIMAGENUMBER, DCM_IS, "Instance Number", 1, DICOM_IS_LENGTH,
      (void*)valInstanceNumber },
    { DCM_IMGPHOTOMETRICINTERP, DCM_CS, "Photometric Interpretation", 1,
      DICOM_CS_LENGTH, (void*)valPhotometInterp },
    { DCM_RELSERIESNUMBER, DCM_IS, "Series Number", 1, DICOM_IS_LENGTH,
      (void*)valSeriesNumber},
    { DCM_IMGPIXELSPACING, DCM_DS, "Pixel Spacing", 2, DICOM_DS_LENGTH*2,
      (void*)valPixelSpacing},
    { DCM_ACQSLICETHICKNESS, DCM_DS, "Slice Thickness", 1, DICOM_DS_LENGTH,
      (void*)valSliceThickness },
    { DCM_RELIMAGEPOSITIONPATIENT, DCM_DS, "Image Position (Patient)", 3,
      DICOM_DS_LENGTH*3, (void*)valImagePosition },
    { DCM_RELIMAGEORIENTATIONPATIENT, DCM_DS, "Image Orientation (Patient)",
      6, DICOM_DS_LENGTH*6, (void*)valImageOrientation },
    { DCM_IMGPLANARCONFIGURATION, DCM_US, "Planar Configuration", 1,
      sizeof(U16), (void*)&valPlanarConfig },
    { DCM_IDSTUDYDATE, DCM_DA, "Study Date", 1, DICOM_DA_LENGTH,
      (void*)valStudyDate },
    { DCM_IDSERIESDESCR, DCM_LO, "Study Description", 1,
      DICOM_LO_LENGTH, (void*)valStudyDesc },
    { DCM_ACQPROTOCOLNAME, DCM_LO, "Protocol Name", 1, DICOM_LO_LENGTH,
      (void*)valProtocol },
    { DCM_PATNAME, DCM_PN, "Patient's Name", 1, DICOM_PN_LENGTH,
      (void*)valPatientName },
    { DCM_PATID, DCM_LO, "Patient ID", 1, DICOM_LO_LENGTH,
      (void*)valPatientID },
    { DCM_IMGLARGESTIMAGEPIXELVALUE, DCM_US, "Largest Image Pixel Value",
      1, sizeof(U16), (void*)&valLargestPix },
    { DCM_IMGSMALLESTIMAGEPIXELVALUE, DCM_US, "Smallest Image Pixel Value",
      1, sizeof(U16), (void*)&valSmallestPix },
    { DCM_ACQSLICESPACING, DCM_DS, "Spacing Between Slices", 1,
      DICOM_DS_LENGTH, (void*)valSliceSpacing },
    { DCM_ACQSEQUENCENAME, DCM_SH, "Sequence Name", 1, DICOM_SH_LENGTH,
      (void*)valSeqName },
    { DCM_IDMANUFACTURER, DCM_LO, "Manufacturer", 1, DICOM_LO_LENGTH,
      (void*)valManufacturer },
    { DCM_IDMANUFACTURERMODEL, DCM_LO, "Manufacturer's Model Name", 1,
      DICOM_LO_LENGTH, (void*)valModelName },
    { DCM_IDMODALITY, DCM_CS, "Modality", 1, DICOM_CS_LENGTH,
      (void*)valModality },
    { DCM_IDINSTITUTIONNAME, DCM_LO, "Institution Name", 1,
      DICOM_LO_LENGTH, (void*)valInstitution }
  };
  int nElems = sizeof(ELEMS_TO_PARSE)/sizeof(DCM_ELEMENT);

  if ((DCM_CreateObject(&obj, 0) != DCM_NORMAL) || !obj)
    CU_FAIL("Unable to create DCM object");

  for (i = 0; i < nElems; i++)
    if (DCM_AddElement(&obj, ELEMS_TO_PARSE+i) != DCM_NORMAL)
      CU_FAIL("Unable to add element to DCM object");

  CU_ASSERT_EQUAL(0, parseParams(&obj, 0, &p));

  CU_ASSERT_EQUAL(valRows, p.rows);
  CU_ASSERT_EQUAL(valCols, p.columns);
  CU_ASSERT_EQUAL(valBitsAllocated, p.bitsAllocated);
  CU_ASSERT_EQUAL(valBitsStored, p.bitsStored);
  CU_ASSERT_EQUAL(valHighBit, p.highBit);
  CU_ASSERT_EQUAL(valPixelRep, p.pixelRepresentation);
  CU_ASSERT_STRING_EQUAL(valSeriesInstanceUID, p.seriesUID);
  CU_ASSERT_EQUAL(atoi(valInstanceNumber), p.instanceNumber);
  CU_ASSERT_STRING_EQUAL(valPhotometInterp, p.photometricInterpretation);
  CU_ASSERT_EQUAL(atoi(valSeriesNumber), p.seriesNumber);
  CU_ASSERT_STRING_EQUAL(valPixelSpacing, p.pixelSpacing);
  CU_ASSERT_STRING_EQUAL(valSliceThickness, p.sliceThickness);
  {
    double px, py, pz;
    int n = sscanf(valImagePosition, "%lf\\%lf\\%lf", &px, &py, &pz);
    CU_ASSERT_EQUAL(3, n);
    CU_ASSERT_EQUAL(px, p.positionX);
    CU_ASSERT_EQUAL(py, p.positionY);
    CU_ASSERT_EQUAL(pz, p.positionZ);
  }
  {
    int i;
    double o[6];
    int n = sscanf(valImageOrientation, "%lf\\%lf\\%lf\\%lf\\%lf\\%lf",
		   &o[0], &o[1], &o[2], &o[3], &o[4], &o[5]);
    CU_ASSERT_EQUAL(6, n);
    for (i = 0; i < 6; i++)
      CU_ASSERT_EQUAL(o[i], p.imageOrientation[i]);
  }

  CU_ASSERT_EQUAL(valPlanarConfig, p.planarConfiguration);
  CU_ASSERT_STRING_EQUAL(valStudyDate, p.studyDate);
  CU_ASSERT_STRING_EQUAL(valStudyDesc, p.seriesDescription);
  CU_ASSERT_STRING_EQUAL(valProtocol, p.protocolName);
  CU_ASSERT_STRING_EQUAL(valPatientName, p.patientName);
  CU_ASSERT_STRING_EQUAL(valPatientID, p.patientID);

  /* These two parameters are calculated rather
   * than read from the DICOM information.
   *  CU_ASSERT_EQUAL(valLargestPix, p.imgLargestPix);
   *  CU_ASSERT_EQUAL(valSmallestPix, p.imgSmallestPix);
   */

  CU_ASSERT_STRING_EQUAL(valSliceSpacing, p.sliceSpacing);
  CU_ASSERT_STRING_EQUAL(valSeqName, p.acqSequenceName);
  CU_ASSERT_STRING_EQUAL(valManufacturer, p.idManufacturer);
  CU_ASSERT_STRING_EQUAL(valModelName, p.idModel);
  CU_ASSERT_STRING_EQUAL(valModality, p.idModality);
  CU_ASSERT_STRING_EQUAL(valInstitution, p.idInstitution);

  /* do one more with a specified slice thickness */
  CU_ASSERT_EQUAL(0, parseParams(&obj, specifiedThickness, &p2));

  /* should be all the same as p except for slice thickness */
  CU_ASSERT_STRING_EQUAL(specifiedThickness, p2.sliceThickness);
  CU_ASSERT_EQUAL(valRows, p2.rows);
  CU_ASSERT_EQUAL(valCols, p2.columns);
  CU_ASSERT_EQUAL(valBitsAllocated, p2.bitsAllocated);
  CU_ASSERT_EQUAL(valBitsStored, p2.bitsStored);
  CU_ASSERT_EQUAL(valHighBit, p2.highBit);
  CU_ASSERT_EQUAL(valPixelRep, p2.pixelRepresentation);
  CU_ASSERT_STRING_EQUAL(valSeriesInstanceUID, p2.seriesUID);
  CU_ASSERT_EQUAL(atoi(valInstanceNumber), p2.instanceNumber);
  CU_ASSERT_STRING_EQUAL(valPhotometInterp, p2.photometricInterpretation);
  CU_ASSERT_EQUAL(atoi(valSeriesNumber), p2.seriesNumber);
  CU_ASSERT_STRING_EQUAL(valPixelSpacing, p2.pixelSpacing);
  {
    double px, py, pz;
    int n = sscanf(valImagePosition, "%lf\\%lf\\%lf", &px, &py, &pz);
    CU_ASSERT_EQUAL(3, n);
    CU_ASSERT_EQUAL(px, p2.positionX);
    CU_ASSERT_EQUAL(py, p2.positionY);
    CU_ASSERT_EQUAL(pz, p2.positionZ);
  }
  {
    int i;
    double o[6];
    int n = sscanf(valImageOrientation, "%lf\\%lf\\%lf\\%lf\\%lf\\%lf",
		   &o[0], &o[1], &o[2], &o[3], &o[4], &o[5]);
    CU_ASSERT_EQUAL(6, n);
    for (i = 0; i < 6; i++)
      CU_ASSERT_EQUAL(o[i], p2.imageOrientation[i]);
  }

  CU_ASSERT_EQUAL(valPlanarConfig, p2.planarConfiguration);
  CU_ASSERT_STRING_EQUAL(valStudyDate, p2.studyDate);
  CU_ASSERT_STRING_EQUAL(valStudyDesc, p2.seriesDescription);
  CU_ASSERT_STRING_EQUAL(valProtocol, p2.protocolName);
  CU_ASSERT_STRING_EQUAL(valPatientName, p2.patientName);
  CU_ASSERT_STRING_EQUAL(valPatientID, p2.patientID);

  /* These two parameters are calculated rather
   * than read from the DICOM information.
   *  CU_ASSERT_EQUAL(valLargestPix, p2.imgLargestPix);
   *  CU_ASSERT_EQUAL(valSmallestPix, p2.imgSmallestPix);
   */

  CU_ASSERT_STRING_EQUAL(valSliceSpacing, p2.sliceSpacing);
  CU_ASSERT_STRING_EQUAL(valSeqName, p2.acqSequenceName);
  CU_ASSERT_STRING_EQUAL(valManufacturer, p2.idManufacturer);
  CU_ASSERT_STRING_EQUAL(valModelName, p2.idModel);
  CU_ASSERT_STRING_EQUAL(valModality, p2.idModality);
  CU_ASSERT_STRING_EQUAL(valInstitution, p2.idInstitution);
}

void test_testRangeText() {
  DCM_OBJECT *obj = 0;
  CONDITION rval;
  DCM_ELEMENT e;
  DCM_TAG tagIS = DCM_RELACQUISITIONNUMBER; /* VR type Integer String */
  /* should also try: DS, FL, FD, SL, SS, UL, US, maybe DA, TM, etc. */

  char buf[1024];

  if (DCM_CreateObject(&obj, 0) != DCM_NORMAL || obj == 0)
    CU_FAIL("Unable to create DCM object");

  /* we reuse this element; okay because AddElement makes deep copy */
  memset(&e, 0, sizeof(e));
  e.tag = tagIS;
  e.representation = DCM_IS;
  e.d.string = buf;
  sprintf(buf, "-5");
  e.length = sizeof(buf) - 1;

  if (DCM_AddElement(&obj, &e) != DCM_NORMAL)
    CU_FAIL("Unable to add element 1 to DCM object");

  CU_ASSERT_EQUAL(0, testRangeText(&obj, tagIS, "-10", "-1"));
  CU_ASSERT_EQUAL(0, testRangeText(&obj, tagIS, "-5", "-5"));
  CU_ASSERT_EQUAL(1, testRangeText(&obj, tagIS, "-10", "-8"));
  CU_ASSERT_EQUAL(1, testRangeText(&obj, tagIS, "-4", "4"));
  CU_ASSERT_EQUAL(1, testRangeText(&obj, tagIS, "-4", "+4"));
  CU_ASSERT_EQUAL(1, testRangeText(&obj, tagIS, "2", "8"));

  DCM_CloseObject(&obj);
  if (DCM_CreateObject(&obj, 0) != DCM_NORMAL || obj == 0)
    CU_FAIL("Unable to create DCM object");

  sprintf(buf, "0");
  if (DCM_AddElement(&obj, &e) != DCM_NORMAL)
    CU_FAIL("Unable to add element 1 to DCM object");
  
  CU_ASSERT_EQUAL(1, testRangeText(&obj, tagIS, "-10", "-1"));
  CU_ASSERT_EQUAL(0, testRangeText(&obj, tagIS, "0", "0"));
  CU_ASSERT_EQUAL(1, testRangeText(&obj, tagIS, "1", "10"));
  CU_ASSERT_EQUAL(0, testRangeText(&obj, tagIS, "-4", "+4"));

  DCM_CloseObject(&obj);
  if (DCM_CreateObject(&obj, 0) != DCM_NORMAL || obj == 0)
    CU_FAIL("Unable to create DCM object");

  sprintf(buf, "42");
  if (DCM_AddElement(&obj, &e) != DCM_NORMAL)
    CU_FAIL("Unable to add element 1 to DCM object");

  CU_ASSERT_EQUAL(1, testRangeText(&obj, tagIS, "-45", "-40"));
  CU_ASSERT_EQUAL(1, testRangeText(&obj, tagIS, "-45", "40"));
  CU_ASSERT_EQUAL(0, testRangeText(&obj, tagIS, "40", "45"));
  CU_ASSERT_EQUAL(0, testRangeText(&obj, tagIS, "+40", "45"));
  CU_ASSERT_EQUAL(0, testRangeText(&obj, tagIS, "+40", "+45"));
}


/* The main() function for setting up and running the tests.
 * Returns a CUE_SUCCESS on successful running, another
 * CUnit error code on failure.
 */
int main() {
   CU_pSuite pSuite = NULL;

   /* initialize the CUnit test registry */
   if (CUE_SUCCESS != CU_initialize_registry())
      return CU_get_error();

   /* add a suite to the registry */
   pSuite = CU_add_suite("Suite_1", init_suite1, clean_suite1);
   if (NULL == pSuite) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* add the tests to the suite */
   if (!CU_add_test(pSuite, "setBase()", test_setBase) ||
       !CU_add_test(pSuite, "testSelectionText()", test_testSelectionText) ||
       !CU_add_test(pSuite, "checkOrientation()", test_checkOrientation) ||
       !CU_add_test(pSuite, "parseParams()", test_parseParams) ||
       !CU_add_test(pSuite, "testRangeText()", test_testRangeText))
     {
       CU_cleanup_registry();
       return CU_get_error();
     }

   /* Run all tests using the CUnit Basic interface */
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   CU_cleanup_registry();
   return CU_get_error();
}
