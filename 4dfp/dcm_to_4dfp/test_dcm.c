/*
 * $Id: f1776966c408c02e9ce8d6e403267fff5c494eda $
 * CUnit unit tests for dcm.c
 * by Kevin A. Archie, karchie@wustl.edu
 */

#include <stdio.h>
#include <CUnit/Basic.h>

#include "dicom.h"
#include "condition.h"
#include "lst.h"
#include "dicom_objects.h"
#include "dicom_uids.h"
#include "dcmprivate.h"

extern CONDITION updateSpecialElements(PRIVATE_OBJECT **object,
				       PRV_ELEMENT_ITEM *item);
extern CONDITION newElementItem(DCM_ELEMENT *src, CTNBOOLEAN allocateData,
				PRV_ELEMENT_ITEM **dst);


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

void test_updateSpecialElements(void) {
  PRIVATE_OBJECT *obj = 0;
  DCM_ELEMENT e;
  PRV_ELEMENT_ITEM *item = 0;
  DCM_TAG tag = DCM_METATRANSFERSYNTAX;
  DCM_VALUEREPRESENTATION vr = DCM_UI;
  CONDITION cond;
  char buf[BUF_SIZE];

  if (DCM_CreateObject((void*)&obj, 0) != DCM_NORMAL || obj == 0)
    CU_FAIL("Unable to create DCM object");
  

  memset(&e, 0, sizeof(e));
  e.tag = tag;
  e.representation = vr;
  e.d.string = buf;
  e.length = sizeof(buf) - 1;

  strcpy(buf, DICOM_TRANSFERLITTLEENDIAN);
  if (newElementItem(&e, FALSE, &item) != DCM_NORMAL || item == 0)
    CU_FAIL("Unable to create DCM Element Item");
  item->element.d.string = buf;
  cond = updateSpecialElements(&obj, item);
  CU_ASSERT_EQUAL(DCM_NORMAL, cond);
  CU_ASSERT_EQUAL(DCM_ORDERLITTLEENDIAN, obj->dataOptions);
  
  strcpy(buf, DICOM_TRANSFERLITTLEENDIANEXPLICIT);
  cond = updateSpecialElements(&obj, item);
  CU_ASSERT_EQUAL(DCM_NORMAL, cond);
  CU_ASSERT_EQUAL(DCM_EXPLICITLITTLEENDIAN, obj->dataOptions);

  strcpy(buf, DICOM_TRANSFERBIGENDIANEXPLICIT);
  cond = updateSpecialElements(&obj, item);
  CU_ASSERT_EQUAL(DCM_NORMAL, cond);
  CU_ASSERT_EQUAL(DCM_EXPLICITBIGENDIAN, obj->dataOptions);

  /* this one will generate warnings that don't interest us. */
  strcpy(buf, DICOM_TRANSFERJPEGEXTENDEDPROC2AND4);
  fprintf(stderr, "\nIgnore the following warning about transfer syntax:\n");
  cond = updateSpecialElements(&obj, item);
  CU_ASSERT_EQUAL(DCM_NORMAL, cond);
  CU_ASSERT_EQUAL(DCM_EXPLICITLITTLEENDIAN, obj->dataOptions);
}

void test_checkOrientation(void) {
  CU_FAIL("test unimplemented");
}


/* The main() function for setting up and running the tests.
 * Returns a CUE_SUCCESS on successful running, another
 * CUnit error code on failure.
 */
int main()
{
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
   if (!CU_add_test(pSuite, "updateSpecialElements()", test_updateSpecialElements))
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
