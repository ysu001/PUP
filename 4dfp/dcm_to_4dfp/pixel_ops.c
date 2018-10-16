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
** Last Update:		$Author: karchie $, $Date: 2007/09/21 00:38:48 $
** Source File:		$RCSfile: pixel_ops.c,v $
** Revision:		$Revision: 1.2 $
** Status:		$State: Exp $
*/

static char rcsid[] = "$Revision: 1.2 $ $RCSfile: pixel_ops.c,v $";

#if HAVE_MALLOC
#include <malloc.h>
#endif
#include <string.h>
#include <stdlib.h>

#include "include/dicom.h"
#include "include/ctnthread.h"
#include "include/condition.h"
#include "include/lst.h"
#include "include/utility.h"
#include "include/dicom_objects.h"

#include "include/dcm_to_analyze.h"

/******************************************************************************/
static void* extractFPixelsMONOCHROME2(DCM_OBJECT** obj, PARAMS *p)
{
  void 		*pixels = NULL;
  void 		*fpixels = NULL;
  
  int 		byteCount;
  int 		pixelCount;
  
  DCM_ELEMENT 	e;
  CONDITION 	cond;

  S16* 		p1;	/* 16 bit pointer for Pixel Rep = 1 (pet)  */
  U16*		p0;	/* 16 bit data Pixel Rep = 0 (mr)  */
  float* 	pf;	/* Pointer to the float data   */
  int 		index;	/* pixel index  */

  pixelCount = p->rows * p->columns;
  byteCount = p->rows * p->columns;
  
  if (p->bitsAllocated > 8)byteCount *= 4; /* 4 bytes are in 32 bits */

  if (pixels == NULL) {
    pixels = malloc(byteCount);  /* 16 bit pixels */
    if (pixels == NULL) {
      fprintf(stderr, "Could not allocated %d bytes for pixel data\n",
	      byteCount);
      exit(7);
    }
  }
  
  memset(&e, 0, sizeof(e));
  e.tag = DCM_PXLPIXELDATA;
  e.representation = DCM_OW;
  e.length = byteCount;
  e.d.ot = pixels;

  cond = DCM_ParseObject(obj, &e, 1, NULL, 0, NULL);
  if (cond != DCM_NORMAL  && cond != DCM_GETINCOMPLETE) {
    fprintf(stderr, "Could not extract pixels: %d %d %d %d\n",
	    p->rows, p->columns, p->bitsAllocated, byteCount);
    COND_DumpConditions();
    exit(8);
  }

  pf = malloc(pixelCount * sizeof(float));
  p1 = pixels;
  p0 = pixels;
  if (p->pixelRepresentation == 1){
      for (index = 0; index < pixelCount; index++) {
           pf[index] = p1[index];
      }
  }else{
      for (index = 0; index < pixelCount; index++) {
           pf[index] = p0[index];
      }
  }
  free(pixels);
  fpixels = pf;		/* fpixels now points to pf */ 
  
  return fpixels;	/* fpixels should be freed in the calling routine */
  
}

/******************************************************************************/

static void* extractPixelsMONOCHROME2(DCM_OBJECT** obj, PARAMS *p)
{
  void *pixels = NULL;
  int byteCount;
  int pixelCount;
  DCM_ELEMENT e;
  CONDITION cond;

  pixelCount = p->rows * p->columns;
  byteCount = p->rows * p->columns;
  if (p->bitsAllocated > 8)
    byteCount *= 2;

  if (pixels == NULL) {
    pixels = malloc(byteCount);
    if (pixels == NULL) {
      fprintf(stderr, "Could not allocated %d bytes for pixel data\n",
	      byteCount);
      exit(7);
    }
  }
  memset(&e, 0, sizeof(e));
  e.tag = DCM_PXLPIXELDATA;
  e.representation = DCM_OW;
  e.length = byteCount;
  e.d.ot = pixels;

  cond = DCM_ParseObject(obj, &e, 1, NULL, 0, NULL);
  if (cond != DCM_NORMAL  && cond != DCM_GETINCOMPLETE) {
    fprintf(stderr, "Could not extract pixels: %d %d %d %d\n",
	    p->rows, p->columns, p->bitsAllocated, byteCount);
    COND_DumpConditions();
    exit(8);
  }

  if (p->bitsAllocated == 8) {
    unsigned char* p1;
    U16* p2;
    int index;

    p2 = malloc(pixelCount * 2);
    p1 = (unsigned char*) pixels;
    for (index = 0; index < pixelCount; index++) {
      p2[index] = p1[index];
    }

    free(pixels);
    pixels = p2;
  }

  return pixels;
}

 /******************************************************************************/

void* extractPixels(DCM_OBJECT** obj, PARAMS *p)
{
  
  void *pixels = 0;
 
 
  
  if (strcmp(p->photometricInterpretation, "MONOCHROME2") == 0) {
       pixels = extractFPixelsMONOCHROME2(obj, p);
  } else {
       fprintf(stderr, "photometricInterpretation ERROR\n Not ready for Floating Point conversion of: %s\n",
       p->photometricInterpretation);
       exit(4);
  }
  return pixels;
  
}
