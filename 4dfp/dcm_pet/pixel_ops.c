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
** Original Author, Date:	Stephen M. Moore, 14-Oct-2001
** Intent:		
** Last Update:		$Author: jon $, $Date: 2008/03/26 20:23:19 $
** Source File:		$RCSfile: pixel_ops.c,v $
** Revision:		$Revision: 1.2 $
** Status:		$State: Exp $
*/

static char rcsid[] = "$Revision: 1.2 $ $RCSfile: pixel_ops.c,v $";

#include <malloc.h>
#include <stdlib.h>
#include <string.h>

#include "include/dicom.h"
#include "include/ctnthread.h"
#include "include/condition.h"
#include "include/lst.h"
#include "include/utility.h"
#include "include/dicom_objects.h"
/*#include "include/analyze.h"*/
#include "include/dcm_pet.h"

/******************************************************************************/
static void* extractFPixelsMONOCHROME2(DCM_OBJECT** obj, PARAMS *p)
{
  void 		*pixels = NULL;
  void 		*fpixels = NULL;
  
  int 		byteCount;
  int 		pixelCount;
  
  DCM_ELEMENT 	e;
  CONDITION 	cond;

  S16* 		p1;	/* 16 bit pointer for Pixel Rep = 1  */
  U16*		p0;	/* 16 bit pointer for Pixel Rep = 0  */
  float* 	pf;	/* Pointer to the float data   */
  int 		index;	/* pixel index  */

  if(p->img_number_frames > 0){
     pixelCount = p->rows * p->columns * p->img_number_frames;
     byteCount =  p->rows * p->columns * p->img_number_frames;
     if (p->bitsAllocated > 8)byteCount *= 2; /* 2 bytes are in 32 bits */
  }else{
     pixelCount = p->rows * p->columns;
     byteCount =  p->rows * p->columns;
     if (p->bitsAllocated > 8)byteCount *= 2;
  }
  
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
  
  printf("/n In extractPixelsMONOCHROME2 /n");
  
  if(p->img_number_frames > 0){
     /*pixelCount = p->rows * p->columns * p->img_number_frames;*/
     /*byteCount =  p->rows * p->columns * p->img_number_frames;*/
     pixelCount = p->rows * p->columns;
     byteCount = p->rows * p->columns;
  }else{
     pixelCount = p->rows * p->columns;
     byteCount = p->rows * p->columns;
  }
  
  if (p->bitsAllocated > 8)
    byteCount *= 2;

  if (pixels == NULL) {
    /*pixels = malloc(byteCount * sizeof(U16));*/
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
  CTNBOOLEAN floatFlag = TRUE;
  void *pixels = 0;
  
  if(floatFlag) {
  
     if (strcmp(p->photometricInterpretation, "MONOCHROME2") == 0) {
       pixels = extractFPixelsMONOCHROME2(obj, p);
     } else {
       fprintf(stderr, "pixelOps: photometricInterpretation ERROR\n Not ready for conversion of: %s\n in File %s",
       p->photometricInterpretation,p->fileName);
       exit(4);
     }
     
  } else if (strcmp(p->photometricInterpretation, "MONOCHROME2") == 0) {
        pixels = extractPixelsMONOCHROME2(obj, p);
  } else {
        fprintf(stderr, "Not ready to handle conversion for: %s\n",
	p->photometricInterpretation);
        exit(4);
  }
  
  return pixels;
  
}
