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
** Module Name(s):	common_disp.h (main())
** Author, Date:	David E. Beecher, 30-Jun-93
** Intent:		Display server for the imaging demonstration
** Last Update:		$Author: smm $, $Date: 2002/12/13 15:17:52 $
** Source File:		$RCSfile: common_disp.h,v $
** Revision:		$Revision: 1.10 $
** Status:		$State: Exp $
*/

#ifndef _COMMMONDISP_IS_IN
#define _COMMMONDISP_IS_IN

#ifdef CTN_MULTIBYTE
#include "dicom_chr.h"
#endif

#define		DEFAULT			0
#define		SLOPE_INTERCEPT		1
#define		WINDOW_LEVEL		2
#define		S_I_W_L			3
#define         MONOCHROME1             4
#define         MONOCHROME2             5
#define         PALETTECOLOR            6
#define         CTNRGB                  7
#define         HSV                     8
#define         CTNRGBA                 9
#define         CTNCMYK                 10
#define         NOT_IMPLEMENTED         11

#define		MAX_STRING_LENGTH	1024

#define RED_COMP        0.299
#define GREEN_COMP      0.587
#define BLUE_COMP       0.114

typedef struct _ImageStruct {
#ifdef CTN_MULTIBYTE
    CHR_ENCODING encoding;      /* Encoding used for patient name, ID .. */
#endif
    int
        sampsperpixel,		/* Samples per pixel             */
        photointerp,		/* The photometric interpretation */
        nrows,			/* Number of rows in image       */
        ncols,			/* Number of columns in image    */
        bitsallocated,		/* Bits allocated per pixel      */
        bitsstored,		/* Bits stored per pixel         */
        highbit,		/* High bit for pixels           */
        pixrep,			/* Pixel representation          */
        xpos,			/* X position on display         */
        ypos,			/* Y position on display         */
        smethod;		/* The scaling method            */
    float
        m,			/* The slope                     */
        b,			/* The y intercept               */
        window,			/* The Window width              */
        level;			/* The Window Level              */
    char
        study[MAX_STRING_LENGTH],
        patname[MAX_STRING_LENGTH],
        date[MAX_STRING_LENGTH],
        modality[MAX_STRING_LENGTH],
        patstring[MAX_STRING_LENGTH];
    int
        error_condition;
    char
        error_message[MAX_STRING_LENGTH];
    unsigned char
       *data;			/* the image data		 */
}   ImageStruct;
#if 0
typedef struct _XImageStruct {
    XImage *x_image;		/* the X image struct		 */
    Window image_window;	/* the window for display	 */
}   XImageStruct;
#endif

/* This flag is used if the user does not want to overide Width and Center
** values when retrieving the image attributes.
*/

#define NO_OVERRIDE	0x80000
/*
 * Function Prototypes
 */
unsigned char
   *ScaleImageData(void *data, ImageStruct * idata);
void
   *
GetDICOMData(const char *filename, unsigned long options, ImageStruct * idata,
	     int overrideWidth, int overrideCenter, CTNBOOLEAN useImport);

#endif
