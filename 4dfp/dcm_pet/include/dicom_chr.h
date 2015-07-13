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
** Author, Date:	Stephen M. Moore, 1-Jan-1999
** Intent:		This file contains definitions and function prototypes
**			for the CHR facility, used to manipulate character
**			sets.
** Last Update:		$Author: smm $, $Date: 1999/11/22 17:40:01 $
** Source File:		$RCSfile: dicom_chr.h,v $
** Revision:		$Revision: 1.3 $
** Status:		$State: Exp $
*/

#ifndef DCM_CHR_IS_IN
#define DCM_CHR_IS_IN 1

#ifdef  __cplusplus
extern "C" {
#endif

typedef enum { CHR_ISO2022JP, CHR_EUC_JP , CHR_EUC_JPROMAJI,
               CHR_DEFAULTISOIR13,
	       CHR_ISO2022KR, CHR_EUC_KR } CHR_ENCODING;

typedef enum { CHR_ASCII,
	       CHR_ISOIR13,  CHR_ISOIR14, CHR_ISOIR87, CHR_ISOIR159,
	       CHR_ISOIR149
	} CHR_CHARACTER;

typedef struct {
  const unsigned char* s;
  int sLength;
  int sOffset;
  int sMode;
  CHR_ENCODING encoding;
} CHR_ITERATOR_CONTEXT;

#define CHR_ESCAPE 0x1B
#define CHR_SS2    0x8E
#define CHR_SS3    0x8F

CONDITION CHR_Translate(const void* src, int srcLength, CHR_ENCODING srcCode,
			void* dest, int destLength, CHR_ENCODING destCode,
			int* encodedLength);
CONDITION CHR_Encode(const void* src, int srcLength, CHR_ENCODING srcCode,
		     void* dest, int destLength, CHR_ENCODING destCode,
		     int* encodedLength);
CONDITION CHR_IterateBegin(const void* src, int length,
			   CHR_ENCODING code, CHR_ITERATOR_CONTEXT* ctx);
CONDITION CHR_NextCharacter(CHR_ITERATOR_CONTEXT* ctx, CHR_CHARACTER* charSet,
			    void* s, int* length);

CONDITION CHR_ValidateEncoding(const void* src, int srcLength,
			       CHR_ENCODING srcCode, CHR_ENCODING destCode,
			       int verbose);

char* CHR_Message(CONDITION condition);


#define CHR_NORMAL		/* Normal return from CHR package */ \
	FORM_COND(FAC_CHR, SEV_SUCC, 1)
#define	CHR_UNSUPPORTEDTRANSLATION \
	FORM_COND(FAC_CHR, SEV_ERROR, 2)
#define	CHR_ILLEGALCHARACTERSET \
	FORM_COND(FAC_CHR, SEV_ERROR, 3)

#ifdef  __cplusplus
}
#endif

#endif
