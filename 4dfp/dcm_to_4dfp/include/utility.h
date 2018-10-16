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
**		     Electronic Radiology Laboratory
**		   Mallinckrodt Institute of Radiology
**		Washington University School of Medicine
**
** Module Name(s):	utility.h
** Author, Date:	David E. Beecher, March 1994
** Intent:		Define typedefs and function prototypes for
**			Utility (UTL) facility (for functions which may
**			generally useful in a number of areas).
** Last Update:		$Author: mohanar $, $Date: 2006/08/25 20:55:50 $
** Source File:		$RCSfile: utility.h,v $
** Revision:		$Revision: 1.1 $
** Status:		$State: Exp $
*/


#ifndef _UTL_IS_IN
#define _UTL_IS_IN 1

#include "lst.h"

#ifdef  __cplusplus
extern "C" {
#endif

#define OFF		0
#define ON		1
#define REGEX_SIZE	128

#if !defined ( LINUX ) && !defined ( IRIX )
char *re_comp(char *);
int re_exec(char *);
#endif

typedef struct {
  void* reserved[2];
  char path[1024];
} UTL_FILEITEM;

char *UTL_ConvertRegex(char *regex);

long UTL_ConvertDatetoLong(const char *date);
double UTL_ConvertTimetoFloat(const char *time);

void UTL_ConvertLongtoDate(long ld, char *date);
void UTL_ConvertFloattoTime(double dt, char *time);
void UTL_SqueezeBlanks(char *s);

void UTL_GetDicomDate(char *date);
void UTL_GetDicomTime(char *time);

CONDITION UTL_RegexMatch(char *regex, char *stm);
CONDITION UTL_DateMatch(char *datestring, char *stm);
CONDITION UTL_TimeMatch(char *timestring, char *stm);

void* UTL_GetTimeStamp();
double UTL_DeltaTime(void* timeStamp);
void UTL_ReleaseTimeStamp(void* timeStamp);

CONDITION UTL_VerifyCreatePath(const char* path);
CTNBOOLEAN UTL_IsDirectory(const char* path);
CTNBOOLEAN UTL_IsFile(const char* path);
CONDITION UTL_DeleteFile(const char* path);
CONDITION UTL_ScanDirectory(const char* path, LST_HEAD** lst);

CONDITION UTL_SetConfigFile(const char* configFile);
CONDITION UTL_TestConfigFile(const char* configFile);
char* UTL_GetConfigParameter(const char* paramName);
char**
UTL_ExpandToPointerArray(const char* inputText,
			 const char* delimiters,
			 int* numberOfEntries);
CONDITION UTL_FileSize(const char* path, U32* size);

#define	UTL_NORMAL		FORM_COND(FAC_UTL, SEV_SUCC, 1)
#define	UTL_UNIMPLEMENTED	FORM_COND(FAC_UTL, SEV_ERROR, 2)
#define UTL_MATCH		FORM_COND(FAC_UTL, SEV_SUCC, 3)
#define UTL_NOMATCH		FORM_COND(FAC_UTL, SEV_SUCC, 4)
#define UTL_PATHNOTDIR		FORM_COND(FAC_UTL, SEV_ERROR, 5)
#define UTL_FILECREATEFAILED	FORM_COND(FAC_UTL, SEV_ERROR, 6)
#define UTL_NO_CTN_TARGET	FORM_COND(FAC_UTL, SEV_ERROR, 7)
#define UTL_DELETEFILEFAILED	FORM_COND(FAC_UTL, SEV_ERROR, 8)
#ifdef  __cplusplus
}
#endif

#endif
