/*$Header: /data/petsun4/data1/src_solaris/dcm_dump_file/RCS/dcm_dump_file.c,v 1.2 2012/12/05 01:45:24 avi Exp $*/
/*$Log: dcm_dump_file.c,v $
 * Revision 1.2  2012/12/05  01:45:24  avi
 * new-style main to suppress compile warnings
 *
 * Revision 1.1  2012/12/04  04:05:53  avi
 * Initial revision
 **/
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
** Module Name(s):	main
**			usageerror
** Author, Date:	Stephen M. Moore, 9-May-93
** Intent:		This program uses the DICOM OBJECTS package to open
**			DICOM files and dump their contents to stdout.  Each
**			argument to the program is expected to be the name
**			of a file containing a DICOM stream.
**   Usage:
**			dcm_dump_file [-b] [-g] [-v] [-z] file [file ...]
** Last Update:		$Author: avi $, $Date: 2012/12/05 01:45:24 $
** Source File:		$RCSfile: dcm_dump_file.c,v $
** Revision:		$Revision: 1.2 $
** Status:		$State: Exp $
*/

static char rcsid[] = "$Revision: 1.2 $ $RCSfile: dcm_dump_file.c,v $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef GCCSUNOS
#include <sys/types.h>
#endif
#ifdef MALLOC_DEBUG
#include <malloc.h>
#endif

#include "include/dicom.h"
#include "include/ctnthread.h"
#include "include/lst.h"
#include "include/condition.h"
#include "include/dicom_objects.h"
#include "include/utility.h"

static void usageerror();
typedef struct {
  char path[1024];
} PATH_STRUCT;

static void
fillFileList(const char* f, LST_HEAD** lst)
{
  if (UTL_IsDirectory(f)) {
    LST_HEAD* l = 0;
    UTL_FILEITEM* item;

    UTL_ScanDirectory(f, &l);
    item = (UTL_FILEITEM*)LST_Dequeue(&l);
    while(item != NULL) {
      char p2[1024];
      if ((strcmp(item->path, ".") != 0) && (strcmp(item->path, "..") != 0)) {
	strcpy(p2, f);
	strcat(p2, "/");
	strcat(p2, item->path);
	fillFileList(p2, lst);
      }
      free(item);
      item = (UTL_FILEITEM*)LST_Dequeue(&l);
    }
  } else {
    UTL_FILEITEM * p;
    p = malloc(sizeof(*p));
    strcpy(p->path, f);
    LST_Enqueue(lst, p);
  }
}

int main (int argc, char **argv)
{
    DCM_OBJECT
	* object;
    CONDITION
	cond;
    CTNBOOLEAN
	verbose = FALSE,
	exitFlag = FALSE,
	formatFlag = FALSE;
    unsigned long
        options = DCM_ORDERLITTLEENDIAN;
    long vmLimit = 0;
    LST_HEAD* fileNames = 0;
    UTL_FILEITEM* p = NULL;

    while (--argc > 0 && (*++argv)[0] == '-') {
	switch (*(argv[0] + 1)) {
	case 'b':
	    options &= ~DCM_ORDERMASK;
	    options |= DCM_ORDERBIGENDIAN;
	    break;
	case 'e':
	    exitFlag = TRUE;
	    break;
	case 'E':
	    options &= ~DCM_FILEFORMATMASK;
	    options |= DCM_PART10FILE;
	    options &= ~DCM_SPECIALFORMATMASK;
	    options |= DCM_EFILM;
	    break;
	case 'f':
	    formatFlag = TRUE;
	    break;
	case 'g':
	    options &= ~DCM_GROUPLENGTHMASK;
	    options |= DCM_NOGROUPLENGTH;
	    break;
	case 'l':
	    options &= ~DCM_LENGTHTOENDMASK;
	    options |= DCM_USELENGTHTOEND;
	    break;
	case 'L':
	    options &= ~DCM_ORDERMASK;
	    options |= DCM_EXPLICITLITTLEENDIAN;
	    break;
	case 'm':
	    argc--;
	    argv++;
	    vmLimit = atoi(*argv);
	    break;
	case 't':
	    options &= ~DCM_FILEFORMATMASK;
	    options |= DCM_PART10FILE;
	    options |= DCM_ACCEPTVRMISMATCH;
	    break;

	case 'v':
	    verbose = TRUE;
	    break;
	case 'w':
	    argc--; argv++;
	    if (argc < 1) {
		usageerror();
	    }
	    if (strcmp(*argv, "REPEAT") == 0)
		options |= DCM_ALLOWREPEATELEMENTS;
	    break;
	case 'z':
	    options |= DCM_FORMATCONVERSION;
	    break;
	default:
	    break;
	}
    }

    THR_Init();
    DCM_Debug(verbose);
    if (argc < 1)
	usageerror();

    fileNames = LST_Create();
    while (argc-- > 0) {
	fillFileList(*argv, &fileNames);
	argv++;
    }

    p = LST_Dequeue(&fileNames);
    while (p != NULL) {
	printf("DICOM File: %s\n", p->path);
	cond = DCM_OpenFile(p->path, options, &object);
	if (cond != DCM_NORMAL && ((options & DCM_PART10FILE) == 0)) {
	    COND_DumpConditions();
	    (void) DCM_CloseObject(&object);
	    (void) COND_PopCondition(TRUE);
	    fprintf(stderr, "Could not open %s as expected.  Trying Part 10 format.\n", p->path);
	    cond = DCM_OpenFile(p->path, options | DCM_PART10FILE, &object);
	}
	if (cond == DCM_NORMAL) {
	    if (formatFlag)
		cond = DCM_FormatElements(&object, vmLimit, "");
	    else
		cond = DCM_DumpElements(&object, vmLimit);
	}
	COND_DumpConditions();
	(void) DCM_CloseObject(&object);
	(void) COND_PopCondition(TRUE);

	if (cond != DCM_NORMAL && exitFlag) {
	    THR_Shutdown();
	    exit(1);
	}
	free(p);
	p = LST_Dequeue(&fileNames);
    }
#ifdef MALLOC_DEBUG
    malloc_verify(0);
    malloc_shutdown();
#endif
    THR_Shutdown();
    return 0;
}

/* usageerror
**
** Purpose:
**	Print the usage string for this application and exit.
**
** Parameter Dictionary:
**	None
**
** Return Values:
**	None
**
** Algorithm:
**	Description of the algorithm (optional) and any other notes.
*/

static void
usageerror()
{
    char msg[] = "\
Usage: dcm_dump_file [-b] [-g] [-l] [-m mult] [-t] [-v] [-w flag] [-z] file [file ...]\n\
\n\
    -b        Input files are stored in big-endian byte order\n\
    -e        Exit on file open error.  Do not process other files\n\
    -g        Remove group length elements\n\
    -l        Use (retired) length-to-end attribute for object length\n\
    -m mult   Change VM limit from 0 to mult\n\
    -t        Part 10 file\n\
    -v        Place DCM facility in verbose mode\n\
    -w        Set open options; flag can be REPEAT \n\
    -z        Perform format conversion (verification) on data in files\n\
\n\
    file      One or more information objects to be dumped\n";

    fprintf(stderr, msg);
    exit(1);
}
