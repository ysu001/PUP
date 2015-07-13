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
** Module Name(s):	SQ_Message
** Author, Date:	Stephen M. Moore, 21-Jul-93
** Intent:		Define the ASCIZ messages that go with each SQ
**			error number and provide a function for looking up
**			the error message.
** Last Update:		$Author: jon $, $Date: 2008/04/03 13:47:31 $
** Source File:		$RCSfile: sqcond.c,v $
** Revision:		$Revision $
** Status:		$State: Exp $
*/

static char rcsid[] = "$Revision: 1.1 $ $RCSfile: sqcond.c,v $";

#include <stdio.h>
#include <sys/types.h>
#include "dicom.h"
#include "lst.h"
#include "dicom_objects.h"
#include "dicom_sq.h"

typedef struct vector {
    CONDITION cond;
    char *message;
}   VECTOR;

static VECTOR messageVector[] = {
    {SQ_NORMAL, "SQ Normal return from SQ routine"},
    {SQ_MALLOCFAILURE, "SQ Function failed to allocate %d bytes (%s)"},
    {SQ_LISTCREATEFAILURE, "SQ Function failed to create a new list (%s)"},
    {SQ_OBJECTCREATEFAILED, "SQ Function failed to create DCM object (%s)"},
    {SQ_MODIFICATIONFAILURE, "SQ Function failed to modify DCM object (%s)"},
    {SQ_LISTFAILURE, "SQ Function failed on list operation (%s)"},
    {SQ_NULLLIST, "SQ Null list passed to function %s"},
    {SQ_EMPTYLIST, "SQ Empty list passed to function %s"},
    {SQ_OBJECTACCESSFAILED, "SQ Failed to access DCM object in %s"},
    {SQ_INTERNALSEQBUILDFAILED,
    "SQ Failed to build internal sequence in %s"},
    {SQ_INTERNALSEQPARSEFAILED,
    "SQ Failed to parse internal sequence in %s"},
    {SQ_BADSEQUENCETYPEELEMENT, "SQ Bad sequence type element found in %s"},
    {0, NULL}
};


/* SQ_Message
**
** Purpose:
**	Find the ASCIZ message that goes with an DCM error number and
**	return a pointer to static memory containing that error message.
**
** Parameter Dictionary:
**	condition	The error condition number
**
** Return Values:
**	Appropriate error message corresponding to the condition. If the
**	condition is invalid, a NULL message is returned.
**
** Algorithm:
**	Description of the algorithm (optional) and any other notes.
*/

char *
SQ_Message(CONDITION condition)
{
    int
        index;

    for (index = 0; messageVector[index].message != NULL; index++)
	if (condition == messageVector[index].cond)
	    return messageVector[index].message;

    return NULL;
}
