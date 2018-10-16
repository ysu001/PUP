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
** Module Name(s):	UTL_RegexMatch, UTL_ConvertRegex,
**			UTL_ConvertDatetoLong, UTL_ConvertLongtoDate,
**			UTL_ConvertTimetoFloat, UTL_ConvertFloattoTime,
**			UTL_SqueezeBlanks, UTL_DateMatch, UTL_TimeMatch
**			UTL_GetDicomDate, UTL_GetDicomTime
**
** Author, Date:	David E. Beecher, March 1994
** Intent:		Miscellaneous functions that may be useful in
**			a number of different areas.
**
** Last Update:		$Author: karchie $, $Date: 2007/09/21 00:38:48 $
** Source File:		$RCSfile: utility.c,v $
** Revision:		$Revision: 1.3 $
** Status:		$State: Exp $
*/

/* The following system includes added for compatability */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <dirent.h>
#include <errno.h>
#include <signal.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#if HAVE_MALLOC
#include <malloc.h>
#endif
#include <fcntl.h>
#include <ctype.h>
#include <unistd.h>

#include "config.h"
#if HAVE_REGEX_H
#include <regex.h>
#endif

#include "include/dicom.h"
#include "include/lst.h"
#include "include/condition.h"
#include "include/utility.h"

#ifndef _MSC_VER
/* UTL_RegexMatch
**
** Purpose:
**	Perform a DICOM regular expression match with the specified string, stm.
**
** Parameter Dictionary:
**	char *regex:
**		The DICOM regular expression to try and match.
**	char *stm:
**		The input string to match.
**
** Return Values:
**	UTL_MATCH:	The input string matched the regular expression.
**	UTL_NOMATCH:	The input string did not match the regular expression.
**
** Algorithm:
**	A simple function to perform a DICOM regular expression match with the
**	specified  string, stm.  The sematics of the DICOM patterns must be altered
**	slightly to work correctly with regex under unix...more information may
**	be found below.
**
*/

CONDITION
UTL_RegexMatch(char *regex, char *stm)
{
#if (HAVE_REGEX_H && HAVE_REGCOMP)
    int
        ret,
        regexReturn;
    char
       *new_rstring;
    regex_t preg;
    char errorBuff[256];
    regmatch_t pmatch;

    new_rstring = UTL_ConvertRegex(regex);

    regexReturn = regcomp(&preg, new_rstring, 0);
    if (regexReturn != 0) {
	regerror(regexReturn, &preg, errorBuff, sizeof(errorBuff));
	fprintf(stderr, "%d\n", regexReturn);
	fprintf(stderr, "%s\n", errorBuff);

	free(new_rstring);
	return (UTL_NOMATCH);
    } else {
	ret = regexec(&preg, stm, 1, &pmatch, 0);

	switch (ret) {
	case 0:
	    free(new_rstring);
	    return (UTL_MATCH);
	    break;
	default:
	    free(new_rstring);
	    return (UTL_NOMATCH);
	    break;
	}
    }
#elif HAVE_RE_COMP
    int
        ret;
    char
       *new_rstring;

    new_rstring = UTL_ConvertRegex(regex);
    if (re_comp(new_rstring) != (char *) 0) {
	free(new_rstring);
	return (UTL_NOMATCH);
    } else {
	ret = re_exec(stm);
	switch (ret) {
	case 0:
	case -1:
	    free(new_rstring);
	    return (UTL_NOMATCH);
	    break;
	case 1:
	    free(new_rstring);
	    return (UTL_MATCH);
	    break;
	}
    }
#else
#error No suitable regular expression library found, sorry.
#endif
}

/* UTL_ConvertRegex
**
** Purpose:
**	This function converts a DICOM "regular expression" to the proper
**	regex semantics under unix.
**
** Parameter Dictionary:
**	char *regex:
**		The DICOM regular expression to convert.
**
** Return Values:
**	char *:	The converted regular expression which expresses DICOM pattern
**		matching in regex semantics.
**
** Notes:
**	This routine needs to return a string of unknown length.  Since we
**	don't want to burden the caller with having to remember to free the
**	string after it has been used, we simply reuse the same piece of storage
**	and realloc it when necessary to increase the size.
**
** Algorithm:
**	Simple function to convert a DICOM "regular expression" to the proper
**	regex semantics under unix.  DICOM has only 2 meta characters, "*" for 0
**	or more occurrences, and "?" for a single character.  The "*" must be
**	converted to ".*" for regex while the "?" must be converted to ".".
**	Other special characters to regex like "[", "]", and "." must also
**	be escaped with the "\".  The DICOM escape character is assumed to be "\".
*/
char *
UTL_ConvertRegex(char *regex)
{

    char
       *new_regex = (char *) NULL;
    int
        malloced_size = 0;
    int
        i,
        j,
        escape_on;

    if (new_regex == (char *) NULL) {
	malloced_size = REGEX_SIZE;
	if ((new_regex = (char *) malloc(malloced_size)) == (char *) NULL) {
	    return ((char *) NULL);
	}
    }
    i = j = 0;
    escape_on = OFF;
    new_regex[j++] = '^';
    while (regex[i] != '\000') {
	switch (regex[i]) {
	case '*':		/* Transform the "*" to ".*" or "\*" if
				 * escaped */
	    switch (escape_on) {
	    case OFF:
		new_regex[j++] = '.';
		break;
	    case ON:
		new_regex[j++] = '\\';
		escape_on = OFF;
		break;
	    }
	    new_regex[j++] = '*';
	    i++;
	    break;
	case '?':		/* Transform the "?" to "." or "?" if escaped */
	    switch (escape_on) {
	    case OFF:
		new_regex[j++] = '.';
		break;
	    case ON:
		new_regex[j++] = '?';
		escape_on = OFF;
		break;
	    }
	    i++;
	    break;
	case '\\':		/* Note that we have seen the escape
				 * character */
	    switch (escape_on) {
	    case OFF:
		escape_on = ON;
		break;
	    case ON:
		escape_on = OFF;
		new_regex[j++] = '\\';
		new_regex[j++] = '\\';
		break;
	    }
	    i++;
	    break;
	case '.':
	case '[':		/* These are special to regex and need to be
				 * escaped */
	case ']':
	    new_regex[j++] = '\\';
	    new_regex[j++] = regex[i++];
	    escape_on = OFF;
	    break;
	default:		/* Leave the "\" in at this juncture */
	    switch (escape_on) {
	    case ON:
		new_regex[j++] = '\\';
		escape_on = OFF;
		break;
	    case OFF:
		break;
	    }
	    new_regex[j++] = regex[i++];
	    break;
	}
	if (j >= (malloced_size - 2)) {
	    malloced_size += REGEX_SIZE;
	    if ((new_regex = (char *) realloc(new_regex, malloced_size)) ==
		(char *) NULL) {
		return ((char *) NULL);
	    }
	}
    }
    new_regex[j++] = '$';
    new_regex[j] = '\000';
    return (new_regex);
}
#endif
/* UTL_ConvertDatetoLong
**	Convert a Dicom date to a long for comparision ease.
*/
long
UTL_ConvertDatetoLong(const char *date)
{

    char
        year[5],
        month[3],
        day[3];

    strncpy(year, date, 4);
    year[4] = '\000';
    strncpy(month, date + 4, 2);
    month[2] = '\000';
    strncpy(day, date + 6, 2);
    day[2] = '\000';

    return ((atol(year) * 10000) + (atol(month) * 100) + atol(day));
}

/* UTL_ConvertLongtoDate
**	Convert a long to a Dicom date.
*/
void
UTL_ConvertLongtoDate(long ld, char *date)
{

    int
        year,
        month,
        day;

    year = ld / 10000;
    ld -= (year * 10000);
    month = ld / 100;
    ld -= (month * 100);
    day = ld;

    sprintf(date, "%04d%02d%02d", year, month, day);

    return;
}

/* UTL_ConvertTimetoFloat
**	Convert a Dicom time to a floating point number for comparision ease.
*/
double
UTL_ConvertTimetoFloat(const char *time)
{

    size_t
    i;
    char
        hour[3],
        minute[3],
        second[3],
        fracsec[7];
    const char *p;
    double
        divisor,
        hh,
        mm,
        ss,
        fs;

    hh = mm = ss = fs = 0.0;
    hour[0] = minute[0] = second[0] = fracsec[0] = '\000';

    p = time;
    /*
     * Just a brute force way to tear down a Dicom time...not very pretty,
     * but it works... We are not guaranteed to have every field present as
     * we are in the date...
     */
    hour[0] = *p++;
    hour[1] = *p++;
    hour[2] = '\000';
    if (isdigit(*p)) {
	minute[0] = *p++;
	minute[1] = *p++;
	minute[2] = '\000';
	if (isdigit(*p)) {
	    second[0] = *p++;
	    second[1] = *p++;
	    second[2] = '\000';
	    if (*p == '.') {
		p++;
		fracsec[0] = *p++;
		if ((*p != '\000') && (isdigit(*p))) {
		    fracsec[1] = *p++;
		    if ((*p != '\000') && (isdigit(*p))) {
			fracsec[2] = *p++;
			if ((*p != '\000') && (isdigit(*p))) {
			    fracsec[3] = *p++;
			    if ((*p != '\000') && (isdigit(*p))) {
				fracsec[4] = *p++;
				if ((*p != '\000') && (isdigit(*p))) {
				    fracsec[5] = *p++;
				    fracsec[6] = '\000';
				} else
				    fracsec[5] = '\000';
			    } else
				fracsec[4] = '\000';
			} else
			    fracsec[3] = '\000';
		    } else
			fracsec[2] = '\000';
		} else
		    fracsec[1] = '\000';
	    }
	}
    }
    hh = atof(hour);
    mm = atof(minute);
    ss = atof(second);
    divisor = 1;
    for (i = 0; i < strlen(fracsec); i++)
	divisor *= 10;
    fs = atof(fracsec) / divisor;

    return ((hh * 3600.0) + (mm * 60.0) + ss + fs);
}

/* UTL_ConvertFloattoTime
**	Convert a floating point number to a Dicom time.
*/
void
UTL_ConvertFloattoTime(double dt, char *time)
{
    int
        hour,
        minute,
        second,
        fracsec;

    hour = (int) (dt / 3600.0);
    dt -= (hour * 3600);

    minute = (int) (dt / 60.);
    dt -= (minute * 60);

    second = (int) dt;
    dt -= second;

    fracsec = (int) ((dt * 1000000) + 0.5);

    sprintf(time, "%02d%02d%02d.%06d", hour, minute, second, fracsec);

    return;
}


/* UTL_SqueezeBlanks
**
*/
void
UTL_SqueezeBlanks(char *s)
{

    char
       *t1,
       *t2;

    t1 = t2 = s;
    while (*t2 != '\000') {
	if (*t2 != ' ') {
	    *t1 = *t2;
	    t1++;
	}
	t2++;
    }
    *t1 = '\000';

    return;
}
/* UTL_DateMatch
**	Match a date range as specified in the Dicom standard
*/
CONDITION
UTL_DateMatch(char *datestring, char *stm)
{

    int
        match;
    char
       *ndate;
    long
        start_date,
        end_date,
        date_in_question;

    if ((ndate = (char *) malloc(strlen(datestring) + 1)) == (char *) NULL)
	return (UTL_NOMATCH);

    strcpy(ndate, datestring);
    UTL_SqueezeBlanks(ndate);
    UTL_SqueezeBlanks(stm);

    match = 0;
    if (strchr(ndate, (int) '-') == (char *) NULL) {
	if (strcmp(ndate, stm) == 0)
	    match = 1;
    } else {
	date_in_question = UTL_ConvertDatetoLong(stm);
	if (ndate[0] == '-') {
	    end_date = UTL_ConvertDatetoLong(ndate + 1);
	    if (date_in_question <= end_date)
		match = 1;
	} else if (ndate[strlen(ndate) - 1] == '-') {
	    start_date = UTL_ConvertDatetoLong(ndate);
	    if (date_in_question >= start_date)
		match = 1;
	} else {
	    start_date = UTL_ConvertDatetoLong(ndate);
	    end_date = UTL_ConvertDatetoLong(strchr(ndate, (int) '-') + 1);
	    if ((date_in_question >= start_date) &&
		(date_in_question <= end_date))
		match = 1;
	}
    }
    free(ndate);
    if (match)
	return (UTL_MATCH);
    else
	return (UTL_NOMATCH);
}
/* UTL_TimeMatch
**	Match a time range as specified in the Dicom standard
*/
CONDITION
UTL_TimeMatch(char *timestring, char *stm)
{

    int
        match;
    char
       *ntime;
    double
        start_time,
        end_time,
        time_in_question;

    if ((ntime = (char *) malloc(strlen(timestring) + 2)) == (char *) NULL)
	return (UTL_NOMATCH);

    strcpy(ntime, timestring);
    UTL_SqueezeBlanks(ntime);
    UTL_SqueezeBlanks(stm);

    match = 0;
    if (strchr(ntime, (int) '-') == (char *) NULL) {
	if (strcmp(ntime, stm) == 0)
	    match = 1;
    } else {
	time_in_question = UTL_ConvertTimetoFloat(stm);
	if (ntime[0] == '-') {
	    end_time = UTL_ConvertTimetoFloat(ntime + 1);
	    if (time_in_question <= end_time)
		match = 1;
	} else if (ntime[strlen(ntime) - 1] == '-') {
	    start_time = UTL_ConvertTimetoFloat(ntime);
	    if (time_in_question >= start_time)
		match = 1;
	} else {
	    start_time = UTL_ConvertTimetoFloat(ntime);
	    end_time = UTL_ConvertTimetoFloat(strchr(ntime, (int) '-') + 1);
	    if ((time_in_question >= start_time) &&
		(time_in_question <= end_time))
		match = 1;
	}
    }
    free(ntime);
    if (match)
	return (UTL_MATCH);
    else
	return (UTL_NOMATCH);
}
/*
** UTL_GetDicomDate
**	Get the current date and store as a Dicom date.
*/
void
UTL_GetDicomDate(char *datestr)
{

    struct tm
       *tf;
    time_t
	loctime;

    loctime = time((time_t *) NULL);
    tf = localtime(&loctime);

    sprintf(datestr, "%04d%02d%02d", (tf->tm_year) + 1900, (tf->tm_mon) + 1, tf->tm_mday);
    return;

}
/*
** UTL_GetDicomTime
**	Get the current time and store as a Dicom time.
*/
void
UTL_GetDicomTime(char *timestr)
{

    struct tm
       *tf;
    time_t
	loctime;

    loctime = time((time_t *) NULL);
    tf = localtime(&loctime);

    sprintf(timestr, "%02d%02d%02d.%06d", (tf->tm_hour), (tf->tm_min), (tf->tm_sec), 0);
    return;
}

#ifdef _MSC_VER
typedef struct {
    char key[10];
    struct _timeb t;
}   UTL_TIMESTRUCTURE;
#else
typedef struct {
    char key[10];
    struct timeval t;
}   UTL_TIMESTRUCTURE;
#endif

void *
UTL_GetTimeStamp()
{
    UTL_TIMESTRUCTURE *t;

    t = calloc(1, sizeof(*t));
    if (t == NULL)
	return NULL;

    strcpy(t->key, "UTL STAMP");

#ifdef _MSC_VER
    _ftime(&t->t);
#elif (TIMEOFDAYARGS == 2)
    gettimeofday(&t->t, NULL);
#else
    /* Correction ? */
    gettimeofday(&t->t, NULL);
#endif

    return t;
}

double
UTL_DeltaTime(void *timeStamp)
{
#ifdef _MSC_VER
    struct _timeb timeNow;
#else
    struct timeval timeNow;
#endif
    UTL_TIMESTRUCTURE *t;
    double delta = 0.;

#ifdef _MSC_VER
    _ftime(&timeNow);
#elif (TIMEOFDAYARGS == 2)
    gettimeofday(&timeNow, NULL);
#else
    /* Correction? */
    gettimeofday(&timeNow, NULL);
#endif

    t = (UTL_TIMESTRUCTURE *) timeStamp;
    if (t == NULL)
	return -1.0;

    if (strcmp(t->key, "UTL STAMP") != 0)
	return -1.0;

#ifdef _MSC_VER
    delta = timeNow.time - t->t.time;
    delta += (timeNow.millitm - t->t.millitm) / 1000.;
#else
    delta = timeNow.tv_sec - t->t.tv_sec;
    delta += (timeNow.tv_usec - t->t.tv_usec) / 1000000.;
#endif

    return delta;
}

void
UTL_ReleaseTimeStamp(void *timeStamp)
{
    UTL_TIMESTRUCTURE *t;

    t = (UTL_TIMESTRUCTURE *) timeStamp;
    if (t == NULL)
	return;

    if (strcmp(t->key, "UTL STAMP") != 0)
	return;

    free(timeStamp);
}

CONDITION
UTL_VerifyCreatePath(const char *path)
{
    int i;
#ifdef _MSC_VER
    struct _stat buf;
#else
    struct stat buf;
#endif
    char
       *p,
        temp[1024];
    int flag = 0;
    static int statCount = 0;

#ifdef _MSC_VER
    statCount++;
    i = _stat(path, &buf);
#else
    i = stat(path, &buf);
#endif


    if (i == 0) {
#ifdef _MSC_VER
	flag = ((buf.st_mode & _S_IFDIR) != 0);
#else
	flag = (S_ISDIR(buf.st_mode));
#endif
	if (flag)
	    return UTL_NORMAL;
	else
	    return UTL_PATHNOTDIR;
    }
    p = temp;

    while (*path != '\0') {
	*p++ = *path++;
	while (*path != '/' && *path != '\\' && *path != '\0') {
#ifdef _MSC_VER
	    if (*path == ':') {
		*p++ = *path++;
		if (*path == '\0')	/* We should not get C:\0, but test
					 * it */
		    break;
	    }
#endif
	    *p++ = *path++;
	}

	*p = '\0';
#ifdef _MSC_VER
	statCount++;
	i = _stat(temp, &buf);
#else
	i = stat(temp, &buf);
#endif

	if (i == 0) {
#ifdef _MSC_VER
	    flag = ((buf.st_mode & _S_IFDIR) != 0);
#else
	    flag = (S_ISDIR(buf.st_mode));
#endif
	    if (!flag)
		return UTL_PATHNOTDIR;
	} else {
#ifdef _MSC_VER
	    int e1;
	    e1 = errno;
	    memset(&buf, 0, sizeof(buf));
	    /*fprintf(stderr, "Stat Count = %d\n", statCount);*/
	    statCount++;
	    i = _stat(temp, &buf);
	    e1 = errno;
	    i = _mkdir(temp);
#else
	    i = mkdir(temp, 0777);
#endif
	    if (i != 0) {
		int e1;
		e1 = errno;
		fprintf(stderr, "Stat Count = %d\n", statCount);
		perror(temp);
		return UTL_FILECREATEFAILED;
	    }
	}
    }
    return UTL_NORMAL;
}

CTNBOOLEAN UTL_IsDirectory(const char* path)
{
    int i;
#ifdef _MSC_VER
    struct _stat buf;
#else
    struct stat buf;
#endif

    int flag = 0;

#ifdef _MSC_VER
    i = _stat(path, &buf);
#else
    i = stat(path, &buf);
#endif


    if (i == 0) {
#ifdef _MSC_VER
	flag = ((buf.st_mode & _S_IFDIR) != 0);
#else
	flag = (S_ISDIR(buf.st_mode));
#endif
	if (flag)
	    return TRUE;
    }
    return FALSE;
}



CONDITION UTL_ScanDirectory(const char* path,
			    LST_HEAD** lst)
{
  UTL_FILEITEM* item = 0;

#ifdef _WIN32
  long hFile = 0;
  struct _finddata_t fileInfo;
  char directoryText[1024];
  *lst = LST_Create();
  strcpy(directoryText, path);
  strcat(directoryText, "/*");
  if( (hFile = _findfirst(directoryText, &fileInfo)) == -1L)
    return 0;

  item = malloc(sizeof(*item));
  strcpy(item->path, fileInfo.name);
  LST_Enqueue(lst, item);

  while(_findnext(hFile, &fileInfo) == 0) {
    item = malloc(sizeof(*item));
    strcpy(item->path, fileInfo.name);
    LST_Enqueue(lst, item);
  }
  _findclose(hFile);

#else
  DIR* dirp;
  struct dirent* dp;

  *lst = LST_Create();
  dirp = opendir(path);
  if (dirp == 0)
    return 0;

  while ((dp = readdir(dirp)) != NULL) {
    item = malloc(sizeof(*item));
    strcpy(item->path, dp->d_name);
    LST_Enqueue(lst, item);
  }
  closedir(dirp);
#endif

  return UTL_NORMAL;
}

static char* UTL_configFile = 0;
static LST_HEAD* UTL_configList = 0;
typedef struct {
  void* reserved[2];
  char *pName;
  char *pValue;
} CONFIG_ITEM;

CONDITION UTL_ReadConfigFile( )
{
  FILE* f;
  char buf[1024];

  if (UTL_configList != 0)
    return UTL_NORMAL;

  UTL_configList = LST_Create();
  if (UTL_configList == NULL)
    return 0;

  if (UTL_configFile == 0)
    return UTL_NORMAL;

  if (UTL_configFile[0] == '\0')
    return UTL_NORMAL;

  f = fopen(UTL_configFile, "r");
  if (f == NULL)
    return 0;

  while (fgets(buf, sizeof(buf), f) != NULL) {
    char* token1;
    char* token2;
    CONFIG_ITEM* item;

    if (buf[0] == '#') continue;
    if (buf[0] == '\n') continue;
    token1 = strtok(buf, " \t\n");
    token2 = strtok(0, " \t\n");
    if (token2 == NULL) continue;

    item = (CONFIG_ITEM*)malloc(sizeof(*item) + strlen(token1) +
				strlen(token2) + 2);
    item->pName = ((char*)item) + sizeof(*item);
    strcpy(item->pName, token1);
    item->pValue = item->pName + strlen(token1) + 1;
    strcpy(item->pValue, token2);

    LST_Enqueue(&UTL_configList, item);
  }

  fclose(f);

  return UTL_NORMAL;
}

CONDITION UTL_SetConfigFile(const char* configFile)
{
  if (UTL_configFile != 0) {
    CTN_FREE(UTL_configFile);
  }

  if (configFile == 0 || configFile[0] == '\0') {
    char* p = getenv("CTN_TARGET");
    if (p == NULL) {
      return UTL_NO_CTN_TARGET;
    }
    UTL_configFile = (char*) malloc(strlen(p) + strlen("/runtime/ctn_cfg.txt") + 1);
    strcpy(UTL_configFile, p);
    strcat(UTL_configFile, "/runtime/ctn_cfg.txt");
  } else {
    UTL_configFile = (char*) malloc(strlen(configFile)+1);
    strcpy(UTL_configFile, configFile);
  }

  return UTL_NORMAL;
}

CONDITION UTL_TestConfigFile(const char* configFile)
{
  return UTL_NORMAL;
}
char* UTL_GetConfigParameter(const char* paramName)
{
  CONDITION cond;
  char nameCopy[256];
  CONFIG_ITEM* item;
  int idx;

  cond = UTL_ReadConfigFile( );
  if (cond != UTL_NORMAL)
    return NULL;

  item = LST_Head(&UTL_configList);
  if (item == NULL)
    return NULL;

  (void) LST_Position(&UTL_configList, item);
  while(item != NULL) {
    if (strcmp(item->pName, paramName) == 0)
      return item->pValue;

    item = LST_Next(&UTL_configList);
  }

  strcpy(nameCopy, paramName);
  idx = strlen(nameCopy) - 1;
  while (idx > 0) {
    if (nameCopy[idx] == '/') {
      nameCopy[idx] = '\0';
      idx = -1;
      break;
    } else {
      idx--;
    }
  }

  if (idx < 0) {
    return UTL_GetConfigParameter(nameCopy);
  } else {
    return NULL;
  }
}

char**
UTL_ExpandToPointerArray(const char* inputText,
			 const char* delimiters,
			 int* numberOfEntries)
{
  int idx;
  int memorySize = 0;
  int arrayIndex = 0;
  char** array;
  char* outputPtr;
  char* token;

  *numberOfEntries = 1;
  for (idx = 0; inputText[idx] != '\0'; idx++) {
    int j;
    for (j = 0; delimiters[j] != '\0'; j++) {
      if (inputText[idx] == delimiters[j]) {
	(*numberOfEntries)++;
	break;
      }
    }
  }

  memorySize = (sizeof(char*)) * (*numberOfEntries);
  memorySize += strlen(inputText) + 1;

  array = (char**)CTN_MALLOC(memorySize);
  outputPtr = ((char*) array) + ((sizeof(char*)) * (*numberOfEntries));
  strcpy(outputPtr, inputText);

  token = strtok(outputPtr, delimiters);
  while(token != NULL) {
    array[arrayIndex++] = token;
    token = strtok(NULL, delimiters);
  }

  return array;
}

CTNBOOLEAN UTL_IsFile(const char* path)
{
  int i;
  CTNBOOLEAN rtnValue = FALSE;

#ifdef _WIN32
  struct _stat buf;

  i = _stat(path, &buf);
  if (i == 0) {
    rtnValue = ((buf.st_mode & _S_IFREG) != 0);
  }
#else
  struct stat buf;
  i = stat(path, &buf);
  if (i == 0) {
    rtnValue = (S_ISREG(buf.st_mode));
  }
#endif

  return rtnValue;
}

CONDITION UTL_DeleteFile(const char* path)
{
  int i = 0;

  i = unlink(path);

  if (i == 0)
    return UTL_NORMAL;

  return COND_PushCondition(UTL_DELETEFILEFAILED, "");
}


CONDITION
UTL_FileSize(const char* path, U32* size)
{
  int status;
  struct stat im_stat;

  status = stat(path, &im_stat);
  if (status < 0) {
    *size = 0;
    return 0;
  } else {
    *size = im_stat.st_size;
    return UTL_NORMAL;
  }
}

