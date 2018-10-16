/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4scale.c,v 1.2 2007/05/04 01:33:55 avi Exp $*/
/*$Log: t4scale.c,v $
 * Revision 1.2  2007/05/04  01:33:55  avi
 * supply missing argument to error report
 *
 * Revision 1.1  2006/09/26  22:50:47  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#define MAXL		256

float t4scale (char *t4file) {
	FILE		*fp;
	float		q;
	char		*ptr, string[MAXL];
	
	q = 1.0;
	if (!(fp = fopen (t4file, "r"))) {
		fprintf (stderr, "t4scale: %s read error\n", t4file);
		exit (-1);
	}
	while (fgets (string, MAXL, fp)) if ((ptr = strstr (string, "scale:"))) q = atof (ptr + 6);
	fclose (fp);
	return q;
}
