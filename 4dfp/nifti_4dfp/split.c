/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/split.c,v 1.1 2009/07/08 01:01:35 coalsont Exp $*/
/*$Log: split.c,v $
 * Revision 1.1  2009/07/08  01:01:35  coalsont
 * Initial revision
 **/
static char* rcsid = "$Id: split.c,v 1.1 2009/07/08 01:01:35 coalsont Exp $";
/* Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009 */
/* borrowed from TRX/asciito4dfp.c */
/*
 *  split.c
 *  
 *
 *  Created by NRG on 6/25/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;
	
	if ((ptr = strchr (string, '#'))) *ptr = '\0';
	i = m = 0;
	while (m < maxp) {
		while (!isgraph ((int) string[i]) && string[i]) i++;
		if (!string[i]) break;
		srgv[m++] = string + i;
		while (isgraph ((int) string[i])) i++;
		if (!string[i]) break;
		string[i++] = '\0';
	}
	return m;
}

