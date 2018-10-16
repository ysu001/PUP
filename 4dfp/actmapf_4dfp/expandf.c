/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/expandf.c,v 1.7 2009/07/30 05:40:54 avi Exp $*/
/*$Log: expandf.c,v $
 * Revision 1.7  2009/07/30  05:40:54  avi
 * extensive checking to prevent overwriting memory
 *
 * Revision 1.6  2007/08/14  01:24:55  avi
 * more informative error message
 *
 * Revision 1.5  2006/10/02  01:49:26  avi
 * complete Solaris 10 #includes
 *
 * Revision 1.4  2005/08/30  21:51:48  avi
 * better error reporting
 *
 * Revision 1.3  1999/03/04  03:52:53  avi
 * minor cleaning
 *
 * Revision 1.2  1997/12/05  04:48:05  avi
 * parse cos and sin epochs
 *
 * Revision 1.1  1997/04/28  00:51:47  yang
 * Initial revision
 **/

#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static	char	rcsid[] = "$Id: expandf.c,v 1.7 2009/07/30 05:40:54 avi Exp $";
int expandf (char *strin, int len) {
	char	*stringt, *string;	/* buffers */
	char	*lp;			/* left delimeter */
	char	*np;			/* ascii digit */
	char	*rp;			/* right delimiter */
	char	ascnum[16];
	char	c2[2];
	int	c, i, k, l;
	int	len2;
	int	level;
	int	num;
	int	status;
	int	debug = 0;

	if (!strlen (strin)) return 0;
	status = 0;
	memset (c2, '\0', 2);
	len2 = len * 2;
	stringt = (char *) calloc (len2 + 1, sizeof (char));
	string  = (char *) calloc (len2 + 1, sizeof (char));
	if (!stringt || !string) {
		fprintf (stderr, "expandf: memory allocation error\n");
		exit (-1);
	}
	if (debug) printf ("expandf: len=%d\n", len);

/***************************/
/* squeeze out white space */
/***************************/
	rp = strin;
	lp = string;
	while (c = *rp++) if (!isspace (c)) *lp++ = c;
	*lp = '\0';	/* safety */
	if (debug) {l = strlen (string); printf ("%s\t%d\n", string, l);}

/****************************/
/* expand sinusoidal epochs */
/****************************/
	rp = stringt;			/* expansion buffer */
	lp = string;			/* input string */
	stringt[len2] = '\xff';		/* write stop barrier */
	while (c = *lp++) {
		switch (c) {
		case 'C': case 'c': case 'S': case 's':
			np = ascnum;
			k = 0;
			while (k < 16 && isdigit (*lp)) {
				*np++ = *lp++;
				k++;
			}
			*np = '\0';
			k = atoi (ascnum) - 1;
			*rp++ = '(';			if (*rp == '\xff') goto ABORT;
			for (i = 0; i < k; i++) {
				*rp++ = c;		if (*rp == '\xff') goto ABORT;
			}
			*rp++ = '~';			if (*rp == '\xff') goto ABORT;
			*rp++ = ')';			if (*rp == '\xff') goto ABORT;
			break;
		default:
			*rp++ = c;			if (*rp == '\xff') goto ABORT;
			break;
		}
	}
	*lp = '\0';	/* safety */
	stringt[len2] = '\0';				/* remove write stop barrier */
	strncpy (string, stringt, len2);
	if (debug) {l = strlen (string); printf ("%s\t%d\n", string, l);}

/**********************/
/* expand parentheses */
/**********************/
	level = 0;
	while (rp = strrchr (string, ')')) {
		*rp = '\0';
		lp = rp;
		while (lp > string && isdigit (c = *(lp - 1))) *--lp = '\0';
		level++;
		while ((level > 0) && (lp > string)) {
			lp--;
			if (*lp == ')') level++;
			if (*lp == '(') level--;
		}
		if (level) {
			printf ("expandf error: unbalanced parentheses\n");
			status = 1; goto DONE;
		}
		*lp = '\0';
		num = 1;
		np = lp;
		while (np > string && isdigit (c = *(np - 1))) np--;
		if (strlen (np) > 0) {
			num = atoi (np);		/* printf ("num=%d\n", num); */
			*np = '\0';
		}
		strncpy (stringt, string, len2);
		if (strlen (stringt) + num*strlen (lp + 1) > len2) goto ABORT;
		for (k = 0; k < num; k++) strcat (stringt, lp + 1);
		if (strlen (stringt) + strlen (rp + 1) > len2) goto ABORT;
		strcat (stringt, rp + 1);
		strncpy (string, stringt, len2);	/* printf ("%s\n", string); */
	}
	if (strrchr (string, '(')) {
		printf ("expandf error: unbalanced parentheses\n");
		status = 1; goto DONE;
	}
	if (debug) {l = strlen (string); printf ("%s\t%d\n", string, l);}

/********************/
/* expand multiples */
/********************/
	while (np = strpbrk (string, "0123456789")) {
		rp = np;
		while (isdigit (c = *++rp));
		strncpy (c2, rp, 1);
		*rp = '\0';
		num = atoi (np);			/* printf ("num=%d\n", num); */
		*np = '\0';
		strncpy (stringt, string, len2);
		if (strlen (stringt) + num > len2) goto ABORT;
		for (k = 0; k < num; k++) strcat (stringt, c2);
		if (strlen (stringt) + strlen (rp + 1) > len2) goto ABORT;
		strcat (stringt, rp + 1);
		strncpy (string, stringt, len2);	/* printf ("%s\n", string); */
	}
	if (debug) {l = strlen (string); printf ("%s\t%d\n", string, l);}

DONE:	l = strlen (string);
	if (l >= len) {
		fprintf (stderr, "expandf error: expanded format length exceeds allocated buffer size (%d)\n", len);
		status = 1;
	}
	string[len - 1] = '\0';
	if (!status) strcpy (strin, string);
	free (string);
	free (stringt);
	return (status);
ABORT:	status = -1;
	goto DONE;
}
