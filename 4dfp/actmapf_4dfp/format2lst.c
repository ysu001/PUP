/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/format2lst.c,v 1.6 2009/07/30 05:42:08 avi Exp $*/
/*$Log: format2lst.c,v $
 * Revision 1.6  2009/07/30  05:42:08  avi
 * linux compliant
 *
 * Revision 1.5  2005/08/30  21:52:22  avi
 * expanded format buffer size increased to 16384 chars
 *
 * Revision 1.4  2005/08/30  07:05:12  avi
 * -e option
 *
 * Revision 1.3  2004/05/16  01:11:54  avi
 * -w option
 *
 * Revision 1.2  1999/04/17  05:22:31  avi
 * correct MAXL (256)
 *
 * Revision 1.1  1999/04/16  05:15:19  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAXF	16384
#define	MAXL	256

extern int	expandf (char *string, int len);

static char rcsid[] = "$Id: format2lst.c,v 1.6 2009/07/30 05:42:08 avi Exp $";
int main (int argc, char *argv[]) {
	char		format[MAXF];

/***********/
/* utility */
/***********/
	char		command[MAXL], *ptr, program[MAXL];
	int		c, i, j, k, n;

/*********/
/* flags */
/*********/
	int		expand = 0;
	int		wei_flag = 0;
	int		debug = 0;

	if (ptr = strrchr (argv[0], '/')) ptr++; else ptr = argv[0];
	strcpy (program, ptr);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'e': expand++;		break;
				case 'w': wei_flag++;		break;
				case 'd': debug++;		break;
			}
		}
		else switch (k) {
			case 0:	strcpy (format, argv[i]);	k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:\t%s <format>\n", program);
		printf ("e.g.,\t%s \"2x3-2+1-2+2-2+1-1+2-1+1-1+1-1+2-1+1-1+1-2+2-1+1-1+2-2+1-2+2-\"\n", program);
		printf ("\toptions\n");
		printf ("\t-w\tconvert {'x' '+' '-'} to {0.0 1.0 -1.0}\n");
		printf ("\t-e\texpand on single line\n");
		exit (1);
	}
	n = strlen (format);
	if (expandf (format, MAXF)) exit (-1);
	n = strlen (format);
	for (i = 0; i < n; i++) {
		if (wei_flag) {
			switch (format[i]) {
				case '+':	printf ("1.0");		break;
				case '-': 	printf ("-1.0");	break;
				default:	printf ("0.0");		break;
			}
		} else {
			printf ("%c", format[i]);
		}
		if (!expand) printf ("\n");
	}
	exit (0);
}
