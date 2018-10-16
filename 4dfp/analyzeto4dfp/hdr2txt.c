/*$Header: /data/petsun4/data1/src_solaris/analyzeto4dfp/RCS/hdr2txt.c,v 1.3 2007/09/21 23:36:38 avi Exp $*/
/*$Log: hdr2txt.c,v $
 * Revision 1.3  2007/09/21  23:36:38  avi
 * linux, X86 compliant
 *
 * Revision 1.2  2002/06/24  18:55:37  chad
 * do not strip off .4dfp .4dint
 *
 * Revision 1.1  2002/06/20  01:49:19  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <endianio.h>

#define MAXL 256

void getrootl (char *filespc, char *imgroot) {
	char	*str;
	strcpy (imgroot, filespc);
	while (str = strrchr (imgroot, '.')) {
			if (!strcmp (str, ".rec"))	*str = '\0';
		else	if (!strcmp (str, ".img"))	*str = '\0';
		else	if (!strcmp (str, ".ifh"))	*str = '\0';
		else	if (!strcmp (str, ".hdr"))	*str = '\0';
		else	break;
	}
}

static char rcsid[] = "$Id: hdr2txt.c,v 1.3 2007/09/21 23:36:38 avi Exp $";
int main (int argc, char *argv[]) {
/*********/
/* flags */
/*********/
	int			debug = 0;
	int			swab_flag = 0;
	int			status = 0;

	FILE			*anafp;
	char			*ptr, command[MAXL], program[MAXL];
	char			imgroot[MAXL], hdrfile[MAXL], outfile[MAXL];
	float			mmppix[3];
	int			xdim, ydim, zdim;
	int			c, i, j, k;
	struct dsr		header;			/* ANALYZE header */

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);
/******************************/
/* get command line arguments */
/******************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
			}
		} else switch (k) {
		 	case 0: getrootl (argv[i], imgroot);	k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage: %s <analyze_image>\n", program);
		printf (" e.g., %s brain_asig[.hdr]]\n", program);
		printf ("\toption\n");
		exit (1);
	}

/***********************/
/* read ANALYZE header */
/***********************/
	sprintf (hdrfile, "%s.hdr", imgroot);
	if (!(anafp = fopen (hdrfile, "rb"))) errr (program, hdrfile);
	printf ("Reading: %s\n", hdrfile);
	k = fread (&header, sizeof (header), 1, anafp);
	if (k != 1) errr (program, hdrfile);
	fclose (anafp);
	if (header.hk.sizeof_hdr != sizeof (struct dsr)) {
		printf ("converting byte order\n");
		swab_hdr (&header);
		swab_flag++;
	}
	printf ("header size=%d\n", header.hk.sizeof_hdr);
	printf ("%s\n", header.hk.db_name);
	printf ("%s\n", header.hist.descrip);
	strcpy (command, header.hist.originator);
	if ((ptr = strchr (command, ' '))) *ptr = '\0';
	printf ("originator=%s\n", command);

	sprintf (outfile, "%s.hdr.txt", imgroot);
	if (!(anafp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);
	xdim = header.dime.dim[1];
	ydim = header.dime.dim[2];
	zdim = header.dime.dim[3];
	fprintf (anafp, "%10d%10d%10d\n", xdim, ydim, zdim);
	for (k = 0; k < 3; k++) mmppix[k] = header.dime.pixdim[k + 1];
	fprintf (anafp, "%10.4f%10.4f%10.4f\n", mmppix[0], mmppix[1], mmppix[2]);
	fclose (anafp);

	sprintf (command, "cat %s\n", outfile); system (command);
        exit (status);
}
