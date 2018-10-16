/*$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/conc_test.c,v 1.4 2006/09/23 05:08:29 avi Exp avi $*/
/*$Log: conc_test.c,v $
 * Revision 1.4  2006/09/23  05:08:29  avi
 * make Solaris 10 compliant
 *
 * Revision 1.3  2006/09/23  04:25:02  avi
 * endian control
 *
 * Revision 1.2  2004/12/28  03:37:04  avi
 * Revision 1.1  2004/11/27  05:45:05  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <rec.h>
#include <conc.h>
#include <Getifh.h>
#include <endianio.h>

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: conc_test.c,v 1.4 2006/09/23 05:08:29 avi Exp avi $";
static char	program[MAXL];
static char	trailer[MAXL] = "xxxx";

void usage (char *program) {
	fprintf (stderr, "Usage:\t%s <concfile>\n", program);
	fprintf (stderr, "e.g.,\t%s test.conc\n", program);
	fprintf (stderr, "\t-a<str>\tspecify output concfile trailer (default = \"%s\")\n", trailer);
	fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default CPU-endian)\n");
	exit (1);
}

int main (int argc, char *argv[]) {
	char		*ptr, command[MAXL];
	CONC_BLOCK	conc_block;
	char		control = 'c', concfile[MAXL] = "";
	int		status, c, i, k, l;
	float		*imgt;
	int		ifile, nfile, imgdim[4], vdim, ivol;

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'a': strcpy (trailer, ptr);	*ptr = '\0'; break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	strcpy (concfile, argv[i]);	k++; break;
		}	
	}
	if (k < 1) usage (program);

/* verify that calling fclose is NOT OK on an uninitialized *FILE
	printf ("OK0\n");
	fclose (conc_block.imgfp);
	fclose (conc_block.imgfp);
*/
	printf ("OK\n");
	conc_init_quiet (&conc_block, program);
	conc_open_quiet (&conc_block, concfile);
if (0) {
	conc_free       (&conc_block);
	printf ("OK1\n");
	conc_init (&conc_block, program);
	conc_open (&conc_block, concfile);
}
	printf ("OK1\n");

	nfile = conc_block.rnfile;
	for (ifile = 0; ifile < nfile; ifile++) {
		imgdim[3] = conc_block.nvol[ifile];
		printf ("nvol[%d]=%d\n", ifile, conc_block.nvol[ifile]);
		printf ("imgfile0[%d]=%s\n", ifile, conc_block.imgfile0[ifile]);
	}
	if (!(imgt = (float *) calloc (conc_block.vdim, sizeof (float)))) errm (program);
	vdim = conc_block.vdim;

/**********************/
/* conc_rewind() test */
/**********************/
	vdim = conc_block.vdim;
	printf ("Reading first volume\n");
	conc_read_vol (&conc_block, imgt);
	for (i = 0; i < 10; i++) printf ("%10.4f", imgt[i*(vdim/10)]); printf ("\n");
	conc_rewind (&conc_block);
	printf ("Reading first volume\n");
	conc_read_vol (&conc_block, imgt);
	for (i = 0; i < 10; i++) printf ("%10.4f", imgt[i*(vdim/10)]); printf ("\n");
	conc_rewind (&conc_block);

if (1) {
/************************/
/* conc_read_vol() test */
/************************/
	printf ("Reading %s volume:", conc_block.lstfile);
	for (l = 0; l < conc_block.imgdim[3]; l++) {printf (" %d", l + 1); fflush (stdout);
		conc_read_vol (&conc_block, imgt);
	} printf ("\n"); fflush (stdout);
	printf ("rifile=%d rivol=%d\n", conc_block.rifile, conc_block.rivol);
}

/********************/
/* conc_seek() test */
/********************/
	ivol = 151;
	conc_seek (&conc_block, ivol);
	printf ("rifile=%d rivol=%d\n", conc_block.rifile, conc_block.rivol);
	printf ("Reading volume %d\n", ivol);
	conc_read_vol (&conc_block, imgt); ivol++;
	for (i = 0; i < 10; i++) printf ("%10.4f", imgt[i*(vdim/10)]); printf ("\n");
	printf ("Reading volume %d\n", ivol);
	conc_read_vol (&conc_block, imgt); ivol++;
	for (i = 0; i < 10; i++) printf ("%10.4f", imgt[i*(vdim/10)]); printf ("\n");
	ivol = 0;
	conc_seek (&conc_block, ivol);
	printf ("Reading volume %d\n", ivol);
	conc_read_vol (&conc_block, imgt); ivol++;
	for (i = 0; i < 10; i++) printf ("%10.4f", imgt[i*(vdim/10)]); printf ("\n");

if (0) {
/*******************/
/* conc_new() test */
/*******************/
	conc_rewind (&conc_block);
	conc_newe (&conc_block, trailer, control);
	printf ("Reading: %s\n", conc_block.lstfile);
	printf ("Writing: %s\n", conc_block.outfile);
	printf ("volume:"); for (l = 0; l < conc_block.imgdim[3]; l++) {printf (" %d", l + 1); fflush (stdout);
		conc_read_vol  (&conc_block, imgt);
		conc_write_vol (&conc_block, imgt);
	} printf ("\n"); fflush (stdout);
	conc_ifh_hdr_rec (&conc_block, argc, argv, rcsid);

	switch (conc_block.control) {
		case 'b': case 'B': conc_block.osbig = 1; break;
		case 'l': case 'L': conc_block.osbig = 0; break;
		default: conc_block.osbig = CPU_is_bigendian(); break;
	}
	startrec (conc_block.outfile, argc, argv, rcsid);
	printrec ((conc_block.osbig) ? "bigendian\n" : "littleendian\n");
	printrec ("This is a test.\n");
	catrec (concfile);
	endrec ();
}

	conc_free (&conc_block);
	free (imgt);
	exit (status);
}
