/*$Header: /data/petsun4/data1/src_solaris/scale_4dfp/RCS/scale_4dfp.c,v 1.14 2007/07/05 05:01:52 avi Exp $*/
/*$Log: scale_4dfp.c,v $
 * Revision 1.14  2007/07/05  05:01:52  avi
 * startrec() -> startrece()
 *
 * Revision 1.13  2006/09/29  23:36:57  avi
 * Solaris 10
 *
 * Revision 1.12  2004/09/23  20:58:31  rsachs
 * Installed 'errm','errr','errw','setprog' & 'get_4dfp_dimo'.
 *
 * Revision 1.11  2004/05/10  02:01:23  avi
 * -E option (preserve 1.0e-37)
 *
 * Revision 1.10  2002/04/24  21:41:01  avi
 * -b (add constant) option
 *
 * Revision 1.9  2001/11/03  00:21:47  avi
 * protect against interpreting switch as second command line argument
 *
 * Revision 1.8  2001/05/18  01:17:57  avi
 * better output recfile logic
 *
 * Revision 1.7  2000/04/06  01:58:55  avi
 * allow negative scale factors
 * extensive revisions to manage overwriting
 * include math.h (potentially a serious problem hitherto)
 *
 * Revision 1.6  1998/10/12  21:59:12  mcavoy
 * Revision 1.5  1998/10/02  06:29:42  avi
 * minor typo corrections
 *
 * Revision 1.4  1998/10/01  23:42:51  mcavoy
 *
 * Revision 1.3  1998/05/07  01:25:05  avi
 * -a (out file name trailer) option
 * bring c source style into conformity with standard
 *
 * Revision 1.2  1998/05/06  23:44:21  tscull
 * Revision 1.1  1997/01/26  00:37:38  avi
 * Initial revision
 **/
/*_________________________________________________________________________
Module:		scale_4dfp.c
Date:		24-Jan-97
Authors:	Avi Snyder
Description:	Scale and write back 4dfp stack.
_________________________________________________________________________*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <endianio.h>				/* getpid () */
#include <Getifh.h>
#include <rec.h>

#define MAXL	256

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[] = "$Id: scale_4dfp.c,v 1.14 2007/07/05 05:01:52 avi Exp $";
int main (int argc, char *argv[]) {
/*************/
/* image I/O */
/*************/
	FILE 		*fp_img, *fp_out;
	char            imgfile[MAXL], recfile[MAXL];
	char            outfile[MAXL];
	char 		imgroot[MAXL];
	char 		outroot[MAXL];
	char		trailer[MAXL];			/* optional output filename trailer */

/**************/
/* processing */
/**************/
	IFH		ifh;
        int  		imgdim[4], dimension, orient, isbig;
        float	 	voxdim[3];	
	float		*image3d;
	float		factor = 1.;			/* multiplicative scaling factor */
	float		addend = 0.;			/* additive constant */
	char		control = '\0';

/***********/
/* utility */
/***********/
	int		c, i, k;
	char 		*ptr, command[4*MAXL], program[MAXL];

/*********/
/* flags */
/*********/
	int		test = 0;
	int		status = 0;
	int		trl_flag = 0;			/* use output filename trailer */
	int		e37_flag = 0;			/* preserve 1.0e-37 values */
	int		exists_oldrec;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-' && i > 2) {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'E': e37_flag++;				break;
				case '@': control = *ptr++;			*ptr = '\0'; break;
				case 'a': strcpy (trailer, ptr); trl_flag++;	*ptr = '\0'; break;
				case 'b': addend = atof (ptr);			*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0:	getroot (argv[i], imgroot);			k++; break;
			case 1:	test = isdigit (argv[i][0]) || argv[i][0] == '.' || argv[i][0] == '-';
				if (test) {factor = atof (argv[i]); k++;}
				break;
		}	
	}
	if (k < 2) {
		printf ("Usage:\tscale_4dfp <image_4dfp> <scale_factor> [options]\n");
		printf ("\toption\n");
		printf ("\t-E\tpreserve 1.0e-37 values (fidl NaN convention)\n");
		printf ("\t-a<str>\tappend trailer to output file name\n");
		printf ("\t-b<flt>\tadd specified constant to each voxel\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("e.g.,\tscale_4dfp b2_xfrm_avg 12\n");
		printf ("e.g.,\tscale_4dfp b2_xfrm_avg 12 -b5 -ax12+5\n");
		printf ("N.B.:\t<image_4dfp> is overwritten unless the trailer option is used\n");
		printf ("N.B.:\t<scale_factor> must be specified for proper operation\n");
		exit (1);
	}
	printf ("%s: factor=%f\n", program, factor);
	printf ("%s: addend=%f\n", program, addend);

        if (trl_flag) {
                sprintf (outroot, "%s_%s", imgroot, trailer);
        } else {
		strcpy (outroot, imgroot);
	}

/*****************************/
/* get input 4dfp dimensions */
/*****************************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	if (Getifh (imgfile, &ifh)) errr (program, imgroot);
	if (!control) control = (isbig) ? 'b' : 'l';
	dimension = imgdim[0] * imgdim[1] * imgdim[2];

/****************/
/* alloc buffer */
/****************/
	if (!(image3d = (float *) malloc (dimension * sizeof (float)))) errm (program);

/***********/
/* process */
/***********/
	if (trl_flag) {
		if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
		if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	} else {
		if (!(fp_img = fopen (imgfile, "r+b"))) errr (program, imgfile);
		fp_out = fp_img;
	}
	fprintf (stdout, "Reading: %s\n", imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);
	for (k = 0; k < imgdim[3]; k++) {
		fseek (fp_img, (long) k * dimension * sizeof (float), SEEK_SET);
		if (eread (image3d, dimension, isbig, fp_img)) errr (program, imgfile);
		for (i = 0; i < dimension; i++) {
			if (e37_flag && (image3d[i] == (float) 1.0e-37)) continue;
			image3d[i] *= factor;
			if (addend != 0.0) image3d[i] += addend;
		}
		fseek (fp_out, (long) k * dimension * sizeof (float), SEEK_SET);
		if (ewrite (image3d, dimension, control, fp_out)) errw (program, outfile);
	}
	fclose (fp_img);
	if (trl_flag) fclose (fp_out);

/*******************/
/* create rec file */
/*******************/
	sprintf (recfile, "%s.rec", imgfile);
	exists_oldrec = !access (recfile, R_OK);
	if (exists_oldrec) {
		sprintf  (command, "/bin/cp %s.rec %s.rec0", imgfile, imgfile);	status |= system (command);
	}
	startrece (outfile, argc, argv, rcsid, control);
		sprintf  (command, "Scaled by %f\n", factor); printrec (command);
	if (addend != 0.0) {
		sprintf  (command, "Offset by %f\n", addend); printrec (command);
	}
	if (e37_flag) {
		sprintf  (command, "1.0e-37 values preserved (fidl convention)\n"); printrec (command);
	}
	if (exists_oldrec) {
		sprintf  (command, "cat %s.rec0 >> %s.rec",  imgfile, outfile);	status |= system (command);
		sprintf  (command, "/bin/rm %s.rec0", imgfile);			status |= system (command);
		if (status) {
			fprintf (stderr, "%s: %s.rec create error\n", program, outfile);
		}
	} else {
		sprintf (command, "%s not found\n", recfile); printrec (command);
	}
	endrec ();

	if (trl_flag) {
/***************/
/* ifh and hdr */
/***************/
		if (Writeifh (program, outfile, &ifh, control)) errw (program, outfile);
		sprintf (command, "ifh2hdr %s", outroot); printf ("%s\n", command);
		status = system (command);
	}

	free (image3d);
	exit (status);
}
