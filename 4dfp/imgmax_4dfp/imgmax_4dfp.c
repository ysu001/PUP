/*$Header: /data/petsun4/data1/src_solaris/imgmax_4dfp/RCS/imgmax_4dfp.c,v 1.11 2006/09/24 01:29:40 avi Exp $*/
/*$Log: imgmax_4dfp.c,v $
 *Revision 1.11  2006/09/24 01:29:40  avi
 *Solaris 10
 *
 * Revision 1.10  2005/05/03  02:55:24  avi
 * report NaNcount
 *
 * Revision 1.9  2005/02/02  21:43:23  avi
 * -r (root sum of squares) option
 *
 * Revision 1.8  2004/12/07  05:32:44  avi
 * option -e (scientific notation)
 *
 * Revision 1.7  2004/09/22  01:05:12  avi
 * use get_4dfp_dimo_quiet () in non-verbose mode
 *
 * Revision 1.6  2004/09/21  20:13:43  rsachs
 * Installed 'errm','errr','errw','setprog'. Also implemented calculation of imgmin.
 *
 * Revision 1.5  1998/12/19  04:58:50  avi
 * *** empty log message ***
 *
 * Revision 1.4  1998/12/19  04:54:22  avi
 * verbose switch
 *
 * Revision 1.3  1998/10/12  19:15:21  mcavoy
 *
 * Revision 1.1  1997/11/13  17:25:17  tscull
 * Initial revision
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <endianio.h>

#define MAXL 256

static char rcsid[] = "$Id: imgmax_4dfp.c,v 1.11 2006/09/24 01:29:40 avi Exp $";

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

int main (int argc, char **argv) {
/*************/
/* image I/O */
/*************/
        FILE            *imgfp;
        char            *ptr, imgfile[MAXL], imgroot[MAXL];

/**************/
/* processing */
/**************/
        int             imgdim[4], dimension, orient, NaNcount, isbig;
        float           voxdim[3];
        float           *imgt;
	float		frame_max, stack_max, frame_min, stack_min;
	double		ssq, frame_ssq;

/*********/
/* flags */
/*********/
	int		verbose = 0;
	int		min = 0;
	int		eflag = 0;
	int		rssq_flag = 0;

/***********/
/* utility */
/***********/
	char		program[MAXL], command[MAXL];
	int		c, i, k, l;

	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'v': verbose++; 	break;
				case 'm': min++;	break;
				case 'e': eflag++;	break;
				case 'r': rssq_flag++;	break;
			}
		}
		else switch (k) {
			case 0:	getroot (argv[i], imgroot); k++; break;
		}	
	}
        if (verbose) printf ("%s\n", rcsid);
	if (k < 1) {
		printf ("Usage:\t%s <my_image[.4dfp.img]>\n", program);
		printf ("\toption\n");
		printf ("\t-m\treport min as well as max\n");
		printf ("\t-e\treport max/min values in scientific notation\n");
		printf ("\t-r\treport root sum of squares\n");
		printf ("\t-v\tverbose (time series) mode\n");
		exit (0);
      	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
/*****************************/
/* get input 4dfp dimensions */
/*****************************/
	if (verbose) {
	        if (get_4dfp_dimoe	 (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	} else {
	        if (get_4dfp_dimoe_quiet (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	}	
        dimension = imgdim[0] * imgdim[1] * imgdim[2];

/*****************/
/* alloc buffers */
/*****************/
        imgt = (float *) malloc (dimension * sizeof (float));
        if (!imgt) errm (program);

/***********/
/* process */
/***********/
        if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
        if (verbose) fprintf (stdout, "Reading: %s\n", imgfile);
	stack_max = -FLT_MAX; stack_min = FLT_MAX;
	ssq = 0.0;
        for (l = 0; l < imgdim[3]; l++) {
                if (eread (imgt, dimension, isbig, imgfp)) errr (program, imgfile);
		frame_max = -FLT_MAX; frame_min = FLT_MAX;
		frame_ssq = 0;
		for (NaNcount = k = 0; k < dimension; k++) {
			if (isnan (imgt[k])) {NaNcount++; continue;}
			if (imgt[k] > frame_max) frame_max = imgt[k];
			if (imgt[k] < frame_min) frame_min = imgt[k];
			frame_ssq += imgt[k]*imgt[k];
		}
		if (verbose) {
			if (eflag) {
					printf ("frame %3d\tmax = %e\n", l + 1, frame_max);
			} else {
					printf ("frame %3d\tmax = %f\n", l + 1, frame_max);
			}
			if (min) {
				if (eflag) {
					printf ("frame %3d\tmin = %e\n", l + 1, frame_min);
				} else {
					printf ("frame %3d\tmin = %f\n", l + 1, frame_min);
				}
			}
			if (rssq_flag) {
				printf ("frame %3d\troot sum of squares = %f\n", l + 1, sqrt (frame_ssq));
			}
			if (NaNcount) printf ("NaNcount=%d\n", NaNcount);
		}
		if (frame_max > stack_max) stack_max = frame_max;
		if (frame_min < stack_min) stack_min = frame_min;
		ssq += frame_ssq;
        }
        fclose (imgfp);
	if (!verbose) {
		if (rssq_flag) {
			printf ("%f", sqrt (ssq));
		} else {
			if (eflag) {
					printf ("%e",   stack_max);
			} else {
					printf ("%f",   stack_max);
			}
			if (min) {
				if (eflag) {
					printf ("\t%e", stack_min);
				} else {
					printf ("\t%f", stack_min);
				}
			}
		}
		printf ("\n");
	}
	free (imgt);
	exit (0);
}
