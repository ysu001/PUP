/*$Header: /data/petsun4/data1/src_solaris/var_4dfp/RCS/dvar_4dfp.c,v 1.3 2008/02/01 06:25:47 avi Exp $*/
/*$Log: dvar_4dfp.c,v $
 * Revision 1.3  2008/02/01  06:25:47  avi
 * output frame-by-frame listing to a file
 *
 * Revision 1.2  2008/02/01  04:17:10  avi
 * linux and Solaris 10 compliant
 * options -m and s
 * correct reading of mask image
 *
 * Revision 1.1  2001/03/21  02:17:54  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef _FEATURES_H		/* defined by math.h in Linux */
	#include <ieeefp.h>	/* isnan() */
#endif
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>

#define MAXL	256

static char rcsid[] = "$Id: dvar_4dfp.c,v 1.3 2008/02/01 06:25:47 avi Exp $";
int main (int argc, char *argv[]) {
	FILE	*imgfp, *outfp, *datfp;
	IFH	ifh;

/*********/
/* flags */
/*********/
	int	mask = 0;
	int	status = 0;
	int	sd0_flag = 0;
	int	debug = 0;

/***********/
/* utility */
/***********/
	int	c, i, j, k;
	char	*ptr, command[MAXL], program[MAXL];

/*******************/
/* imge processing */
/*******************/
	char	imgroot[MAXL], imgfile[MAXL], datfile[MAXL];
	char	mskroot[MAXL], mskfile[MAXL];
	char	outroot[MAXL], outfile[MAXL];
	float	*imgt, *imgu, *imgm, *imgv;
	int	*imgn;
	float	voxdim[2][3];
	float	threshold = 0.;
	int	imgdim[2][4], dimension, nf_anat = 0, nvox;
	int	isbig[2], orient[2];
	char	control = '\0';
	double	q, qsum;

	fprintf (stdout, "%s\n", rcsid);
	if (ptr = strrchr (argv[0], '/')) ptr++; else ptr = argv[0];
	strcpy (program, ptr);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 's': sd0_flag++;				break;
				case 'm': getroot (ptr, mskroot); mask++;	*ptr = '\0'; break;
				case 'n': nf_anat = atoi (ptr);			*ptr = '\0'; break;
				case 't': threshold = atof (ptr);		*ptr = '\0'; break;
				case '@': control = *ptr++;			*ptr = '\0'; break;
				default:					break;
			}
		} else switch (k) {
			case 0:	getroot (argv[i], imgroot); k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:\t%s [options] <stack_4dfp>\n", program);
		printf (" e.g.,\t%s -n3 test_b1_rmsp_dbnd -mtest_anat_ave_mskt\n", program);
		printf ("\toption\n");
		printf ("\t-m<(4dfp) mask>\tuse specified 4dfp mask\n");
		printf ("\t-n<int>\t\tspecify number of pre-functional (anatomy) frames\n");
        	printf ("\t-t<flt>\t\tspecify maskfile threshold (default = 0.0)\n");
		printf ("\t-s\toutput sqrt(dvar) (default dvar)\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		exit (1);
	}
	printf ("nf_anat=%d\n", nf_anat);

/***********************************/
/* get input 4dfp image dimensions */
/***********************************/
	if (get_4dfp_dimoe (imgroot, imgdim[0], voxdim[0], orient+0, isbig+0)) errr (program, imgroot);
	dimension = imgdim[0][0] * imgdim[0][1] * imgdim[0][2];
	printf ("image dimensions %10d%10d%10d%10d\n",   imgdim[0][0], imgdim[0][1], imgdim[0][2], imgdim[0][3]);
	printf ("image voxdim     %10.6f%10.6f%10.6f\n", voxdim[0][0], voxdim[0][1], voxdim[0][2]);
	if (!control) control = (isbig[0]) ? 'b' : 'l';
	if (Getifh (imgroot, &ifh)) errr (program, imgroot);

	if (!(imgt = (float *) calloc (dimension, sizeof (float)))
	||  !(imgu = (float *) calloc (dimension, sizeof (float)))
	||  !(imgv = (float *) calloc (dimension, sizeof (float)))
	||  !(imgm = (float *) calloc (dimension, sizeof (float)))
	||  !(imgn = (int *)   calloc (dimension, sizeof (int)))) errm (program);

	if (mask) {
		if (get_4dfp_dimoe (mskroot, imgdim[1], voxdim[1], orient+1, isbig+1)) errr (program, mskfile);
		status = orient[1] - orient[0];
		for (k = 0; k < 3; k++) {
			status |= imgdim[1][k] - imgdim[0][k];
			status |= (fabs (voxdim[1][k] - voxdim[0][k]) > 1.0e-4);
		}
		if (status) {
			fprintf (stderr, "%s: %s %s dimension/orientation mismatch\n", program, imgroot, mskroot);
			exit (-1);
		}
		sprintf (mskfile, "%s.4dfp.img", mskroot);
		fprintf (stdout, "Reading: %s\n", mskfile);
		if (!(imgfp = fopen (mskfile, "rb"))
		|| eread (imgm, dimension, isbig[1], imgfp)
		|| fclose (imgfp)) errr (program, mskfile);
		for (nvox = i = 0; i < dimension; i++) if (imgm[i] > threshold) nvox++;
		printf ("Within-mask voxel count = %d\n", nvox);
	}

/***********/
/* compute */
/***********/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	fprintf (stdout, "Reading: %s\n", imgfile);
	if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
	sprintf (datfile, "%s.ddat", imgroot);
	fprintf (stdout, "Writing: %s\n", datfile);
	if (!(datfp = fopen (datfile, "w"))) errw (program, datfile);
	for (j = 0; j < imgdim[0][3]; j ++) {
		if (eread (imgt, dimension, isbig[0], imgfp)) errr (program, imgfile);
		for (qsum = nvox = i = 0; i < dimension; i++) {
			if (mask && imgm[i] <= threshold) continue;
			if (imgt[i] == (float) 1.e-37 || isnan (imgt[i]) || imgt[i] == 0.0) continue;
			if (imgu[i] == (float) 1.e-37 || isnan (imgu[i]) || imgu[i] == 0.0) {if (j) continue;}
			q = imgt[i] - imgu[i];
			nvox++;
			imgu[i] = imgt[i];
			qsum += q*q;
			if (j > nf_anat) {
				imgv[i] += q*q;
				imgn[i]++;
			}
		}
		if (j) {
			fprintf (datfp, "frame %3d sqrt(<(dI/dt)^2>): %10.4f\n", j + 1, sqrt (qsum/nvox));
		} else {
			fprintf (datfp, "frame %3d sqrt(<(dI/dt)^2>): %10.4f\n", j + 1, 500.);
		}
	}
	fclose (imgfp);
	fclose (datfp);

/****************************/
/* normalize summed vaiance */
/****************************/
	for (i = 0; i < dimension; i++) {
		if (imgn[i]) {
			imgv[i] /= imgn[i];
			if (sd0_flag) imgv[i] = sqrt (imgv[i]);
		} else {
			imgv[i] = 1.e-37;
		}
	}
				
/***********************/
/* write output volume */
/***********************/
	sprintf (outroot, "%s_%s", imgroot, ((sd0_flag) ? "dsd0" : "dvar"));
	sprintf (outfile, "%s.4dfp.img", outroot);
	fprintf (stdout, "Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "wb"))
	|| ewrite (imgv, dimension, control, outfp)
	|| fclose (outfp)) errw (program, outfile);

/***************/
/* ifh and hdr */
/***************/
	imgdim[0][3] = 1;
	writeifhmce (program, outfile, imgdim[0],
				ifh.scaling_factor, ifh.orientation, ifh.mmppix, ifh.center, control);
	sprintf (command, "ifh2hdr %s", outroot); status |= system (command);

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, rcsid, control);
	sprintf (command, "Differentiated time series %s\n", ((sd0_flag) ? "standard deviation" : "variance"));
	printrec (command);
	sprintf (command, "Number of skipped frames = %d\n", nf_anat);
	printrec (command);
	if (mask) {sprintf (command, "Within-mask voxel count = %d\n", nvox); printrec (command);}
	catrec (imgfile);
	if (mask) catrec (mskfile);
	endrec ();

	free (imgt); (imgm); free (imgu); free (imgv); free (imgn);
	exit (0);
}
