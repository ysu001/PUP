/*$Header: /data/petsun4/data1/src_solaris/mprtir/RCS/2Dhist.c,v 1.13 2007/08/08 23:21:55 avi Exp $*/
/*$Log: 2Dhist.c,v $
 * Revision 1.13  2007/08/08  23:21:55  avi
 * protect against zero binwidth (binwid1 and binwid2)
 *
 * Revision 1.12  2007/08/07  20:46:59  avi
 * endian compliant
 * separate imgo from img1
 *
 * Revision 1.11  2007/08/07  00:44:51  avi
 * -r (voxel value range) option
 *
 * Revision 1.10  2003/05/05  23:37:22  avi
 * -g (suppress grid burn-in) option
 *
 * Revision 1.9  1999/03/01  02:23:28  avi
 * CUM_D -> 0.000005
 * -d
 *
 * Revision 1.8  1998/07/11  02:36:50  avi
 * grid
 * Revision 1.7  1998/07/09  04:27:03  avi
 * limit checks on k in mrg1[k]++; mrg1[k]++;
 * rec file
 * Revision 1.6  1998/07/01  05:39:01  avi
 * in range modification loop discount first bin for mrg2 as well as mrg1
 * Revision 1.5  1998/03/17  08:51:35  avi
 * correct range tayloring
 * allow negative voxel values
 * Revision 1.4  1998/03/16  08:26:13  avi
 * rcs $Id: 2Dhist.c,v 1.13 2007/08/08 23:21:55 avi Exp $ and $Log: 2Dhist.c,v $
 * Revision 1.13  2007/08/08  23:21:55  avi
 * protect against zero binwidth (binwid1 and binwid2)
 *
 * Revision 1.12  2007/08/07  20:46:59  avi
 * endian compliant
 * separate imgo from img1
 *
 * Revision 1.11  2007/08/07  00:44:51  avi
 * -r (voxel value range) option
 *
 * Revision 1.10  2003/05/05  23:37:22  avi
 * -g (suppress grid burn-in) option
 *
 * Revision 1.9  1999/03/01  02:23:28  avi
 * CUM_D -> 0.000005
 * -d
 *
 * Revision 1.8  1998/07/11  02:36:50  avi
 * grid
 * Revision 1.7  1998/07/09  04:27:03  avi
 * limit checks on k in mrg1[k]++; mrg1[k]++;
 * rec file
 * Revision 1.6  1998/07/01  05:39:01  avi
 * in range modification loop discount first bin for mrg2 as well as mrg1
 * Revision 1.5  1998/03/17  08:51:35  avi
 * correct range tayloring
 * allow negative voxel values
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <endianio.h>
#include <rec.h>

#define MAXL	256
#define NBIN	256
#define CUM_D	0.000005

void getrangei (char *string, int *minval, int *maxval) {
	char	*str;

	str = strstr (string, "to");
	if (str) {
		*str = '\0';
		*minval = atoi (string);
		*maxval = atoi (str + 2);
	} else {
		*minval = 0;
		*maxval = atoi (string);
	}
}

void erri (char* program, char* filespc) {
	fprintf (stderr, "%s: %s not short int\n", program, filespc);
	exit (-1);
}

void usage (char *program) {
	printf ("usage:\t%s <img1> <img2> <2Dhist_out>\n", program);
	printf ("\tooptions\n");
	printf ("\t-g\tsuppress grid burn-in\n");
	printf ("\t-r<1|2>:[<int>to]<int>\tspecify usable voxel value range for <img1> or <img2>\n");
	printf ("\t-@<b|l>\toutput big or little endian (default CPU endian)\n");
	printf ("e.g.,\t%s va1701_tir_ecatt_ANALYZE va1701_mpr_S va1701_tir_mpr_2Dhist\n", program);
	printf ("e.g.,\t%s va1701_tir_ecatt_ANALYZE va1701_mpr_S va1701_tir_mpr_2Dhist -r1:0to2500\n", program);
	printf ("e.g.,\t%s va1701_tir_ecatt_ANALYZE va1701_mpr_S va1701_tir_mpr_2Dhist -r1:2500\n", program);
	printf ("N.B.:\t%s operates on short int (ANALYZE) format images\n", program);
	exit (1);
}

static char rcsid[] = "$Id: 2Dhist.c,v 1.13 2007/08/08 23:21:55 avi Exp $";
int main (int argc, char *argv[]) {
	FILE 		*anafp;
	char            in1root[MAXL];
	char            in2root[MAXL];
	char		outroot[MAXL];
	char		imgfile[MAXL];
	char            program[MAXL], string[MAXL], *ptr;

	struct dsr	hdr;
	short		*img1, *img2, *imgo;
	int		xdim, ydim, zdim, dimension;
	int		img1max = 0, img2max = 0;
	int		img1min = 0, img2min = 0;
	int		swab_img1 = 0, swab_img2 = 0, swab_imgo = 0;
	char		orient, control = '\0';

	int		nbin = NBIN;
	int		*hist, *mrg1, *mrg2;
	int		histmax = 0;
	int		range1, range2;
	int		binwid1, binwid2;
	int		grid1, grid2;
	short		valmax = 0;

	float		ftmp;
	int		c, i, j, k;

/*********/
/* flags */
/*********/
	int		status = 0;
	int		debug = 0;
	int		do_grid = 1;

	if (ptr = strrchr (argv[0], '/')) ptr++; else ptr = argv[0];
	strcpy (program, ptr);
	printf ("%s\n", rcsid);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) if (*argv[i] == '-') {
		strcpy (string, argv[i]); ptr = string;
		while (c = *ptr++) switch (c) {
			case 'd': debug++;		break;
			case 'g': do_grid = 0;		break;
			case '@': control = *ptr++;	*ptr = '\0'; break;
			case 'r': j = atoi (ptr++);
				ptr++;		/* skip over ":" */
				switch (j) {
					case 1: getrangei (ptr, &img1min, &img1max);	break;
					case 2: getrangei (ptr, &img2min, &img2max);	break;
					default: usage (program);			break;
				}			*ptr = '\0'; break;
		}
	} else switch (k) {
		case 0: getroot (argv[i], in1root);	k++; break;
		case 1: getroot (argv[i], in2root);	k++; break;
		case 2: getroot (argv[i], outroot);	k++; break;
	}
	if (k < 3) usage (program);
	fprintf (stdout, "img1: %s\nimg2: %s\n", in1root, in2root);

	sprintf (imgfile, "%s.hdr", in1root);
	if (!(anafp = fopen (imgfile, "rb")) || fread (&hdr, sizeof (struct dsr), 1, anafp) != 1
	|| fclose (anafp)) errr (program, imgfile);
	if (hdr.hk.sizeof_hdr != sizeof (struct dsr)) {
		printf ("converting %s byte order\n", in1root);
		swab_hdr (&hdr);
		swab_img1++;
	}
	if (hdr.dime.bitpix != 16) erri (program, imgfile);
	xdim = hdr.dime.dim[1];
	ydim = hdr.dime.dim[2];
	zdim = hdr.dime.dim[3];
	orient = hdr.hist.orient;
	dimension = xdim * ydim * zdim;

	sprintf (imgfile, "%s.hdr", in2root);
	if (!(anafp = fopen (imgfile, "rb")) || fread (&hdr, sizeof (struct dsr), 1, anafp) != 1
	|| fclose (anafp)) errr (program, imgfile);
	if (hdr.hk.sizeof_hdr != sizeof (struct dsr)) {
		printf ("converting %s byte order\n", in2root);
		swab_hdr (&hdr);
		swab_img2++;
	}
	if (hdr.dime.bitpix != 16) erri (program, imgfile);
	if (xdim != hdr.dime.dim[1] || ydim != hdr.dime.dim[2] || zdim != hdr.dime.dim[3] || orient != hdr.hist.orient) {
		fprintf (stderr, "%s: %s %s dimension/orientation mismatch\n", program, in1root, in2root);
		exit (-1);
	}

	img1 = (short *) malloc (dimension * sizeof (short));
	img2 = (short *) malloc (dimension * sizeof (short));
	mrg1 = (int *) calloc (nbin, sizeof (int));
	mrg2 = (int *) calloc (nbin, sizeof (int));
	hist = (int *) calloc (nbin*nbin, sizeof (int));
	imgo = (short *) malloc (nbin*nbin * sizeof (short));
	if (!img1 || !img2 || !hist || !mrg1 || !mrg2 || !imgo) errm (program);

	sprintf (imgfile, "%s.img", outroot); startrece (imgfile, argc, argv, rcsid, control);
	sprintf (imgfile, "%s.img", in1root); catrec (imgfile);
	printf ("Reading: %s\n", imgfile);
	if (!(anafp = fopen (imgfile, "rb"))
	|| fread (img1, sizeof (short), dimension, anafp) != dimension || fclose (anafp)) errr (program, imgfile);
	if (swab_img1) for (i = 0; i < dimension; i++) swab2 ((char *) &img1[i]);

	sprintf (imgfile, "%s.img", in2root); catrec (imgfile);
	printf ("Reading: %s\n", imgfile);
	if (!(anafp = fopen (imgfile, "rb"))
	|| fread (img2, sizeof (short), dimension, anafp) != dimension || fclose (anafp)) errr (program, imgfile);
	if (swab_img2) for (i = 0; i < dimension; i++) swab2 ((char *) &img2[i]);

	if (!img1max && !img1min) {
		img1max -32768; img1min = 32767;
		for (i = 0; i < dimension; i++) {
			if (img1[i] > img1max) img1max = img1[i];
			if (img1[i] < img1min) img1min = img1[i];
		}
	} else {
		for (i = 0; i < dimension; i++) {
			if (img1[i] > img1max) img1[i] = img1max;
			if (img1[i] < img1min) img1[i] = img1min;
		}
	}
	range1 = img1max - img1min;
	if (!img2max && !img2min) {
		img2max -32768; img2min = 32767;
		for (i = 0; i < dimension; i++) {
			if (img2[i] > img2max) img2max = img2[i];
			if (img2[i] < img2min) img2min = img2[i];
		}
	} else {
		for (i = 0; i < dimension; i++) {
			if (img2[i] > img2max) img2[i] = img2max;
			if (img2[i] < img2min) img2[i] = img2min;
		}
	}
	range2 = img2max - img2min;

	fprintf (stdout, "img1: min=%6d\tmax=%6d\trange=%6d\n", img1min, img1max, range1);
	fprintf (stdout, "img2: min=%6d\tmax=%6d\trange=%6d\n", img2min, img2max, range2);
	sprintf (string, "before auto range adjust\n"); printrec (string);
	sprintf (string, "img1: min=%6d\tmax=%6d\trange=%6d\n", img1min, img1max, range1); printrec (string);
	sprintf (string, "img2: min=%6d\tmax=%6d\trange=%6d\n", img2min, img2max, range2); printrec (string);

	for (i = 0; i < dimension; i++) {
		k = (nbin * (img1[i] - img1min)) / range1;
		if (k < nbin && k > 0) mrg1[k]++;
		k = (nbin * (img2[i] - img2min)) / range2;
		if (k < nbin && k > 0) mrg2[k]++;
	}
	for (k = 2; k < nbin; k++) {
		mrg1[k] += mrg1[k - 1];
		mrg2[k] += mrg2[k - 1];
	}
	if (debug) for (k = 0; k < nbin; k++) {
		printf ("%6d  %10.6f  %10.6f\n", k,
			(float) mrg1[k] / (float) mrg1[nbin - 1], (float) mrg2[k] / (float) mrg2[nbin - 1]);
	}
	for (i = 1; i < nbin; i++) if (((float) mrg1[i] / (float) mrg1[nbin - 1]) > CUM_D) break;
	for (j = 1; j < nbin; j++) if (((float) mrg1[j] / (float) mrg1[nbin - 1]) > 1.0-CUM_D) break;
	if (debug) printf ("img1: first_bin=%d\tlast_bin=%d\n", i, j);
	img1min += (float) (range1 * i) / nbin;
	range1  *= (float) (j - i) / nbin;
	binwid1 = (float) range1 / nbin;
	if (!binwid1) binwid1++;
	range1 = binwid1 * nbin;
	img1min += binwid1 - (img1min + range1) % binwid1;
	for (i = 1; i < nbin; i++) if (((float) mrg2[i] / (float) mrg2[nbin - 1]) > CUM_D) break;
	for (j = 1; j < nbin; j++) if (((float) mrg2[j] / (float) mrg2[nbin - 1]) > 1.0-CUM_D) break;
	if (debug) printf ("img2: first_bin=%d\tlast_bin=%d\n", i, j);
	img2min += (float) (range2 * i) / nbin;
	range2  *= (float) (j - i) / nbin;
	binwid2 = (float) range2 / nbin;
	if (!binwid2) binwid2++;
	range2 = binwid2 * nbin;
	img2min += binwid2 - (img2min + range2) % binwid2;
	fprintf (stdout, "img1: min=%6d\tmax=%6d\trange=%6d\tbinwidth=%d\n", img1min, img1min + range1, range1, binwid1);
	fprintf (stdout, "img2: min=%6d\tmax=%6d\trange=%6d\tbinwidth=%d\n", img2min, img2min + range2, range2, binwid2);
	sprintf (string, "after auto range adjust\n");
		printrec (string);
	sprintf (string, "img1: min=%6d\tmax=%6d\trange=%6d\tbinwidth=%d\n", img1min, img1min + range1, range1, binwid1);
		printrec (string);
	sprintf (string, "img2: min=%6d\tmax=%6d\trange=%6d\tbinwidth=%d\n", img2min, img2min + range2, range2, binwid2);
		printrec (string);
	if (!range1 || !range2) exit (-1);

	for (k = 0; k < dimension; k++) {
		i = (img1[k] - img1min) / binwid1;
		j = (img2[k] - img2min) / binwid2;
		if (i < 0 || i >= nbin || j < 0 || j >= nbin) continue;
		hist[nbin * j + i]++;
		if (!hist[nbin * j + i]) {
			fprintf (stderr, "2Dhist: histogram bin overflow\n");
			exit (-1);
		}
		if (hist[nbin * j + i] > histmax) histmax = hist[nbin * j + i];
	}
	fprintf (stdout, "maximum bin count=%d\n", histmax);

/***************************/
/* write histogram to imgo */
/***************************/
	for (k = 0; k < nbin*nbin; k++) {
		if (hist[k]) {
			ftmp = 100 * log ((double) hist[k]);
			imgo[k] = (short) ftmp;
		} else imgo[k] = 0;
		if (imgo[k] > valmax) valmax = imgo[k];
	}

	if (do_grid) {
		grid1 = range1 / 10; grid1 /= 10; grid1 *= 10;
		for (i = k = 0; i < nbin; k += grid1) {
			i = (k - img1min) / binwid1;
			if (i >= 0 && i < nbin) for (j = 0; j < nbin; j++) imgo[nbin * j + i] = valmax / 2;
		}
		sprintf (string, "img1: grid interval=%d\n", grid1); printrec (string);

		grid2 = range2 / 10; grid2 /= 10; grid2 *= 10;
		for (j = k = 0; j < nbin; k += grid2) {
			j = (k - img2min) / binwid2;
			if (j >= 0 && j < nbin) for (i = 0; i < nbin; i++) imgo[nbin * j + i] = valmax / 2;
		}
		sprintf (string, "img2: grid interval=%d\n", grid2); printrec (string);
	}

	swab_imgo = ((CPU_is_bigendian() != 0) && (control == 'l' || control == 'L'))
		 || ((CPU_is_bigendian() == 0) && (control == 'b' || control == 'B'));
	if (swab_imgo) for (i = 0; i < nbin*nbin; i++) swab2 ((char *) &imgo[i]);
	sprintf (imgfile, "%s.img", outroot); fprintf (stdout, "Writing: %s\n", imgfile);
	if (!(anafp = fopen (imgfile, "wb")) || fwrite (imgo, sizeof (short), nbin*nbin, anafp) != nbin*nbin
	|| fclose (anafp)) errw (program, imgfile);

	hdr.dime.dim[0] = 2;
	hdr.dime.dim[1] = nbin;
	hdr.dime.dim[2] = nbin;
	hdr.dime.dim[3] = 1;
	hdr.dime.datatype = 4;
	hdr.dime.pixdim[1] = 1;
	hdr.dime.pixdim[2] = 1;
	hdr.dime.pixdim[3] = 1;
	hdr.dime.glmax = valmax;
	hdr.dime.glmin = 0;
	hdr.hist.orient = 0;
	if (swab_imgo) swab_hdr (&hdr);
	sprintf (imgfile, "%s.hdr", outroot);
	fprintf (stdout, "Writing: %s\n", imgfile);
	if (!(anafp = fopen (imgfile, "wb"))
	|| fwrite (&hdr, sizeof (struct dsr), 1, anafp) != 1 || fclose (anafp)) errw (program, imgfile);

	free (img1); free (img2); free (imgo);
	free (mrg1); free (mrg2);
	free (hist);
	endrec ();

	exit (status);
}
