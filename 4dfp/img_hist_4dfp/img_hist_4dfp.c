/*$Header: /data/petsun4/data1/src_solaris/img_hist_4dfp/RCS/img_hist_4dfp.c,v 1.23 2011/01/10 05:09:33 avi Exp $*/
/*$Log: img_hist_4dfp.c,v $
 * Revision 1.23  2011/01/10  05:09:33  avi
 * oprion -C
 *
 * Revision 1.22  2010/12/27  08:03:49  avi
 * toption -f causes selected volume to be reported in filename of -{hpx} created files
 *
 * Revision 1.21  2010/02/11  01:16:37  avi
 * make mask count in auto-range finder
 * (default) NBIN 101->100
 *
 * Revision 1.20  2006/09/23  21:54:37  avi
 * fread() -> eread()
 *
 * Revision 1.19  2006/09/23  21:24:12  avi
 * Solaris 10
 *
 * Revision 1.18  2006/08/11  08:19:59  avi
 * convert moment calculation from histogram to voxel based
 *
 * Revision 1.17  2006/08/11  05:38:38  avi
 * correct output of normalized histograms
 * moments
 * command line saving in output files
 *
 * Revision 1.16  2006/08/09  08:02:18  avi
 * options -s and -u
 *
 * Revision 1.15  2006/06/28  01:04:10  avi
 * new endianio calls (get_4dfp_dimo -> get_4dfp_dimoe)
 * correct bug that cause -r option to fail
 *
 * Revision 1.14  2005/11/11  01:32:42  avi
 * correct bug automatically finding multi-volume voxel value range
 *
 * Revision 1.13  2005/10/11  05:36:40  avi
 * option -x (percentile output)
 *
 * Revision 1.12  2005/03/16  22:00:52  avi
 * tolerate NaN values in input
 *
 * Revision 1.11  2005/03/15  02:33:00  avi
 * mask (-m) option
 *
 * Revision 1.10  2004/09/17  21:15:18  rsachs
 * Replaced 'Get4dfpDimN' with 'get_4dfp_dimo".
 *
 * Revision 1.9  2004/09/17  20:40:53  rsachs
 * Minute cosmetics.
 *
 * Revision 1.8  2003/05/19  23:37:16  avi
 * modernize subroutines and allow histogram builiding over negative voxel values
 *
 * Revision 1.7  1999/08/28  23:15:50  avi
 * correct getrange ()
 *
 * Revision 1.6  1999/07/26  19:36:08  avi
 * -p option
 *
 * Revision 1.5  1999/07/13  02:48:22  avi
 * option -t and -r (threshold and range)
 *
 * Revision 1.4  1999/07/10  02:21:58  avi
 * xmgr_flag
 *
 * Revision 1.3  1998/10/08  22:24:35  mcavoy
 * -b and -f options
 *
 * Revision 1.2  1998/10/07  19:31:15  mcavoy
 * worked on it
 *
 * Revision 1.1  1998/05/08  22:32:59  tscull
 * Revision 1.2  1997/11/26  20:01:56  tscull
 * Revision 1.1  1997/11/13  20:15:58  tscull
 * Initial revision
*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h> 
#include <math.h>
#include <float.h>
#include <assert.h>
#include <endianio.h>

#define MAXL 256
#define NBIN 100	

extern void	find_range (float *vector, float *mask, int size, float *fmin, float *fmax);
extern void	getrange (char *string, float *fmin, float *fmax);

/********************/
/* global variables */
/********************/
static char	rcsid[] = "$Id: img_hist_4dfp.c,v 1.23 2011/01/10 05:09:33 avi Exp $";
static char	program[MAXL];

void write_command_line (FILE *outfp, int argc, char **argv) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
}

int main (int argc, char **argv) {
/*************/
/* image I/O */
/*************/
	FILE		*fp_img, *fp_msk, *fp_out;
	char		imgroot[MAXL], imgfile[MAXL], mskroot[MAXL] = "", mskfile[MAXL], outroot[MAXL], outfile[MAXL];

/**************/
/* processing */
/**************/
	int		*hist, *dist;
	int		imgdim[4], mskdim[4], nbin = NBIN, dimension, nvox, NaNcount, ntot;
	int		orient, orientm, isbig, isbigm;
	float		voxdim[3], voxdimm[3];
	float		*imag, *mask;
	float		fmin, fmax, maxval, minval = 0.0, thresh;
	double		mu[5];		/* moments */

/*********/
/* flags */
/*********/
	int		iframe = 0;
	int		xmgr_flag = 0;
	int		proc_flag = 0;
	int		xtile_flag = 0;
	int		unit_flag = 0;
	int		range_flag = 0;
	int		thresh_flag = 0;
	int		cwd_flag = 0;
	int		moments = 0;
	int		status = 0;

/***********/
/* utility */
/***********/
	char		*ptr, *str, command[MAXL];
	int		c, i, j, k, ii;
	double		q;

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
				case 'b': nbin = atoi (ptr);					*ptr = '\0'; break;
				case 'f': iframe = atoi (ptr);					*ptr = '\0'; break;
				case 't': thresh = atof (ptr);  thresh_flag++;			*ptr = '\0'; break;
				case 'r': getrange (ptr, &minval, &maxval); range_flag++;	*ptr = '\0'; break;
				case 'm': getroot (ptr, mskroot);				*ptr = '\0'; break;
				case 'h': xmgr_flag++;		break;
				case 'p': proc_flag++;		break;
				case 'x': xtile_flag++;		break;
				case 'u': unit_flag++;		break;
				case 's': moments++;		break;
				case 'C': cwd_flag++;		break;
			}
		} else switch (k) {
			case 0:	getroot (argv[i], imgroot); k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:\t%s <(4dfp) image>\n", program);
		printf ("\toption\n");
		printf ("\t-b<int>\tspecify number of bins (default = %d)\n", NBIN);
		printf ("\t-f<int>\tselect volume (counting from 1) of 4dfp stack (default analyze all volumes)\n");
		printf ("\t-t<flt>\tspecify image intensity threshold\n");
		printf ("\t-r<flt>[to<flt>]\tspecify histogram range\n");
		printf ("\t-m<(4dfp) mask>\tmask input using (non-zero voxels of) specified mask\n");
		printf ("\t-h\tcreate <image>.hist file suitable for plotting, e.g., with xmgr\n");
		printf ("\t-p\tcreate <image>.dat  file suitable for input to numerical procedures\n");
		printf ("\t-x\tcreate <image>.xtile percentile listing\n");
		printf ("\t-C\tcreate output files in $cwd (default parallel to <(4dfp) image>)\n");
		printf ("\t-u\tnormalize output .hist and .dat distributions to unit area\n");
		printf ("\t-s\treport moments\n");
		printf ("N.B.:\tonly first frame of mask is used\n");
		printf ("N.B.:\toption -f causes selected volume to be reported in filename of -{hpx} created files\n");
		exit (1);
	} 

/*****************************/
/* get input 4dfp dimensions */
/*****************************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
	nvox = dimension = imgdim[0] * imgdim[1] * imgdim[2];

/*****************/
/* alloc buffers */
/*****************/
	if (!(imag = (float *) malloc (dimension * sizeof (float)))) errm (program);
	if (!(mask = (float *) malloc (dimension * sizeof (float)))) errm (program);
	if (!(hist = (int *)   calloc (nbin, sizeof (int)))) errm (program);

	if (strlen (mskroot)) {
		sprintf (mskfile, "%s.4dfp.img", mskroot);
		if (get_4dfp_dimoe (mskroot, mskdim, voxdimm, &orientm, &isbigm) < 0) errr (program, mskroot);
		status = orientm != orient;
		for (k = 0; k < 3; k++) {
			status |= mskdim[k] - imgdim[k];
			status |= (fabs ((double) (voxdimm[k] - voxdim[k])) > 0.0001);
		}
		if (status) {
			fprintf (stderr, "%s: %s %s dimension mismatch\n", program, imgroot, mskroot);
			exit (-1);
		}
		printf ("Reading: %s\n", mskfile);
		if (!(fp_msk = fopen (mskfile, "rb"))
		|| eread (mask, dimension, isbigm, fp_msk)
		|| fclose (fp_msk)) errr (program, mskfile);
		for (nvox = i = 0; i < dimension; i++) if (mask[i]) nvox++;
	} else {
		for (i = 0; i < dimension; i++) mask[i] = 1.;
	}

/***********/
/* process */
/***********/
	printf ("voxel count = %d\n", nvox);
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (iframe) {
		printf ("Reading: %s frame %d\n", imgfile, iframe);
		if (fseek (fp_img, (long) (iframe - 1) * dimension * sizeof (float), SEEK_SET)
		|| eread (imag, dimension, isbig, fp_img)) errr (program, imgfile);
		if (!range_flag) find_range (imag, mask, dimension, &minval, &maxval);
		for (NaNcount = i = 0; i < dimension; i++) {
			if (mask[i] == 0.0) continue;
			if (isnan (imag[i])) {NaNcount++; continue;}
			if (thresh_flag && imag[i] < thresh) continue;
			j = (int) ((imag[i]-minval)*nbin/(maxval-minval));
			if (j >= nbin) j = nbin - 1;
			j = (j < 0) ? 0 : j;
			hist[j]++;
		}
	} else {
		printf ("Reading: %s\n", imgfile);
		if (!range_flag) {
			minval = +FLT_MAX;
			maxval = -FLT_MAX;
			for (k = 0; k < imgdim[3]; k++) {
				if (eread (imag, dimension, isbig, fp_img)) errr (program, imgfile);
				find_range (imag, mask, dimension, &fmin, &fmax);
				maxval = (fmax > maxval) ? fmax : maxval;
				minval = (fmin < minval) ? fmin : minval;
			}
			rewind (fp_img);
		}
		if (!thresh_flag) thresh = minval - FLT_EPSILON;
		for (NaNcount = k = 0; k < imgdim[3]; k++) {
			if (eread (imag, dimension, isbig, fp_img)) errr (program, imgfile);
			for (i = 0; i < dimension; i++) {
				if (mask[i] == 0.0) continue;
				if (isnan (imag[i])) {NaNcount++; continue;}
				if (thresh_flag && imag[i] < thresh) continue;
				j = (int) ((imag[i]-minval)*nbin/(maxval-minval));
				if (j >= nbin) j = nbin - 1;
				j = (j < 0) ? 0 : j;
				hist[j]++;
			}
		}
	}

	printf ("histogram range (%.4f, %.4f) NaNcount=%d\n", minval, maxval, NaNcount);
	for (i = 0; i < nbin; i++) {
		printf ("%14.6f%10d\n", minval+(i+0.5)*(maxval-minval)/nbin, hist[i]);
	}

/***************************/
/* unit area normalization */
/***************************/
	for (ntot = i = 0; i < nbin; i++) ntot += hist[i];
	q = (unit_flag) ? (maxval-minval)*ntot/(float)nbin : 1.0;

/***********************/
/* compute output root */
/***********************/
	ptr = imgroot;
	if (cwd_flag && (str = strrchr (imgroot, '/'))) ptr = ++str;
	strcpy (outroot, ptr);

	if (xmgr_flag) {
		if (iframe) {
			sprintf (outfile, "%s_vol%d.hist", outroot, iframe);
		} else {
			sprintf (outfile, "%s.hist", outroot);
		}
		if (!(fp_out = fopen (outfile, "w"))) errw (program, outfile);
		printf ("Writing: %s\n", outfile);
		write_command_line (fp_out, argc, argv);
		fprintf (fp_out, "%10.2f%10d\n", minval, 0);
		for (i = 0; i < nbin; i++) {
			if (unit_flag) {
				fprintf (fp_out, "%14.6f %14.6f\n", minval+(i+0)*(maxval-minval)/nbin, hist[i]/q);
				fprintf (fp_out, "%14.6f %14.6f\n", minval+(i+1)*(maxval-minval)/nbin, hist[i]/q);
			} else {
				fprintf (fp_out, "%14.6f %10d\n",   minval+(i+0)*(maxval-minval)/nbin, hist[i]);
				fprintf (fp_out, "%14.6f %10d\n",   minval+(i+1)*(maxval-minval)/nbin, hist[i]);
			}
		}
		fprintf (fp_out, "%10.2f%10d\n", maxval, 0);
		fclose (fp_out);
	}

	if (proc_flag) {
		if (iframe) {
			sprintf (outfile, "%s_vol%d.dat", outroot, iframe);
		} else {
			sprintf (outfile, "%s.dat", outroot);
		}
		if (!(fp_out = fopen (outfile, "w"))) errw (program, outfile);
		printf ("Writing: %s\n", outfile);
		write_command_line (fp_out, argc, argv);
		for (i = 0; i < nbin; i++) {
			if (unit_flag) {
				fprintf (fp_out, "%14.6f %14.6f\n", minval+(i+0.5)*(maxval-minval)/nbin, hist[i]/q);
			} else {
				fprintf (fp_out, "%10.2f %10d\n",   minval+(i+0.5)*(maxval-minval)/nbin, hist[i]);
			}
		}
		fclose (fp_out);
	}

	if (xtile_flag) {
		if (!(dist = (int *) calloc (nbin, sizeof (int)))) errm (program);
		if (iframe) {
			sprintf (outfile, "%s_vol%d.xtile", outroot, iframe);
		} else {
			sprintf (outfile, "%s.xtile", outroot);
		}
		if (!(fp_out = fopen (outfile, "w"))) errw (program, outfile);
		printf ("Writing: %s\n", outfile);
		write_command_line (fp_out, argc, argv);
		fprintf (fp_out, "#    xtile     value\n");
		dist[0] = hist[0]; for (i = 1; i < nbin; i++) dist[i] = dist[i - 1] + hist[i];
		if (0) for (i = 0; i < nbin; i++) {
			printf ("%10.2f%10d%10d%10.2f\n",
				minval+(i+0.0)*(maxval-minval)/nbin, hist[i], dist[i], 100.0*dist[i]/ntot);
		}
		assert (dist[nbin-1] == ntot);
		for (j = 0; j < 100; j++) {
			for (i = 0; i < nbin; i++) if ((float) dist[i]/ntot >= (float) j/100.0) break;
			if (!i) continue;
			for (ii = i - 1; ii > 0; ii--) if (dist[ii] != dist[i]) break;
			q = ((float) ntot*(j/100.) - dist[ii])/((float) dist[i] - dist[ii]);
			if (0) printf ("q=%f\ti=%d\tii=%d\tdist[i]=%d\tdist[ii]=%d\n", q, i, ii, dist[i], dist[ii]);
			fprintf (fp_out, "%10d%10.2f\n", j, minval+(i-1+q)*(maxval-minval)/nbin);
		}
		fclose (fp_out);
		free (dist);
	}

	if (moments) {
		mu[0] = mu[1] = mu[2] = mu[3] = mu[4] = 0.0;
		if (iframe) {
			if (fseek (fp_img, (long) (iframe - 1) * dimension * sizeof (float), SEEK_SET)
			|| eread (imag, dimension, isbig, fp_img)) errr (program, imgfile);
			for (i = 0; i < dimension; i++) {
				if (mask[i] == 0.0) continue;
				if (isnan (imag[i])) continue;
				if (thresh_flag && imag[i] < thresh) continue;
				mu[1] += imag[i];
			}
		} else {
			rewind (fp_img);
			for (k = 0; k < imgdim[3]; k++) {
				if (eread (imag, dimension, isbig, fp_img)) errr (program, imgfile);
				for (i = 0; i < dimension; i++) {
					if (mask[i] == 0.0) continue;
					if (isnan (imag[i])) continue;
					mu[1] += imag[i];
				}
			}
		}
		mu[1] /= ntot;
		if (iframe) {
			if (fseek (fp_img, (long) (iframe - 1) * dimension * sizeof (float), SEEK_SET)
			|| eread (imag, dimension, isbig, fp_img)) errr (program, imgfile);
			for (i = 0; i < dimension; i++) {
				if (mask[i] == 0.0) continue;
				if (isnan (imag[i])) continue;
				if (thresh_flag && imag[i] < thresh) continue;
				q = imag[i] - mu[1];
				mu[2] += q*q;
				mu[3] += q*q*q;
				mu[4] += q*q*q*q;
			}
		} else {
			rewind (fp_img);
			for (k = 0; k < imgdim[3]; k++) {
				if (eread (imag, dimension, isbig, fp_img)) errr (program, imgfile);
				for (i = 0; i < dimension; i++) {
					if (mask[i] == 0.0) continue;
					if (isnan (imag[i])) continue;
					q = imag[i] - mu[1];
					mu[2] += q*q;
					mu[3] += q*q*q;
					mu[4] += q*q*q*q;
				}
			}
		}
		mu[2] /= ntot;
		mu[3] /= ntot;
		mu[4] /= ntot;
		for (k = 1; k <= 4; k++) printf ("u[%d]=%.6f\t", k, mu[k]); printf ("\n");
	}

	fclose (fp_img);
	free (mask);
	free (imag);
	free (hist);
	exit (0);
}

void find_range (float *vector, float *mask, int size, float *minval, float *maxval) {
	int	i;

	*minval =  FLT_MAX;
	*maxval = -FLT_MAX;
	for (i = 0; i < size; i++) {
		if (isnan (vector[i]) || mask[i] == 0.) continue;
		if (vector[i] < *minval) *minval = vector[i];
		if (vector[i] > *maxval) *maxval = vector[i];
	}
}

void getrange (char *string, float *minval, float *maxval) {
	char	*str;

	str = strstr (string, "to");
	if (str) {
		*str = '\0';
		*minval = atof (string);
		*maxval = atof (str + 2);
	} else {
		*minval = 0.0;
		*maxval = atof (string);
	}
}

