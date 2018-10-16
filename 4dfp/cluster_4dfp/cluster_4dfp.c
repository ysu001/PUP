/*$Header: /data/petsun4/data1/src_solaris/cluster_4dfp/RCS/cluster_4dfp.c,v 1.14 2012/06/14 06:01:42 avi Exp $*/
/*$Log: cluster_4dfp.c,v $
 * Revision 1.14  2012/06/14  06:01:42  avi
 * option -a
 *
 * Revision 1.13  2011/06/19  03:39:20  avi
 * trap u2 < 0. from rounding error
 *
 * Revision 1.12  2011/01/08  01:44:16  avi
 * include ROI mean and s.d. in output list (.dat) file
 *
 * Revision 1.11  2011/01/08  00:59:07  avi
 * correct output ifh and hdr to always have nvol = 1
 * disambiguate usage regarding outputs using -n and -R
 *
 * Revision 1.10  2010/12/22  06:37:17  avi
 * option -f
 *
 * Revision 1.9  2010/07/12  06:13:09  avi
 * ensure that all voxels in cluster have same sign
 *
 * Revision 1.8  2010/04/14  02:44:56  avi
 * include total voxel count in output rec and dat files
 *
 * Revision 1.7  2009/08/06  01:37:42  avi
 * correct appearance of NaN in s.d. listing attributable to single vs. double arithmetic
 *
 * Revision 1.6  2009/07/02  05:48:15  avi
 * types in XYZN UCHAR -> USHORT
 * analyzeto4dfp_flip on read and write (index2atl -af compliant)
 *
 * Revision 1.5  2007/09/06  04:57:58  avi
 * option -R
 * include command line in -l (dat file) output
 *
 * Revision 1.4  2007/08/23  04:27:32  avi
 * eliminate hard coded MAXR limit (use realloc)
 * vebose test dependency changed to region number (was region count)
 * MAMB 4096 -> 128 (major memory savings)
 * iv0 strategy -> huge speed improvement
 * imgo type short -> unsigned int
 *
 * Revision 1.3  2007/07/20  03:51:50  avi
 * options -A -t
 *
 * Revision 1.2  2006/11/04  06:24:19  avi
 * Solaris 10
 * option -l (output cluster center of mass in floating indices)
 *
 * Revision 1.1  2005/05/20  06:31:53  avi
 * Initial revision
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define UCHAR		unsigned char
#define USHORT		unsigned short
#define MAXR		4096		/* regions */
#define MR		24		/* gfc regions */
#define MAXL		256		/* string length */
#define MAXD		65536		/* maximum image dimension (dictated by XYZN) */
#define MEMB		128		/* region buffer memory increment */

typedef struct {
	USHORT	ix;
	USHORT	iy;
	USHORT	iz;
	USHORT	nn;
} XYZN;

typedef struct {
	float		u1, u2;
	float		fndex[3];
	int		count;
	int		ne;
	int		nxyzn;
	XYZN		*xyzn;
} REGION;

/********************/
/* global variables */
/********************/
static char		rcsid[] = "$Id: cluster_4dfp.c,v 1.14 2012/06/14 06:01:42 avi Exp $";
static REGION		*region;
static int		nr = 0, maxr = 0;
static char		program[MAXL];
static IFH		ifh;
static int		dimension, xdim, ydim, zdim, imgdim[4], orient, isbig;
static char		control = '\0';
static float		voxdim[3];
static float		*img1;
static unsigned int	*imgo;
static float		thresh = 0.0;		/* definition of cluster boundary */
static int		absflag = 0;		/* threshold on absolute value */
static int		sizecrit = 0;
static int		verbose = 0;

/*************/
/* externals */
/*************/
void flipx (float *imgf, int *pnx, int* pny, int *pnz);		/* cflip_4dfp.c */
void flipy (float *imgf, int *pnx, int* pny, int *pnz);		/* cflip_4dfp.c */
void flipz (float *imgf, int *pnx, int* pny, int *pnz);		/* cflip_4dfp.c */

/*******************************/
/* perform analyzeto4dfp flips */
/*******************************/
void analyzeto4dfp_flip (char *imgroot) {
	switch (orient) { 
		case 4:	flipx (img1, imgdim+0, imgdim+1, imgdim+2);
		case 3:	flipz (img1, imgdim+0, imgdim+1, imgdim+2);
		case 2:	flipy (img1, imgdim+0, imgdim+1, imgdim+2);
			break;					
		default:
			fprintf (stderr, "%s: invalid %s orientation (%d)\n", program, imgroot, orient);
			exit (-1);
			break;
	}
}

void alloc_region (void) {
	int	maxr0;

	maxr0 = maxr;
	maxr += MAXR;
	if (!maxr0) {
		if (!(region = (REGION *) malloc  (MAXR * sizeof (REGION)))) errm (program);
	} else {
		if (!(region = (REGION *) realloc (region, maxr * sizeof (REGION)))) errm (program);
	}
	memset (region + maxr0, '\0', MAXR * sizeof (REGION));
	printf ("alloc_region: maxr=%d\n", maxr);
}

static int above_thresh (float val) {
	if (absflag) {
		return (fabs (val) > thresh) ? 1 : 0;
	} else {
		return (val > thresh) ? 1 : 0;
	}
}

static int rcompare (const void *ptr1, const void *ptr2) {
	REGION *r1, *r2;

	r1 = (REGION *) ptr1;
	r2 = (REGION *) ptr2;
	return r2->count - r1->count;
}

/************************************/
/* compute region mean and variance */
/************************************/
void getmean () {
	XYZN		xyzn;
	double		q, u1, u2;
	int		ir, i, m, iv;

	for (ir = 0; ir < nr; ir++) {
		for (u2 = u1 = m = 0; m < region[ir].count; m++) {
			xyzn = region[ir].xyzn[m]; iv = xyzn.ix + xdim * (xyzn.iy + ydim * xyzn.iz);
			q = img1[iv];
			u1 += q;
			u2 += q*q;
		}
		u1 /= region[ir].count; u2 /= region[ir].count; u2 -= u1*u1;
		u2 = (u2 < 0.) ? 0.0 : u2;
		region[ir].u1 = u1;
		region[ir].u2 = (region[ir].count > 1) ? u2 : 0.0;
	}
}

void convert_to_ROI () {
	XYZN		xyzn;
	int		ir, m, iv;

	for (iv = 0; iv < dimension; iv++) img1[iv] = 0.;
	for (ir = 0; ir < nr; ir++) {
		for (m = 0; m < region[ir].count; m++) {
			xyzn = region[ir].xyzn[m]; iv = xyzn.ix + xdim * (xyzn.iy + ydim * xyzn.iz);
			img1[iv] = ir + 2;
		}
	}
}

static void assign (int ir, XYZN xyzn) {
	int		i, j, k, iv;
	int		debug = 0;

	if (!(region[ir].count < region[ir].nxyzn)) {
		region[ir].nxyzn += MEMB;
		if (region[ir].xyzn) {
			region[ir].xyzn = (XYZN *) realloc (region[ir].xyzn, region[ir].nxyzn * sizeof (XYZN));
		} else {
			region[ir].xyzn = (XYZN *)  malloc (region[ir].nxyzn * sizeof (XYZN));
		}
		if (!region[ir].xyzn) errm (program);
	}
	region[ir].xyzn[region[ir].count++] = xyzn;
	if (debug) {
		iv = xyzn.ix + xdim*(xyzn.iy + ydim*xyzn.iz);
		printf ("assign ir=%6d count=%4d examined=%4d indices=%4d%4d%4d value=%10.4f\n",
		ir, region[ir].count, region[ir].ne, xyzn.ix + 1, xyzn.iy + 1, xyzn.iz + 1, img1[iv]);
	}
}

/**************************/
/* list region statistics */
/**************************/
static void listreg (FILE *fp) {
	int		ir, it;

	fprintf (fp, "region   voxels      mean      s.d.\n");
	for (ir = 0; ir < nr; ir++) {
		it = verbose || region[ir].count >= sizecrit;
		if (it) fprintf (fp, "%-5d%10d%10.4f%10.4f\n", ir + 1, region[ir].count, region[ir].u1, sqrt(region[ir].u2));
	}
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

int main (int argc, char *argv[]) {
	FILE 		*imgfp, *outfp;
	char            imgroot[MAXL], imgfile[MAXL], hdrfile[MAXL], outroot[MAXL], outfile[MAXL], trailer[MAXL] = "";
	char		string[MAXL], command[MAXL], *ptr;

	USHORT		ix, iy, iz;
	XYZN		xyzn;
	int		iv, iv0, jv, ir;
	int		c, i, j, k, m;
	int		nzer = 0, nvox;
	float		amax, amin;

/*********/
/* flags */
/*********/
	int		vflag = 0, iframe = 1;
	int		ROI_flag = 0;
	int		make_list = 0;
	int		status = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) if (*argv[i] == '-') {
		strcpy (string, argv[i]); ptr = string;
		while (c = *ptr++) switch (c) {
			case 'l': make_list++;			break;
			case 'v': verbose++;			break;
			case 'A': absflag++;			break;
			case 'R': ROI_flag++;			break;
			case 'f': iframe = atoi (ptr); vflag++;	*ptr = '\0'; break;
			case 't': thresh = atof (ptr);		*ptr = '\0'; break;
			case 'a': strncpy (trailer, ptr, MAXL);	*ptr = '\0'; break;
			case 'n': sizecrit = atoi (ptr);	*ptr = '\0'; break;
			case '@': control = *ptr++;		*ptr = '\0'; break;
		}
	} else switch (k) {
		case 0: getroot (argv[i], imgroot); k++;  break;
	}
	if (k < 1) {
		fprintf (stderr, "Usage:\t%s <(4dfp) root>\n", program);
		fprintf (stderr, "e.g.,\t%s my_timage -At3.5\n", program);
		fprintf (stderr, "\toption\n");
		fprintf (stderr, "\t-n<int>	zero out clusters with voxel count below specified criterion\n");
		fprintf (stderr, "\t-f<int>	address specified volume (counting from 1) of multi-volume stack\n");
		fprintf (stderr, "\t-t<flt>	specify image value threshold (default = 0)\n");
		fprintf (stderr, "\t-a<str>	append specified string (preceded by \"_\") to all output filenames\n");
		fprintf (stderr, "\t-@<b|l>	output big or little endian (default input endian)\n");
		fprintf (stderr, "\t-A	apply threshold test to image absolute value\n");
		fprintf (stderr, "\t-R	convert clusters to (fidl compliant) ROI image\n");
		fprintf (stderr, "\t-l	create list file of region center of mass indices\n");
		fprintf (stderr, "\t-v	verbose mode\n");
		fprintf (stderr, "N.B.:\tif the input is multi-volume %s addresses the first volume by default\n", program);
		fprintf (stderr, "N.B.:\t-n (without -R) creates output 4dfp image with fileroot trailer \"clus\"\n");
		fprintf (stderr, "N.B.:\t-R              creates output 4dfp image with fileroot trailer \"ROI\"\n");
		fprintf (stderr, "N.B.:\t-l center of mass indices can be converted to atlas coordinates using index2atl -af\n");
		exit (1);
	}

/*****************************/
/* get 4dfp input dimensions */
/*****************************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0
	||  Getifh (imgfile, &ifh)) errr (program, imgfile);
	if (!control) control = (isbig) ? 'b' : 'l';
	dimension = imgdim[0] * imgdim[1] * imgdim[2];
	xdim = imgdim[0]; ydim = imgdim[1]; zdim = imgdim[2];

/*************************************************/
/* check image index limits against XYZN typedef */
/*************************************************/
	if (xdim > MAXD || ydim > MAXD || zdim > MAXD) {
		fprintf (stderr, "%s: %s dimension exceeds %d\n", program, imgroot, MAXD);
		exit (-1);
	}
	if (!(img1 = (float *) malloc (dimension * sizeof (float)))) errm (program);
	if (!(imgo = (unsigned int *) calloc (dimension,  sizeof (unsigned int)))) errm (program);

/************************/
/* read selected volume */
/************************/
	printf ("Reading: %s volume %d\n", imgfile, iframe);
	if (!(imgfp = fopen (imgfile, "rb"))
	|| fseek (imgfp, (long) (iframe - 1) * dimension * sizeof (float), SEEK_SET)
	|| eread  (img1, dimension, isbig, imgfp)
	|| fclose (imgfp)) errr (program, imgfile);
	analyzeto4dfp_flip (imgroot);

/****************************/
/* zero subthreshold voxels */
/****************************/
	for (i = 0; i < dimension; i++) if (!above_thresh (img1[i])) img1[i] = 0.0;

/**********************************************/
/* initially allocate and clear region buffer */
/**********************************************/
	alloc_region ();

	iv0 = ir = 0;
/***************************/
/* find contiguous regions */
/***************************/
	while (1) {
		for (i = 0; i < dimension; i++) {
			iv = (iv0 + i) % dimension;
			if (!imgo[iv] && above_thresh (img1[iv])) break;
		}
		if (i == dimension) break;
		iv0 = iv;

		imgo[iv] = ir + 1;
		k = iv;
		j = dimension;
		j /= zdim; iz = k / j; k -= iz * j;
		j /= ydim; iy = k / j; k -= iy * j;
		j /= xdim; ix = k / j; k -= ix * j;
		xyzn.ix = ix; xyzn.iy = iy; xyzn.iz = iz; xyzn.nn = 0;
		assign (ir, xyzn);

		do {
			m = region[ir].ne++;
			iv = region[ir].xyzn[m].ix + xdim * (region[ir].xyzn[m].iy + ydim * region[ir].xyzn[m].iz);
			for (j = 0; j < 6; j++) {
				xyzn = region[ir].xyzn[m];
				switch (j) {
				case 0:	k =  -1;	if (xyzn.ix < 1)	continue; xyzn.ix--; break;
				case 1:	k *= -1;	if (xyzn.ix > xdim - 2)	continue; xyzn.ix++; break;
				case 2:	k *= -xdim;	if (xyzn.iy < 1)	continue; xyzn.iy--; break;
				case 3:	k *= -1;	if (xyzn.iy > ydim - 2)	continue; xyzn.iy++; break;
				case 4:	k *= -ydim;	if (xyzn.iz < 1)	continue; xyzn.iz--; break;
				case 5:	k *= -1;	if (xyzn.iz > zdim - 2)	continue; xyzn.iz++; break;
				}
				jv = iv + k;
				if (!imgo[jv] && above_thresh (img1[jv]) && img1[jv]*img1[iv0] > 0.) {
					imgo[jv] = ir + 1;
					assign (ir, xyzn);
				}
			}
			if (verbose && !((ir + 1) % 100)) {
				xyzn = region[ir].xyzn[m];
				printf ("region=%8d count=%8d unexamined=%8d coords=%4d%4d%4d\n",
				ir + 1, region[ir].count, region[ir].count - m, xyzn.ix, xyzn.iy, xyzn.iz);
			}
		} while (region[ir].count > region[ir].ne);
		if (++ir == maxr) alloc_region ();
	}

/************************/
/* sort regions by size */
/************************/
	qsort ((void *) region, ir, sizeof (REGION), rcompare);
	nr = ir;
	getmean ();
	listreg (stdout);

/****************************************************************/
/* zero out regions that do not meet a specified size criterion */
/****************************************************************/
	for (ir = 0; ir < nr; ir++) if (region[ir].count < sizecrit) {
		for (m = 0; m < region[ir].count; m++) {
			xyzn = region[ir].xyzn[m]; iv = xyzn.ix + xdim * (xyzn.iy + ydim * xyzn.iz);
			img1[iv] = 0.0;
		}
		nzer++;
		region[ir].count = 0;
	}

/**************************************/
/* count voxels in remaining clusters */
/**************************************/
	for (nvox = ir = 0; ir < nr; ir++) nvox += region[ir].count;

	if (!vflag) {
		sprintf (outroot, "%s_%s", imgroot, ((ROI_flag) ? "ROI" : "clus"));
	} else {
		sprintf (outroot, "%s_vol%d_%s", imgroot, iframe, ((ROI_flag) ? "ROI" : "clus"));
	}
	if (strlen (trailer)) {
		strcat (outroot, "_"); strcat (outroot, trailer);
	}

if (make_list) {
/***********************************/
/* compute regional center of mass */
/***********************************/
	sprintf (outfile, "%s.dat", outroot);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);
	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
	fprintf (outfp, "#%9s%10s%10s%10s%10s%10s%10s\n", "index_x", "index_y", "index_z", "region", "voxels", "mean", "s.d.");
	for (ir = 0; ir < nr; ir++) if (region[ir].count) {
		for (m = 0; m < region[ir].count; m++) {
			xyzn = region[ir].xyzn[m];
			region[ir].fndex[0] += xyzn.ix + 1;
			region[ir].fndex[1] += xyzn.iy + 1;
			region[ir].fndex[2] += xyzn.iz + 1;
		}
		for (k = 0; k < 3; k++) region[ir].fndex[k] /= region[ir].count;
		fprintf (outfp, "%10.4f%10.4f%10.4f    #%5d%10d%10.4f%10.4f\n",
			region[ir].fndex[0], region[ir].fndex[1], region[ir].fndex[2], ir + 1,
			region[ir].count, region[ir].u1, sqrt(region[ir].u2));
	}
	fprintf (outfp, "#total voxel count = %d\n", nvox);
	if (fclose (outfp)) errw (program, outfile);
}

	if (!sizecrit && !ROI_flag) goto FREE;
	if (ROI_flag) convert_to_ROI ();

/*****************************/
/* compute voxel value range */
/*****************************/
	amin = FLT_MAX; amax = -FLT_MAX;
	for (i = 0; i < dimension; i++) {
		if (img1[i] < amin) amin = img1[i];
		if (img1[i] > amax) amax = img1[i];
	}
	amin -= 0.4999; amax += 0.4999;

/************************/
/* write modified image */
/************************/
	analyzeto4dfp_flip (imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);
	fprintf (stdout, "Writing: %s\n", outfile);
	if (!(imgfp = fopen (outfile, "wb"))
	|| ewrite (img1, dimension, control, imgfp)
	|| fclose (imgfp)) errw (program, imgfile);

/*******/
/* ifh */
/*******/
	ifh.matrix_size[3] = 1;		/* all output images are 1 volume */
	if (Writeifh (program, outroot, &ifh, control)) errw (program, outroot);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s -r%dto%d", outroot, (int) amin, (int) amax);
	printf ("%s\n", command);
	status |= system (command);

/*******/
/* rec */
/*******/
	startrece (outfile, argc, argv, rcsid, control);
	sprintf (string, "Threshold = %.4f", thresh);				printrec (string);
	if (absflag) printrec (" applied to absolute value");			printrec ("\n");
	sprintf (string, "Starting number of clusters = %d\n", nr);		printrec (string);
	sprintf (string, "Cluster size criterion = %d voxels\n", sizecrit);	printrec (string);
	sprintf (string, "%d clusters zeroed\n", nzer);				printrec (string);
	sprintf (string, "Final number of clusters = %d\n", nr - nzer);		printrec (string);
	sprintf (string, "Final voxel count = %d\n", nvox);			printrec (string);
	catrec (imgfile);
	endrec ();

FREE:	free (img1);
	free (imgo);
	for (ir = 0; ir < nr; ir++) free (region[ir].xyzn);
	exit (status);
}
