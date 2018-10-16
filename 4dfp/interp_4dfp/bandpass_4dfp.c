/* $Header: /data/petsun4/data1/src_solaris/interp_4dfp/RCS/bandpass_4dfp.c,v 1.10 2013/02/05 23:19:19 avi Exp $*/
/* $Log: bandpass_4dfp.c,v $
 * Revision 1.10  2013/02/05  23:19:19  avi
 * default nskip 4 -> 0
 * always allocate scratch arrays a,b and used new subroutine butt1dbs_() instead of butt1db_()
 * remove option -c
 * add option -r (retain linear trend)
 *
 * Revision 1.9  2011/02/06  08:27:08  avi
 * option -B
 *
 * Revision 1.8  2007/11/20  03:16:10  avi
 * linux compatible (remove f_init() and f_exit())
 *
 * Revision 1.7  2006/09/25  16:53:34  avi
 * Solaris 10
 *
 * Revision 1.6  2006/08/07  03:25:58  avi
 * correct 1.e-37 test
 *
 * Revision 1.5  2004/12/29  00:42:05  avi
 * conc functionality
 * unconfuse low/high half-frequencies vs. low/high pass
 *
 * Revision 1.4  2004/11/22  22:15:59  rsachs
 * Installed 'setprog'. Replaced 'Get4dfpDimN' with 'get_4dfp_dimo'.
 *
 * Revision 1.3  2004/08/31  04:55:38  avi
 * improved algorithm for computing DC shift and linear trend
 * -a now controls only retaining DC shift
 *
 * Revision 1.2  2004/05/26  20:55:42  avi
 * improve rec file filter characteristic listing
 *
 * Revision 1.1  2004/05/26  20:31:21  avi
 * Initial revision
 *
 * Revision 1.4  2002/06/26  05:37:44  avi
 * better usage
 *
 * Revision 1.3  2002/06/26  05:14:20  avi
 * -a (keepDC) option
 *
 * Revision 1.2  2002/06/25  05:32:02  avi
 * impliment Butterworth highpasss filter
 *
 * Revision 1.1  2002/06/25  02:39:29  avi
 * Initial revision
 **/
/*************************************************************
Purpose:	time series filter  
Date:		6/22/02
Author:		Avi Snyder
**************************************************************/
#include	<stdlib.h>
#include	<stdio.h>
#include	<math.h>
#include	<stdlib.h>
#include	<unistd.h>		/* R_OK */
#include	<string.h>
#include	<assert.h>
#include	<Getifh.h>
#include	<endianio.h>
#include        <rec.h>
#include 	<conc.h>

#define MAXL		256
#define MARGIN		16
#define ORDER_LO	0
#define ORDER_HI	0
#define NSKIP		0
#define TRAILER		"bpss"

/*************/
/* externals */
/*************/
extern int npad_ (int *tdim, int *margin);	/* FORTRAN librms */
extern void  butt1dbs_ (float *data,          int *n, float *delta, float *fhalf_lo, int *iorder_lo, float *fhalf_hi, int *iorder_hi, float *a, float *b);
extern void pbutt1dcs_ (float *data, int *n0, int *n, float *delta, float *fhalf_lo, int *iorder_lo, float *fhalf_hi, int *iorder_hi, float *a, float *b);

void usage (char* program) {
	printf ("Usage:\t%s <(4dfp|conc) input> <TR_vol>\n", program);
	printf ("e.g.:\t%s qst1_b1_rmsp_dbnd_xr3d[.4dfp.img] 2.36 -bl0.01 -ol1 -bh0.15 -oh2\n", program);
	printf ("\toption\n");
	printf ("\t-b[l|h]<flt>\tspecify low end or high end half frequency in hz\n");
	printf ("\t-o[l|h]<int>\tspecify low end or high end Butterworth filter order\n");
	printf ("\t-n<int>\tspecify number of pre-functional frames (default = %d)\n", NSKIP);
	printf ("\t-t<str>\tchange output filename trailer (default=\"_%s\")\n", TRAILER);
	printf ("\t-a\tretain DC (constant) component\n");
	printf ("\t-r\tretain linear trend\n");
	printf ("\t-E\tcode undefined voxels as 1.e-37\n");
	printf ("\t-B\tcompute gain using correct Butterworth formula (default squared Butterworth gain)\n");
	printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
	printf ("N.B.:\tinput conc files must have extension \"conc\"\n");
	printf ("N.B.:\toutput filnename root is <input>_<trailer>\n", TRAILER);
	printf ("N.B.:\tomitting low  end order specification disables high pass component\n");
	printf ("N.B.:\tomitting high end order specification disables low  pass component\n");
	exit (1);
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[] = "$Id: bandpass_4dfp.c,v 1.10 2013/02/05 23:19:19 avi Exp $";
int main (int argc, char *argv[]) {
	char		imgroot[MAXL], imgfile[MAXL], outroot[MAXL], outfile[MAXL];
	char		trailer[MAXL] = TRAILER;
	FILE		*fp;
	CONC_BLOCK	conc_block;		/* conc i/o control block */
	IFH		ifh;			/* for non-conc input */

	char		control = '\0';
	int		c, i, j, k, ix, iy, iz, ivox, orient, isbig, osbig;
	int		nfile, ifile;
	int		imgdim[4], vdim, dimension, tdim, tdim_pad, margin = MARGIN, padlen;
	int		order_lo = ORDER_LO, order_hi = ORDER_HI, nskip = NSKIP;

	float		fhalf_lo = 0.0, fhalf_hi = 0.0, TR_vol;
	float		*imgt, *imga;
	float		*a, *b;		/* scratch buffers used by Butterworth filters */
	float		voxdim[3];
	short		*mask;		/* set to 1 at undefined (1.e-37) voxels */
	short		*imgn;		/* denominator for average volume */

/****************************/
/* linear trend computation */
/****************************/
	float		*x, sy, sxy, sxx, a0, a1;
	float		*tpad, q;

/***********/
/* utility */
/***********/
	char		*ptr, program[MAXL], command[MAXL];

/*********/
/* flags */
/*********/
	int		status;
	int		conc_flag = 0;
	int		keepDC = 0;
	int		keepramp = 0;
	int		E_flag = 0;
	int		B_flag = 0;

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
/******************************/
/* get command line arguments */
/******************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'a': keepDC++;			break;
				case 'r': keepramp++;			break;
				case 'E': E_flag++;			break;
				case 'B': B_flag++;			break;
				case '@': control = *ptr++;		*ptr = '\0'; break;
				case 'n': nskip = atoi (ptr);		*ptr = '\0'; break;
				case 't': strcpy (trailer, ptr);	*ptr = '\0'; break;
				case 'o': switch (*ptr++) {
					case 'l': order_lo = atoi (ptr);	break;
					case 'h': order_hi = atoi (ptr);	break;
					default:  usage (program);		break;
					}				*ptr = '\0'; break;
				case 'b': switch (*ptr++) {
					case 'l': fhalf_lo = atof (ptr);	break;
					case 'h': fhalf_hi = atof (ptr);	break;
					default:  usage (program);      	break;
					}				*ptr = '\0'; break;
			}
		} else switch (k) {
		 	case 0: getroot (argv[i], imgroot);
				conc_flag = (strstr (argv[i], ".conc") == argv[i] + strlen (imgroot));
								k++; break;
		 	case 1: TR_vol = atof (argv[i]);	k++; break;
		}
	}
	if (k < 2) usage (program);
	printf ("fhalf_lo %.4f order_lo %d fhalf_hi %.4f order_hi %d\n", fhalf_lo, order_lo, fhalf_hi, order_hi);
	if (fhalf_lo <= 0.0 || order_lo < 0) order_lo = 0;
	if (fhalf_hi <= 0.0 || order_hi < 0) order_hi = 0;

	if (conc_flag) {
		conc_init (&conc_block, program);
		conc_open (&conc_block, imgroot);
		for (k = 0; k < 3; k++) imgdim[k] = conc_block.imgdim[k];
		vdim = conc_block.vdim;
		nfile = conc_block.rnfile;
		isbig = conc_block.isbig;
	} else {
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgfile);
		if (Getifh (imgfile, &ifh)) errr (program, imgroot);
		for (vdim = 1, k = 0; k < 3; k++) vdim *= imgdim[k];
		nfile = 1;
		sprintf (outroot, "%s_%s", imgroot, trailer);
		sprintf (outfile, "%s.4dfp.img", outroot);
	}

	if (!control) control = (isbig) ? 'b' : 'l';
	if (conc_flag) conc_newe (&conc_block, trailer, control);
	sprintf (outroot, "%s_%s", imgroot, trailer);
	if (!(mask = (short *) calloc (vdim, sizeof (short)))) errm (program);
	if (!(imgn = (short *) calloc (vdim, sizeof (short)))) errm (program);	/* denominator for imga */
	if (!(imga = (float *) calloc (vdim, sizeof (float)))) errm (program);	/* average volume */

/***********************/
/* loop on input files */
/***********************/
for (ifile = 0; ifile < nfile; ifile++) {	/* begin ifile loop */
	if (conc_flag) imgdim[3] = conc_block.nvol[ifile];
	dimension = vdim * imgdim[3];
	if (!(imgt = (float *) malloc (dimension * sizeof (float)))) errm (program);
	if (conc_flag) strcpy (imgfile, conc_block.imgfile0[ifile]);
	printf ("Reading: %s\n", imgfile);
	if (!(fp = fopen (imgfile, "rb")) || eread (imgt, dimension, isbig, fp)
	|| fclose (fp)) errr (program, imgfile);

/********************************/
/* allocate time series buffers */
/********************************/
	tdim = imgdim[3] - nskip;
	tdim_pad = npad_ (&tdim, &margin);
	padlen = tdim_pad - tdim;
	printf ("original time series length %d padded to %d\n", tdim, tdim_pad);
	if (!(tpad = (float *) malloc (tdim_pad * sizeof (float)))) errm (program);
	if (!(x =    (float *) malloc (tdim     * sizeof (float)))) errm (program);
	for (i = 0; i < tdim; i++) x[i] = -1. + 2.*i/(tdim - 1);
	sxx = ((float) tdim*(tdim+1)) / (3.*(tdim-1));

	if (!(a = (float *) calloc (tdim_pad/2 + 1, sizeof (float)))) errm (program);
	if (!(b = (float *) calloc (tdim_pad/2 + 1, sizeof (float)))) errm (program);
/*********************************/
/* process all voxels of one run */
/*********************************/
	for (ivox = 0; ivox < vdim; ivox++) mask[ivox] = 0;
	ivox = 0;
	printf ("processing slice");
	for (iz = 0; iz < imgdim[2]; iz++) {printf(" %d", iz + 1); fflush (stdout);
	for (iy = 0; iy < imgdim[1]; iy++) {
	for (ix = 0; ix < imgdim[0]; ix++) {
		for (k = i = 0; i < imgdim[3]; i++) {
			if ((imgt + ivox)[i * vdim] == (float) 1.e-37) mask[ivox] = 1;
			if (i >= nskip) tpad[k++] = (imgt + ivox)[i * vdim];
		}

/******************************/
/* remove DC and linear trend */
/******************************/
		for (sy = sxy = k = 0; k < tdim; k++) {
			sy  += tpad[k];
			sxy += tpad[k]*x[k];
		}
		a0 = sy/tdim;
		a1 = (keepramp) ? 0. : sxy/sxx;
		for (k = 0; k < tdim; k++) tpad[k] -= (a0 + x[k]*a1);

/*********************************/
/* circularly connect timeseries */
/*********************************/
		q = (tpad[0] - tpad[tdim - 1]) / (padlen + 1);
		for (j = 1, k = tdim; k < tdim_pad; k++) tpad[k] = tpad[tdim - 1] + q*j++;

/**********/
/* filter */
/**********/
		if (B_flag) {
			pbutt1dcs_ (tpad, &tdim, &tdim_pad, &TR_vol, &fhalf_lo, &order_lo, &fhalf_hi, &order_hi, a, b);
		} else {
			 butt1dbs_ (tpad,        &tdim_pad, &TR_vol, &fhalf_lo, &order_lo, &fhalf_hi, &order_hi, a, b);
		}

/********************************************************************************/
/* force unpadded timeseries to zero mean and put fltered results back in image */
/********************************************************************************/
		for (q = k = 0; k < tdim; k++) q += tpad[k];
		q /= tdim;
		for (k = 0, i = nskip; i < imgdim[3]; i++) (imgt + ivox)[i * vdim] = tpad[k++] - q;

		for (i = 0; i < imgdim[3]; i++) {
			if (E_flag && mask[ivox]) {				(imgt + ivox)[i * vdim] = 1.e-37;
			} else {
				if (i <  nskip && (conc_flag || !keepDC))	(imgt + ivox)[i * vdim] -= a0;
				if (i >= nskip && !conc_flag &&  keepDC)	(imgt + ivox)[i * vdim] += a0;
			}
		}
		if (!E_flag || !mask[ivox]) {
			imga[ivox] += a0;
			imgn[ivox]++;
		}
		ivox++;
	}}}
	printf("\n");

/**************************************/
/* write bandpass filtered 4dfp stack */
/**************************************/
	if (conc_flag) strcpy (outfile, conc_block.imgfile1[ifile]);
	printf ("Writing: %s\n", outfile);
	if (!(fp = fopen (outfile, "wb")) || ewrite (imgt, dimension, control, fp)
	|| fclose (fp)) errw (program, outfile);

	free (a); free (b);
	free (tpad); free (x); free (imgt);
}	/* end ifile loop */

/************************/
/* restore DC component */
/************************/
	switch (control) {
		case 'b': case 'B': osbig = 1; break;
		case 'l': case 'L': osbig = 0; break;
		default: osbig = CPU_is_bigendian(); break;
	}
	if (conc_flag && keepDC) {
		for (ivox = 0; ivox < vdim; ivox++) if (imgn[ivox]) imga[ivox] /= imgn[ivox];
		if (!(imgt = (float *) malloc (vdim * sizeof (float)))) errm (program);
		for (ifile = 0; ifile < nfile; ifile++) {
			printf ("Adding back DC %s frame", conc_block.imgfile1[ifile]);
			if (!(fp = fopen (conc_block.imgfile1[ifile], "r+b"))) errr (program, conc_block.imgfile1[ifile]);
			for (k = 0; k < conc_block.nvol[ifile]; k++) {printf(" %d", k + 1); fflush (stdout);
				if (fseek (fp, (long) k*vdim*sizeof (float), SEEK_SET)
				||  eread (imgt, vdim, osbig, fp))
					errr (program, conc_block.imgfile1[ifile]);
				for (ivox = 0; ivox < vdim; ivox++) {
					if (E_flag && imgt[ivox] == (float) 1.e-37) continue;
					imgt[ivox] += imga[ivox];
				}
				if (fseek (fp, (long) k*vdim*sizeof (float), SEEK_SET)
				||  ewrite (imgt, vdim, control, fp))
					errw (program, conc_block.imgfile1[ifile]);
			}
			printf ("\n"); fflush (stdout);
			if (fclose (fp)) errw (program, conc_block.imgfile1[ifile]);
		}
		free (imgt);
	}

/***************/
/* ifh and hdr */
/***************/
	if (conc_flag) {
		conc_ifh_hdr_rec (&conc_block, argc, argv, rcsid);
		conc_free (&conc_block);
		sprintf (outfile, "%s.conc", outroot);
		sprintf (imgfile, "%s.conc", imgroot);
		printf ("Writing: %s\n", outfile);
	} else {
		if (Writeifh (program, outfile, &ifh, control)) errw (program, outroot);
		sprintf (command, "ifh2hdr %s", outroot); status |= system (command);
	}

/*******/
/* rec */
/*******/
	startrece (outfile, argc, argv, rcsid, control);
	sprintf (command, "TR_vol (sec) %.4f\n", TR_vol);					printrec (command);
	if (!conc_flag) {
		sprintf (command, "Original time series length %d padded to %d\n", tdim, tdim_pad);
		printrec (command);
	}
	if (keepDC) {
		sprintf (command, "DC shift retained\n");
	} else {
		sprintf (command, "DC shift removed\n");
	}											printrec (command);
	if (keepramp) {
		sprintf (command, "Linear trend retained\n");
	} else {
		sprintf (command, "Linear trend removed\n");
	}											printrec (command);
	if (order_lo) {
		sprintf (command, "High pass fhalf=%.4f (Hz)  order=%d\n", fhalf_lo, order_lo);	printrec (command);
	}
	if (order_hi) {
		sprintf (command, "Low  pass fhalf=%.4f (Hz)  order=%d\n", fhalf_hi, order_hi);	printrec (command);
	}
	catrec (imgfile);
	endrec ();

	free (mask); free (imgn); free (imga);
	exit (0);
}
