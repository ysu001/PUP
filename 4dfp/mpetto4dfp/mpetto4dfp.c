/*$Id: mpetto4dfp.c,v 1.12 2011/08/07 03:59:21 avi Exp $*/
/*$Log: mpetto4dfp.c,v $
 * Revision 1.12  2011/08/07  03:59:21  avi
 * apply calibration_factor to PET data
 * convert CT data
 *
 * Revision 1.11  2008/07/30  03:24:09  avi
 * option -w
 *
 * Revision 1.10  2008/02/22  03:46:15  avi
 * X86 and linux compliant
 *
 * Revision 1.9  2008/01/09  04:08:50  avi
 * option -o; finish all switch looops with default
 *
 * Revision 1.8  2007/08/16  23:06:20  avi
 * enable conversion of file type 8 (Crystal efficiency data)
 *
 * Revision 1.7  2007/05/26  01:18:59  avi
 * get slice thickness form axial_crystal_pitch instead of axial_plane_size
 *
 * Revision 1.6  2007/04/28  00:36:30  avi
 * endianio.h Getifh.h
 * not yet CPU endian compliant
 *
 * Revision 1.5  2006/03/21  04:47:22  avi
 * comment out checking for file type == 5 (image) to work with transmission images
 *
 * Revision 1.4  2006/03/07  06:02:10  avi
 * debug max/min reporting
 *
 * Revision 1.3  2006/03/07  05:45:39  avi
 * in rec translate enumerated fields
 *
 * Revision 1.2  2006/03/07  02:37:30  avi
 * basic inclusion of scan-specific info in rec file
 *
 * Revision 1.1  2006/03/07  01:44:18  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MAXL	256

typedef struct {
	int		event_type;
	int		gate;
	int		bed;
	float		bed_offset;		/* cm */
	float		ending_bed_offset;	/* cm */
	float		vertical_bed_offset;	/* cm */
	unsigned long	data_file_pointer[2];
	float		frame_start;		/* sec */
	float		frame_duration;		/* sec */
	float		scale_factor;
	float		minimum, maximum;
} FRAMEDAT;

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;

	if (ptr = strchr (string, '#')) *ptr = '\0';
	i = m = 0;
	while (m < maxp) {
		while (!isgraph ((int) string[i]) && string[i]) i++;
		if (!string[i]) break;
		srgv[m++] = string + i;
		while (isgraph ((int) string[i])) i++;
		if (!string[i]) break;
		string[i++] = '\0';
	}
	return m;
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

extern void flipx (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */
extern void flipy (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */
extern void flipz (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */

/***********/
/* globals */
/***********/
static FRAMEDAT	*framedat;
static int	imgdim[4];

static char	*file_type[] = {
		"Unknown data file type",
		"List mode data file",
		"Sinogram data file",
		"Normalization data file",
		"Attenuation correction data file",
		"Image data file Blank data file",
		"Mu map data file",
		"Scatter correction data file",
		"Crystal efficiency data",
		"Crystal interference correction",
		"Transaxial geometric correction"};

static char	*acquisition_mode[] = {
		"Unknown acquisition mode",
		"Blank acquisition",
		"Emission acquisition",
		"Dynamic acquisition",
		"Gated acquisition",
		"Continuous bed motion acquisition",
		"Singles transmission acquisition",
		"Windowed coincidence transmission acquisition",
		"Non-windowed coincidence transmission acquisition"};

static char	*bed_motion[] = {
		"Static or unknown bed motion",
		"Continuous bed motion",
		"Multiple bed positions, i.e. step and shoot"};

static char	*tx_src_type[] = {
		"Unknown TX source type",
		"TX point source",
		"TX line source"};

static char	*data_type[] = {
		"Unknown data type",
		"Byte (8-bits) data type",
		"2-byte integer - Intel style",
		"4-byte integer - Intel style",
		"4-byte float - Intel style",
		"4-byte float - Sun style",
		"2-byte integer - Sun style",
		"4-byte integer - Sun style"};

static char	*xyz_filter[] = {
		"No filter",
		"Ramp filter (backprojection) or no filter",
		"First-order Butterworth window",
		"Hanning window",
		"Hamming window",
		"Parzen window",
		"Shepp filter",
		"Second-order Butterworth window"};

static char	*rebinning_type[] = {
		"Unknown, or no, algorithm type",
		"Full 3D binning (span and ring difference)",
		"Single-Slice Rebinning",
		"Fourier Rebinning"};

static char	*recon_algorithm[] = {
		"Unknown, or no, algorithm type",
		"Filtered Backprojection",
		"OSEM2d",
		"?", "?", "?",
		"OSEM3D followed by MAP"};

static char	*deadtime_correction_applied[] = {
		"No deadtime correction applied",
		"Global estimate based on singles",
		"CMS estimate based on singles"};

static char	*normalization_applied[] = {
		"No normalization applied",
		"Point source inversion",
		"Point source component based",
		"Cylinder source inversion",
		"Cylinder source component based"};

static char	*attenuation_applied[] = {
		"No attenuation applied",
		"Point source in windowed TX coincidence",
		"Point source singles based TX",
		"Segmented point source in TX coincidence",
		"Segmented point source singles based TX",
		"Calculated by geometry",
		"Non-positron source singles based TX",
		"Point source in non-windowed TX coincidence"};

static char	*scatter_correction[] = {
		"No scatter correction applied",
		"Fit of emission tail",
		"Monte Carlo of emission and transmission data",
		"Direct calculation from analytical formulas",
		"Model-based scatter for singles TX",
		"TX off-window windowed coincidence subtraction"};

static char	*calibration_units[] = {
		"Unknown calibration units",
		"nCi/cc",
		"Bq/cc"};

static char	*dose_units[] = {
		"Unknown dose units",
		"mCi",
		"MBq"};

static char	*subject_orientation[] = {
		"Unknown subject orientation",
		"Feet first, prone",
		"Head first, prone",
		"Feet first, supine",
		"Head first, supine",
		"Feet first, right",
		"Head first, right",
		"Feet first, left",
		"Head first, left"};

static char	*subject_length_units[] = {
		"Unknown length units",
		"millimeters",
		"centimeters",
		"inches"};

static char	*subject_weight_units[] = {
		"Unknown weight units",
		"grams",
		"ounces",
		"kilograms",
		"pounds"};

static char	*event_type[] = {
		"Unknown event type",
		"Singles",
		"Prompt events (coincidences)",
		"Delay events",
		"Trues (prompts - delays)"};

static char	*keywords[] = {
		"version",
		"model",
		"institution",
		"study",
		"file_name",
		"file_type",
		"acquisition_mode",
		"bed_motion",
		"total_frames",
		"isotope",
		"isotope_half_life",
		"isotope_branching_fraction",
		"transaxial_crystals_per_block",
		"axial_crystals_per_block",
		"intrinsic_crystal_offset",
		"transaxial_blocks",
		"axial_blocks",
		"transaxial_crystal_pitch",
		"axial_crystal_pitch",
		"radius",
		"radial_fov",
		"src_radius",
		"src_cm_per_rev",
		"src_steps_per_rev",
/*		"tx_src_type",		*/
		"default_projections",
		"default_transaxial_angles",
		"crystal_thickness",
		"depth_of_interaction",
		"transaxial_bin_size",
		"axial_plane_size",
		"lld",
		"uld",
/*		"data_type",
		"data_order",		*/
		"span",
		"ring_difference",
/*		"number_of_dimensions",
		"x_dimension",
		"y_dimension",
		"z_dimension",
		"w_dimension",
		"x_filter",
		"y_filter",
		"z_filter",		*/
		"histogram_version",
		"rebinning_version",
/*		"rebinning_type",
		"recon_algorithm",	*/
		"recon_version",
/*		"deadtime_correction_applied",	*/
		"decay_correction_applied",
/*		"normalization_applied",	*/
		"normalization_filename",
/*		"attenuation_applied",		*/
		"attenuation_filename",
/*		"scatter_correction",		*/
		"scatter_version",
		"arc_correction_applied",
		"x_offset",
		"y_offset",
		"zoom",
/*		"pixel_size",
		"calibration_units",		*/
		"calibration_factor",
		"calibration_branching_fraction",
		"number_of_singles_rates",
		"investigator",
		"operator",
		"study_identifier",
		"scan_time",
		"injected_compound",
/*		"dose_units",			*/
		"dose",
		"injection_time",
		"injection_decay_correction",
		"subject_identifier",
		"subject_genus",
/*		"subject_orientation",
		"subject_length_units",		*/
		"subject_length",
/*		"subject_weight_units",		*/
		"subject_weight",
		"subject_phenotype",
		"study_model",
		"anesthesia",
		"analgesia",
		"other_drugs",
		"food_access",
		"water_access",
		"end_of_header"		/* termination signal */
	};
/*
		"frame",
		"event_type",
		"gate",
		"bed",
		"bed_offset",
		"ending_bed_offset",
		"vertical_bed_offset",
		"data_file_pointer",
		"frame_start",
		"frame_duration",
		"scale_factor",
		"deadtime_correction",
		"decay_correction",
		"prompts",
		"delays",
		"trues",
		"prompts_rate",
		"delays_rate"
*/

static char	program[MAXL];
static char	rcsid[] = "$Id: mpetto4dfp.c,v 1.12 2011/08/07 03:59:21 avi Exp $";

void write_framing (FILE *fpout) {
	int	i;

	fprintf (fpout, "%10s%10s%10s%12s%10s%10s\n", "frame", "start", "duration", "scale", "minimum", "maximum");
	for (i = 0; i < imgdim[3]; i++) {
		fprintf (fpout, "%10d%10.2f%10.2f%12.4e%10.4f%10.4f\n", i + 1, framedat[i].frame_start,
			framedat[i].frame_duration, framedat[i].scale_factor, framedat[i].minimum, framedat[i].maximum);
	}
}

int main (int argc, char *argv[]) {
	FILE		*fpimg, *fphdr, *fpout;
	char		imgroot[MAXL], outroot[MAXL] = "", hdrfile[MAXL];
	char		imgfile[MAXL], outfile[MAXL], recfile[MAXL];

/***********************/
/* microPET parameters */
/***********************/
	int		idata_type, idata_order, ifile_type, iacquisition_mode;
	float		calibration_factor = 0;

/****************/
/* image arrays */
/****************/
	short		*imgs;
	int		*imgd;
	float		*imgf;
	char		*iptr;	/* generic read pointer */
	float		voxdim[3], center[3], mmppix[3];
	int		vdim, bytepix;
	float		scale = 1., glmax = -FLT_MAX, glmin = FLT_MAX;
	char		control = '\0';

/***********/
/* utility */
/***********/
	char		*str, command[MAXL], string[MAXL], *srgv[MAXL];
	int 		c, i, k, m;
	float		q;

/*********/
/* flags */
/*********/
	int		modality;	/* -1: Unknown; 0: PET; 1: CT; 2: SPECT */
	int		status = 0;
	int		debug = 0;
	int		weifile = 0;
	int		CPUisbig, isbig, swab_flag = 0;
	int		xflag = 1, yflag = 0, zflag = 1;

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case 'd': debug++;		break;
				case 'w': weifile++;		break;
				case 'x': xflag ^= 1;		break;
				case 'y': yflag ^= 1;		break;
				case 'z': zflag ^= 1;		break;
				case 'o': getroot (str, outroot);	*str = '\0'; break;
				case 'c': scale = atof (str);		*str = '\0'; break;
				case '@': control = *str++;		*str = '\0'; break;
				default: break;
			}
		} else switch (k) {
			case 0:	getroot (argv[i], imgroot);	k++; break;
			default: break;
		}
	}
	if (k < 1) {
		printf ("Usage:	%s <microPET_data>\n", program);
		printf ("e.g.:	%s m1042-cft1_v1\n", program);
		printf ("\toption\n");
		printf ("\t-x	flip x\n");
		printf ("\t-y	flip y\n");
		printf ("\t-z	flip z\n");
		printf ("\t-w	create frame duration listing for use with actmapf_4dfp -w\n");
		printf ("\t-c<flt>	scale all voxel values by specified factor\n");
		printf ("\t-o<str>	name 4dfp output using specified string (default same as input)\n");
		printf ("\t-@<b|l> output big or little endian (default input endian)\n");
		exit (1);
	}

	sprintf (hdrfile, "%s.img.hdr", imgroot);
	printf ("Reading: %s\n", hdrfile);
	if (!(fphdr = fopen (hdrfile, "r"))) errr (program, hdrfile);
	while (fgets (string, MAXL, fphdr)) {
		if (m = split (string, srgv, MAXL) != 2) continue;
		if (!strcmp (srgv[0], "modality"))		modality = atoi (srgv[1]);
		if (!strcmp (srgv[0], "calibration_factor"))	calibration_factor = atof (srgv[1]);
		if (!strcmp (srgv[0], "file_type"))		ifile_type = atoi (srgv[1]);
		if (!strcmp (srgv[0], "data_type"))		idata_type = atoi (srgv[1]);
		if (!strcmp (srgv[0], "data_order"))		idata_order = atoi (srgv[1]);
		if (!strcmp (srgv[0], "acquisition_mode"))	iacquisition_mode = atoi (srgv[1]);
		if (!strcmp (srgv[0], "total_frames"))		imgdim[3] = atoi (srgv[1]);
		if (!strcmp (srgv[0], "x_dimension"))		imgdim[0] = atoi (srgv[1]);
		if (!strcmp (srgv[0], "y_dimension"))		imgdim[1] = atoi (srgv[1]);
		if (!strcmp (srgv[0], "z_dimension"))		imgdim[2] = atoi (srgv[1]);
		if (!strcmp (srgv[0], "pixel_size"))		voxdim[0] = voxdim[1] = atof (srgv[1])*10;	/* PET */
		if (!strcmp (srgv[0], "axial_crystal_pitch"))	voxdim[2] = atof (srgv[1])*5;			/* PET */
		if (!strcmp (srgv[0], "pixel_size_x"))		voxdim[0] = atof (srgv[1]);	/* CT */
		if (!strcmp (srgv[0], "pixel_size_y"))		voxdim[1] = atof (srgv[1]);	/* CT */
		if (!strcmp (srgv[0], "pixel_size_z"))		voxdim[2] = atof (srgv[1]);	/* CT */
	}
	printf ("image dimensions    %10d%10d%10d%10d\n",   imgdim[0], imgdim[1], imgdim[2], imgdim[3]);
	printf ("voxel dimensions    %10.4f%10.4f%10.4f\n", voxdim[0], voxdim[1], voxdim[2]);
	printf ("file type %s\n", file_type[ifile_type]);
	if (ifile_type != 5 && ifile_type != 8) {
		printf ("%s: %s file type conversion not supported\n", program, imgroot);
		exit (-1);
	}
	printf ("data type %s\n", data_type[idata_type]);
	if (idata_type < 0 || idata_type > 7) {
		printf ("%s: %s has illegal data type (%d)\n", program, imgroot, idata_type);
		exit (-1);
	}

	if (!(framedat = (FRAMEDAT *) malloc (imgdim[3] * sizeof (FRAMEDAT)))) errm (program);
	rewind (fphdr);
	while (fgets (string, MAXL, fphdr)) {
		if (m = split (string, srgv, MAXL) < 2) continue;
		if (!strcmp (srgv[0], "frame"))	{
			k = atoi (srgv[1]);
			if (k < 0 || k >= imgdim[3]) errf (hdrfile);
		}
		if (!strcmp (srgv[0], "event_type"))		framedat[k].event_type = atoi (srgv[1]);
		if (!strcmp (srgv[0], "gate"))			framedat[k].gate = atoi (srgv[1]);
		if (!strcmp (srgv[0], "bed"))			framedat[k].bed = atoi (srgv[1]);
		if (!strcmp (srgv[0], "bed_offset"))		framedat[k].bed_offset = atof (srgv[1]);
		if (!strcmp (srgv[0], "ending_bed_offset"))	framedat[k].ending_bed_offset = atof (srgv[1]);
		if (!strcmp (srgv[0], "vertical_bed_offset"))	framedat[k].vertical_bed_offset = atof (srgv[1]);
		if (!strcmp (srgv[0], "data_file_pointer")) {
				framedat[k].data_file_pointer[0] = atoi (srgv[1]);
				framedat[k].data_file_pointer[1] = atoi (srgv[2]);
		}
		if (!strcmp (srgv[0], "frame_start"))		framedat[k].frame_start = atof (srgv[1]);
		if (!strcmp (srgv[0], "frame_duration"))	framedat[k].frame_duration = atof (srgv[1]);
		if (!strcmp (srgv[0], "scale_factor"))		framedat[k].scale_factor = atof (srgv[1]);
		if (!strcmp (srgv[0], "minimum"))		framedat[k].minimum = atof (srgv[1]);
		if (!strcmp (srgv[0], "maximum"))		framedat[k].maximum = atof (srgv[1]);
	}
	write_framing (stdout);
	vdim = imgdim[0]*imgdim[1]*imgdim[2];
	if (!(imgf = (float *) malloc (vdim * sizeof (float)))) errm (program);
	if (!(imgs = (short *) malloc (vdim * sizeof (short)))) errm (program);
	if (!(imgd = (int *)   malloc (vdim * sizeof (int))))   errm (program);

	printf ("data_type=%d %s\n", idata_type, data_type[idata_type]);
	switch (idata_type) {
		case 2:	bytepix = 2;	isbig = 0;	iptr = (char *) imgs;	break;
		case 3:	bytepix = 4;	isbig = 0;	iptr = (char *) imgd;	break;
		case 4:	bytepix = 4;	isbig = 0;	iptr = (char *) imgf;	break;
		case 5:	bytepix = 4;	isbig = 1;	iptr = (char *) imgf;	break;
		case 6:	bytepix = 2;	isbig = 1;	iptr = (char *) imgs;	break;
		case 7:	bytepix = 4;	isbig = 1;	iptr = (char *) imgd;	break;
		default: fprintf (stderr, "%s: %s data type not supported\n", program, imgroot);
			exit (-1);
			break;
	}
	CPUisbig = (CPU_is_bigendian()) ? 1 : 0;
	swab_flag = CPUisbig ^ isbig;
	if (!control) control = (isbig) ? 'b' : 'l';
	printf ("CPUisbig=%d, swab_flag=%d\n", CPUisbig, swab_flag);

/**********************/
/* process all frames */
/**********************/
	q = scale;
	if (calibration_factor > 0.) q *= calibration_factor;
	sprintf (imgfile, "%s.img", imgroot);
	printf ("Reading: %s\n", imgfile);
	if (!(fpimg = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!strlen (outroot)) strcpy (outroot, imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);
	printf ("Writing: %s\n", outfile);
	if (!(fpout = fopen (outfile, "wb"))) errw (program, outfile);
	printf ("frame");
	if (!debug) for (k = 0; k < imgdim[3]; k++) {printf (" %d", k + 1); fflush (stdout);
		if (fread (iptr, bytepix, vdim, fpimg) != vdim) errr (program, imgfile);
		for (i = 0; i < vdim; i++) {
			switch (idata_type) {
			case 2:	if (swab_flag) swab2 ((char *) &imgs[i]); imgf[i] =     q*framedat[k].scale_factor*imgs[i]; break;
			case 3:	if (swab_flag) swab4 ((char *) &imgd[i]); imgf[i] =     q*framedat[k].scale_factor*imgd[i]; break;
			case 4:	if (swab_flag) swab4 ((char *) &imgf[i]); imgf[i] = scale*framedat[k].scale_factor*imgf[i]; break;
			case 5:	if (swab_flag) swab4 ((char *) &imgf[i]); imgf[i] = scale*framedat[k].scale_factor*imgf[i]; break;
			case 6:	if (swab_flag) swab2 ((char *) &imgs[i]); imgf[i] =     q*framedat[k].scale_factor*imgs[i]; break;
			case 7:	if (swab_flag) swab4 ((char *) &imgd[i]); imgf[i] =     q*framedat[k].scale_factor*imgd[i]; break;
			}
			if (imgf[i] > glmax) glmax = imgf[i];
			if (imgf[i] < glmin) glmin = imgf[i];
		}
		if (xflag) flipx (imgf, imgdim + 0, imgdim + 1, imgdim + 2);
		if (yflag) flipy (imgf, imgdim + 0, imgdim + 1, imgdim + 2);
		if (zflag) flipz (imgf, imgdim + 0, imgdim + 1, imgdim + 2);
		if (gwrite ((char *) imgf, sizeof (float), vdim, fpout, control)) errw (program, outfile);
	} printf ("\n"); fflush (stdout);
	if (fclose (fpimg)) errr (program, imgfile);
	if (fclose (fpout)) errw (program, outfile);
	
/***************/
/* ifh and hdr */
/***************/
	if (writeifhe (program, outroot, imgdim, voxdim, 2, control)) errw (program, outroot);
	sprintf (command, "ifh2hdr %s -r%dto%d", outroot, (int) (glmin - 0.5), (int) (glmax + 0.5));
	printf ("%s\n", command);
	status |= system (command);

/*******/
/* rec */
/*******/
	startrec (outfile, argc, argv, rcsid);
	sprintf (command, "Voxel values scaled by %.2f\n", scale); printrec (command);
	rewind (fphdr);	while (fgets (string, MAXL, fphdr)) {
		if (m = split (string, srgv, MAXL) < 2) continue;
		k = atoi (srgv[1]); if (k < 0) continue;
		if (0) {	/* initiate else if chain */
		} else if (!strcmp (srgv[0], "file_type") && k < 11) {
			sprintf (command, "%s: %s\n", srgv[0], file_type[k]); printrec (command);
		} else if (!strcmp (srgv[0], "acquisition_mode") && k < 9) {
			sprintf (command, "%s: %s\n", srgv[0], acquisition_mode[k]); printrec (command);
		} else if (!strcmp (srgv[0], "bed_motion") && k < 3) {
			sprintf (command, "%s: %s\n", srgv[0], bed_motion[k]); printrec (command);
		} else if (!strcmp (srgv[0], "tx_src_type") && k < 3) {
			sprintf (command, "%s: %s\n", srgv[0], tx_src_type[k]); printrec (command);
		} else if (!strcmp (srgv[0], "data_type") && k < 8) {
			sprintf (command, "%s: %s\n", srgv[0], data_type[k]); printrec (command);
		} else if (!strcmp (srgv[0] + 1, "_filter") && k < 8) {
			sprintf (command, "%s: %s cutoff %s/cm\n", srgv[0], xyz_filter[k], srgv[2]); printrec (command);
		} else if (!strcmp (srgv[0], "rebinning_type") && k < 4) {
			sprintf (command, "%s: %s\n", srgv[0], rebinning_type[k]); printrec (command);
		} else if (!strcmp (srgv[0], "recon_algorithm") && k < 7) {
			sprintf (command, "%s: %s\n", srgv[0], recon_algorithm[k]); printrec (command);
		} else if (!strcmp (srgv[0], "deadtime_correction_applied") && k < 3) {
			sprintf (command, "%s: %s\n", srgv[0], deadtime_correction_applied[k]); printrec (command);
		} else if (!strcmp (srgv[0], "normalization_applied") && k < 5) {
			sprintf (command, "%s: %s\n", srgv[0], normalization_applied[k]); printrec (command);
		} else if (!strcmp (srgv[0], "attenuation_applied") && k < 8) {
			sprintf (command, "%s: %s\n", srgv[0], attenuation_applied[k]); printrec (command);
		} else if (!strcmp (srgv[0], "scatter_correction") && k < 6) {
			sprintf (command, "%s: %s\n", srgv[0], scatter_correction[k]); printrec (command);
		} else if (!strcmp (srgv[0], "calibration_units") && k < 3) {
			sprintf (command, "%s: %s\n", srgv[0], calibration_units[k]); printrec (command);
		} else if (!strcmp (srgv[0], "dose_units") && k < 3) {
			sprintf (command, "%s: %s\n", srgv[0], dose_units[k]); printrec (command);
		} else if (!strcmp (srgv[0], "subject_orientation") && k < 9) {
			sprintf (command, "%s: %s\n", srgv[0], subject_orientation[k]); printrec (command);
		} else if (!strcmp (srgv[0], "subject_length_units") && k < 4) {
			sprintf (command, "%s: %s\n", srgv[0], subject_length_units[k]); printrec (command);
		} else if (!strcmp (srgv[0], "subject_weight_units") && k < 5) {
			sprintf (command, "%s: %s\n", srgv[0], subject_weight_units[k]); printrec (command);
		}
	}
	k = 0; while (strcmp (keywords[k], "end_of_header")) {
		rewind (fphdr);
		while (fgets (string, MAXL, fphdr)) {
			strcpy (command, string);
			if (m = split (string, srgv, MAXL) < 2) continue;
			if (!strcmp (srgv[0], keywords[k])) printrec (command);
		}
		k++;
	}
	if (modality == 0) { /* PET */
		sprintf (recfile, "%s.rec", outfile);
		if (!(fpout = fopen (recfile, "a"))) errw (program, outfile);
		write_framing (fpout);
		fclose (fpout);
	}
	endrec ();

	fclose (fphdr);

	if (weifile) {
/**************************************/
/* create fraem duration listing file */
/**************************************/
		sprintf (outfile, "%s.framedur", outroot);
		printf ("Writing: %s\n", outfile);
		if (!(fpout = fopen (outfile, "w"))) errw (program, outfile);
		for (i = 0; i < imgdim[3]; i++) fprintf (fpout, "%10.2f\n", framedat[i].frame_duration);
		if (fclose (fpout)) errw (program, outfile);
	}

/*****************/
/* free and exit */
/*****************/
	free (imgf);
	free (imgs);
	free (imgd);
	free (framedat);
	exit (status);
}


