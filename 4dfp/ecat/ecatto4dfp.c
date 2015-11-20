/*$Header: /data/petsun4/data1/src_solaris/ecat/RCS/ecatto4dfp.c,v 1.6 2011/08/31 16:24:35 larsc Exp $*/
/*$Log: ecatto4dfp.c,v $
 * Revision 1.6  2011/08/31  16:24:35  larsc
 * Changed Start column in rec file to Midpoint and added Ecat_Frame. Eliminated frame factor stout.
 *
 * Revision 1.5  2011/07/08  17:21:01  larsc
 * Made frame number in rec file refer to 4dfp, not ecat.
 *
 * Revision 1.4  2011/07/07  16:01:37  larsc
 * Added frame sorting.
 *
 * Revision 1.3  2011/01/24  21:50:05  larsc
 * Tweaked rec file info
 *
 * Revision 1.2  2011/01/12  19:51:54  larsc
 * Added header information to recfile
 *
 * Revision 1.1  2010/09/01  17:38:26  larsc
 * Initial revision
 **/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "load_volume.h"

#include "ifh.h"
#include "endianio.h"
#include "Getifh.h"
#include "rec.h"

#define MAXL		256
#define TRANSVERSE	2
#define CUBIC		0

static char	program[MAXL];
static char	rcsid[] = "$Id: ecatto4dfp.c,v 1.6 2011/08/31 16:24:35 larsc Exp $";

void setprog (char *program, char **argv)
{
	char *ptr;
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0];
	else ptr++;
	strcpy (program, ptr);
}

int frame_cmp (const void *a, const void *b)
{
	const MatDirNode	*ia = *((const MatDirNode **) a);
	const MatDirNode	*ib = *((const MatDirNode **) b);
	struct Matval		matvala, matvalb;
	mat_numdoc (ia->matnum, &matvala);
	mat_numdoc (ib->matnum, &matvalb);
	return matvala.frame - matvalb.frame;
}

int main (int argc, char *argv[])
{
	char		ecatspec[MAXL], ecatname[MAXL];
	char		imgroot[MAXL], imgname[MAXL];
	FILE		*fptr;
	struct Matval	matval;
	MatrixFile	*mptr;
	MatrixData	*matrix;
	MatDirNode	*entry, **entries;
	Image_subheader	*ih;

	float		global_max, global_min;
	int		volume, plane, line, nframes, iframe, jframe, matnum, *missing, nmissing;
	float		scale_factor, calibration_factor/*, scan_time*/;
	int		cubic = CUBIC;

	float		*img;
	short		*ecat;
	
	int		i, j, k;
	char		string[MAXL], *ptr, c;
	char		control = '\0';
	int		isbig;

	IFH		ifh;
	
	setprog (program, argv);

	for (j = 0, i = 1; i < argc; i++)
	if (*argv[i] == '-')
	{
		strncpy (string, argv[i], MAXL); ptr = string;
		while (c = *ptr++) switch (c)
		{
			case '@': control = *ptr++; *ptr = '\0'; break;
		}
	}
	else switch (j)
	{
		case 0: strncpy (ecatspec, argv[i], MAXL);	j++; break;
		case 1: getroot (argv[i], imgroot);		j++; break;
	}
	if (j != 2)
	{
		fprintf (stderr, "Usage:\t%s <ecat> <(4dfp)image>\n", program);
		fprintf (stderr, "\t-@<b|l>\toutput big or little endian (default input endian)\n");
		exit (1);
	}

/******** Get filenames and open ********/
	matspec (ecatspec, ecatname, &matnum);	/* matnum is unused */
	if (!(mptr = matrix_open (ecatname, MAT_READ_ONLY, MAT_UNKNOWN_FTYPE)))
	{
		fprintf (stderr, "%s: cannot open %s as an ECAT image file\n", program, ecatname);
		exit (-1);
	}
	if (mptr->mhptr->file_type != InterfileImage
	 && mptr->mhptr->file_type != PetImage
	 && mptr->mhptr->file_type != PetVolume
	 && mptr->mhptr->file_type != ByteImage
	 && mptr->mhptr->file_type != ByteVolume)
	{
		fprintf (stderr, "%s: filetype not supported\n", program);
		exit (-1);
	}
	sprintf (imgname, "%s.4dfp.img", imgroot);
	if (!(fptr = fopen (imgname, "w"))) errw (program, imgname);
	
	if (!control) control = (CPU_is_bigendian ()) ? 'b' : 'l';
	
	calibration_factor = mptr->mhptr->calibration_factor;
	fprintf (stdout, "Calibration Factor := %10e\n", calibration_factor);

	if (!(missing = (int *)         calloc (mptr->mhptr->num_frames,  sizeof (int)))
	 || !(entries = (MatDirNode **) malloc (mptr->mhptr->num_frames * sizeof (MatDirNode *)))) errm (program);

	sprintf (string, "%s.4dfp.img.rec", imgroot);
	startrece (string, argc, argv, rcsid, control);
	printrec ("Frame     \t    Length\t  Midpoint\t     Start\t Frame_Min\t Frame_Max\t Decay_Fac\tEcat_Frame\n");

	for (nframes = 0, entry = mptr->dirlist->first; entry; entry = entry->next) entries[nframes++] = entry;
	qsort (entries, nframes, sizeof (MatDirNode *), frame_cmp);

	for (iframe = 1, jframe = 0, nmissing = 0/*, scan_time = 0*/; jframe < nframes; iframe++, jframe++)
	{
		entry = entries[jframe];
		mat_numdoc (entry->matnum, &matval);
		while (iframe < matval.frame) {
			fprintf (stderr, "%s: %s frame %d not included\n", program, ecatname, iframe);
			missing[nmissing++] = iframe++;
		}
		if (!(matrix = load_volume (mptr, matval.frame, cubic)))
		{
			fprintf (stderr, "%s: ecat frame %d not found\n", program, matval.frame);
			exit (-1);
		}
		if ( matrix->data_type != SunShort && matrix->data_type != VAX_Ix2)
		{
			fprintf (stderr, "%s: only integer 2 images are currently supported\n", program);
			exit(-1);
		}

		scale_factor = matrix->scale_factor * calibration_factor;
/*		fprintf (stdout, "Scale Factor := %10e\tTotal Factor := %10e\n", matrix->scale_factor, scale_factor);
*/
		if (jframe == 0)
		{
			ifh.scaling_factor[0] = ifh.scaling_factor[1] = 10 * matrix->pixel_size;
			ifh.scaling_factor[2] = 10 * matrix->z_size;
			ifh.matrix_size[0] = matrix->xdim;
			ifh.matrix_size[1] = matrix->ydim;
			ifh.matrix_size[2] = matrix->zdim;
			line   = ifh.matrix_size[0];
			plane  = ifh.matrix_size[1] * line;
			volume = ifh.matrix_size[2] * plane;
			if (!(img = (float *) malloc (volume * sizeof (float)))) errm (program);
			global_min = matrix->data_min;
			global_max = matrix->data_max;
		}
		else
		{
			if (matrix->data_min < global_min) global_min = matrix->data_min;
			if (matrix->data_max > global_max) global_max = matrix->data_max;
		}

/******** Flip Z, and write to output file ********/
		ecat = (short *) matrix->data_ptr;
		for (i = volume-plane; i >= 0;     i -= plane)
/*		for (j = 0;            j <  plane; j += line)
		for (k = 0;            k <  line;  k++, ecat++) img[i+j+k] = (*ecat) * scale_factor;
*/		for (j = 0;            j <  plane; j++, ecat++) img[i+j] = (*ecat) * scale_factor;
		if (ewrite (img, volume, control, fptr)) errw (program, imgname);

		ih = (Image_subheader*) matrix->shptr;
		sprintf (string, "Frame_%-4d\t%10d\t%10.2f\t%10d\t%10d\t%10d\t%10.5f\t%10d\n", jframe+1,
			ih->frame_duration,
			(ih->frame_start_time + ih->frame_duration / 2.0) / 1000/*scan_time += ih->frame_duration / 1000*/,
			ih->frame_start_time,
			ih->image_min,
			ih->image_max,
			ih->decay_corr_fctr,
			iframe);
		printrec (string);
	}
	
	sprintf (string, "%s Missing Frames:", ecatname);
	printrec (string);
	for (i = 0; i < nmissing; i++)
	{
		sprintf (string, " %d", missing[i]);
		printrec (string);
	}
	printrec ("\n");
	endrec ();

/******** Write ifh, hdr, and rec files ********/
	ifh.matrix_size[3]    = nframes;
	ifh.orientation       = TRANSVERSE;
	sprintf (string, "%s.4dfp.ifh", imgroot);
	if (writeifhe (program, string, ifh.matrix_size, ifh.scaling_factor, ifh.orientation, control)) errw (program, string);

	sprintf (string,  "ifh2hdr %s -r%fto%f",  imgroot, global_min, global_max);
	system (string);
	
/******** Free, close and quit ********/
	matrix_close (mptr);
	fclose (fptr);
	free (img);
	exit (0);
}
