/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/4dfp-format.c,v 1.15 2013/01/11 06:37:47 avi Exp $*/
/*$Log: 4dfp-format.c,v $
 * Revision 1.15  2013/01/11  06:37:47  avi
 * assume signed float 4 byte 4dfp input in parse_4dfp()
 *
 * Revision 1.14  2012/02/09  00:28:40  coalsont
 * changed write() to fwrite() to get around 2GB limitation
 *
 * Revision 1.13  2011/12/02  05:09:26  coalsont
 * changed file writing to not use mmap
 *
 * Revision 1.12  2011/09/13  03:42:42  avi
 * add argument to control saving of mmppix and center fields
 *
 * Revision 1.11  2010/08/28  04:40:50  coalsont
 * fixed revorder to cause matrix size fields to write correctly
 *
 * Revision 1.10  2010/08/16  18:18:44  coalsont
 * used 64bit integers to handle file sizes
 *
 * Revision 1.9  2010/07/08  15:27:21  coalsont
 * used getroot to strip extensions
 *
 * Revision 1.8  2009/08/20  23:15:28  coalsont
 * voxdim is always positive
 *
 * Revision 1.7  2009/08/17  22:57:44  coalsont
 * ground truth from S2T_4dfp
 *
 * Revision 1.6  2009/08/07  22:11:43  coalsont
 * check file size to return null instead of crashing
 *
 * Revision 1.4  2009/07/29  20:02:28  coalsont
 * minor typos
 *
 * Revision 1.3  2009/07/24  18:53:51  coalsont
 * fix for odd axis order nift, was flipping wrong axes to preserve conventions
 *
 * Revision 1.2  2009/07/14  21:41:09  coalsont
 * fix: use root, not name of data file
 *
 * Revision 1.1  2009/07/08  01:01:35  coalsont
 * Initial revision
 **/
static char* rcsid = "$Id: 4dfp-format.c,v 1.15 2013/01/11 06:37:47 avi Exp $";
/* Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009 */
/* Kevin P. Barry [ta0kira@users.berlios.de], BMCLAB, 22 May 2009 */

#include "4dfp-format.h"

#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <libgen.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include <ifh.h>
#include <Getifh.h>
#include "t4_io.h"
#include <rec.h>
#include "endianio.h"

#include "common-format.h"
#include "transform.h"
#include "nifti-format.h"
#include "4dfp-format.h"

#define HEADER_EXT ".4dfp.ifh"
#define IMAGE_EXT ".4dfp.img"
/* used for output only */
#define T4_EXT ".4dfp.img_to_atlas_t4"

static void get_file_names(char *fFile, char **hHeader, char **iImage, char **tT4)
{
	size_t original_length = strlen(fFile);
	/*int extension_index = original_length;
	if (original_length > strlen(HEADER_EXT) &&
		strcmp(HEADER_EXT, fFile + original_length - strlen(HEADER_EXT)) == 0)
		extension_index = original_length - strlen(HEADER_EXT);//*/
	
	char* temp = calloc(1, original_length + 1);
	getroot(fFile, temp);
	int extension_index = strlen(temp);
	free(temp);//use getroot only for length, the rest of the code already works as is
	
	*hHeader = calloc(1, extension_index + strlen(HEADER_EXT) + 1);
	strncpy(*hHeader, fFile, extension_index);
	strncpy(*hHeader + extension_index, HEADER_EXT, strlen(HEADER_EXT));
	
	*iImage = calloc(1, extension_index + strlen(IMAGE_EXT) + 1);
	strncpy(*iImage, fFile, extension_index);
	strncpy(*iImage + extension_index, IMAGE_EXT, strlen(IMAGE_EXT));

	*tT4 = calloc(1, extension_index + strlen(T4_EXT) + 1);
	strncpy(*tT4, fFile, extension_index);
	strncpy(*tT4 + extension_index, T4_EXT, strlen(T4_EXT));
	
}


static void free_file_names(char *hHeader, char *iImage, char *tT4)
{
	free(hHeader);
	free(iImage);
	free(tT4);
}


common_format *parse_4dfp(char *fFile, char* t4name)
{
	char *header_name = NULL, *image_name = NULL, *junk = NULL;
	get_file_names(fFile, &header_name, &image_name, &junk);
	int i, j, warn_flag;
	float tempf;
	IFH header_data;

	if (Getifh(header_name, &header_data) != 0)
	{
		fprintf(stderr, "error opening header '%s'\n", header_name);
		free_file_names(header_name, image_name, junk);
		return NULL;
	}
	
	common_format *file_content = calloc(1, sizeof(common_format));
	if (!file_content)
	{
		fprintf(stderr, "allocation error: %s\n", strerror(errno));
		free_file_names(header_name, image_name, junk);
		return NULL;
	}
	
	vindex axes[4] = AXES_4DFP;
	memcpy(file_content->order, axes, 4 * sizeof(vindex));
	
	file_content->f_offset = 0;//file position tracking for eread
	
/***********************************************/
/* signed 4 byte float is assumed in 4dfp data */
/***********************************************/
	file_content->isFloat = 1;
	file_content->isSigned = 1;
	file_content->bperpix = 4;
/**********************************************************************/
/* check ifh for non-standard values and report as warnings to stderr */
/**********************************************************************/
	warn_flag = 0;
	if (!strstr(header_data.number_format, "float")) warn_flag++;
	if (header_data.number_of_bytes_per_pixel != 4) warn_flag++;
	if (warn_flag) {
		fprintf (stderr, "Warning: nonstandard 4dfp conventions in %s\n", header_name);
	}
/*
	if (file_content->isFloat && (file_content->bperpix != 4 && file_content->bperpix != 8 && file_content->bperpix != 16))
	{
		fprintf(stderr, "unsupported float width: %i\n", file_content->bperpix);
		free_file_names(header_name, image_name, junk);
	}
	if (!file_content->isFloat && (file_content->bperpix != 1 && file_content->bperpix != 2 && file_content->bperpix != 4 && file_content->bperpix != 8))
	{
		fprintf(stderr, "unsupported integer width: %i\n", file_content->bperpix);
		free_file_names(header_name, image_name, junk);
	}
*/
	for (i = 0; i < 3; ++i)
		for (j = 0; j < 4; ++j)
			file_content->sform[i][j] = 0.0;
	
	file_content->haveCenter = 1;
	file_content->scale = 1.0;
	file_content->offset = 0.0;
	/* find axis order first in case we need it for the sform due to nonexistant t4 */
	
	vindex tempv;/* for axes swaps */
	switch (header_data.orientation)
	{/* do only virtual axes swaps for now, sform row order is changed below */
		case 4: /* swaps anatomical x and z first, second swap makes it the correct order */
			tempv = file_content->order[1];
			file_content->order[1] = file_content->order[0];
			file_content->order[0] = tempv;
		case 3:
			tempv = file_content->order[2];
			file_content->order[2] = file_content->order[1];
			file_content->order[1] = tempv;
		case 2: /* nothing to be done, but to identify bad headers */
			break;
			
		default:
			fprintf(stderr, "%s image orientation not recognized\n", header_name);
			free_file_names(header_name, image_name, junk);
			return NULL;
	}
	
	int t4file;
	float t4trans[4][4], scale, *antiwarning;
	antiwarning = (float*)t4trans;
	double temp, temp2, det = 0.0;
	if (t4name)
	{
		t4file = t4_read(t4name, antiwarning);
		if (t4file)
		{
			fprintf(stderr, "error opening t4 file '%s': %s\n", t4name,
					strerror(errno));
			free_file_names(header_name, image_name, junk);
			return NULL;
		}
		iget_t4_scale(t4name, &scale);
		file_content->scale = scale;
		for (i = 0; i < 3; ++i)
		{/* transpose so the columns are the first index */
			for (j = i + 1; j < 4; ++j)
			{
				tempf = t4trans[i][j];
				t4trans[i][j] = t4trans[j][i];
				t4trans[j][i] = tempf;
			}
		}
		for (i = 0; i < 3; ++i)
		{
			printf("t4: %f\t%f\t%f\t%f\n", t4trans[i][0], t4trans[i][1], t4trans[i][2], t4trans[i][3]);
		}
		file_content->offset = 0.0;/* invert the t4 to get sform */
		for (i = 0; i < 3; ++i)
		{/* determinant, can use diagonals trick for 3x3 */
			temp = 1.0;
			temp2 = 1.0;
			for (j = 0; j < 3; ++j)
			{
				temp *= t4trans[j][(i + j) % 3];
				temp2 *= t4trans[j][(i - j + 3) % 3];
			}
			det += temp - temp2;
		}/* adjugate */
		int a, b, c, d;
		for (i = 0; i < 3; ++i)
		{
			a = (i + 1) % 3;/* wraparound avoids +/- alternation */
			b = (i + 2) % 3;
			for (j = 0; j < 3; ++j)
			{
				c = (j + 1) % 3;
				d = (j + 2) % 3;
				file_content->sform[j][i] = t4trans[a][c] * t4trans[b][d] - t4trans[a][d] * t4trans[b][c];
			}
		}/* divide by determinant */
		for (i = 0; i < 3; ++i)
		{
			for (j = 0; j < 3; ++j)
			{
				file_content->sform[i][j] /= det;
			}
		}/* multiply offset by inverse t4 and negate to start building the nifti sform */
		for (i = 0; i < 3; ++i)
		{
			for (j = 0; j < 3; ++j)
			{
				file_content->sform[i][3] -= file_content->sform[i][j] * t4trans[j][3];
			}
		}
	} else {/* take coords as they are from mmppix and center, adjusting for the right orientation */
		for (i = 0; i < 3; ++i)
			file_content->sform[file_content->order[i]][i] = 1.0;
	}
	
	file_content->file = header_name;
	header_name = NULL;
	
	file_content->dims = header_data.number_of_dimensions;
	
	file_content->length[0] = header_data.matrix_size[0];
	file_content->length[1] = header_data.matrix_size[1];
	file_content->length[2] = header_data.matrix_size[2];
	file_content->length[3] = header_data.matrix_size[3];
	
	double spacing[3], center[3];/* adjust spacing/center to something sane before adding them to the sform */
	spacing[0] = header_data.mmppix[0];
	spacing[1] = header_data.mmppix[1];
	spacing[2] = header_data.mmppix[2];
	file_content->timestep = (header_data.scaling_factor[3] ? header_data.scaling_factor[3] : 1.0);
	
	center[0] = header_data.center[0];
	center[1] = header_data.center[1];
	center[2] = header_data.center[2];
	
	file_content->endian = strcmp(header_data.imagedata_byte_order, "littleendian") != 0;
	
	for (i = 0; i < 3; ++i)
	{
		printf("input:\t%f\t%f\n", center[i], spacing[i]);
	}
	
	switch (header_data.orientation)
	{
		case 2:
		case 4://unlike vrtflip, ground truth from S2T says x and y flip (i and k), happens to match i and k flip for ori 2
			center[2] = -center[2];
			spacing[2] = -spacing[2];
			center[2] = spacing[2] * (file_content->length[2] + 1) - center[2];
		case 3:
			center[0] = -center[0];
			spacing[0] = -spacing[0];
			center[0] = spacing[0] * (file_content->length[0] + 1) - center[0];
			break;
		/* no default needed, unknown caught above */
	}
	
	file_content->size = 1;
	if (file_content->length[ file_content->order[0] ]) file_content->size *= file_content->length[ file_content->order[0] ];
	if (file_content->length[ file_content->order[1] ]) file_content->size *= file_content->length[ file_content->order[1] ];
	if (file_content->length[ file_content->order[2] ]) file_content->size *= file_content->length[ file_content->order[2] ];
	if (file_content->length[ file_content->order[3] ]) file_content->size *= file_content->length[ file_content->order[3] ];
	
	/* quirk: center and spacing have same sign, so the math is strange here */
	/* adjust for fortran 1-indexing */
	for (i = 0; i < 3; ++i)
	{
		center[i] -= spacing[i];
		/* adjust center to be the coords of index 0, 0, 0, like nifti notation */
		center[i] = -center[i];
	}	
	/* half voxel correction, corner to center */
	/*for (i = 0; i < 3; ++i)
	{
		*//* FIX: via t4imgs_4dfp/imgvalx.f, the interpolation method uses coordinates specified from first spatial sample */
		/* therefore, for correct interpolation, origin is treated as voxel center despite stated corner of voxel convention (which really isn't used for anything) */
		/*center[i] += spacing[i] / 2.0;
	}*/
	/* spacing and center are now sane, add them to the sform */
	/* can't be done in place or information would be lost */
	/* first copy it, axis order should already be applied */
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			t4trans[i][j] = file_content->sform[i][j];
		}
	}	
	for (i = 0; i < 3; ++i)
	{/* apply spacing */
		for (j = 0; j < 3; ++j)
		{
			file_content->sform[i][j] = t4trans[i][j] * spacing[j];
		}/* center is already copied, add ifh center */
		for (j = 0; j < 3; ++j)
		{
			file_content->sform[i][3] += center[j] * t4trans[i][j];
		}
	}
	
	//image_name = header_data.name_of_data_file;
	//name_of_data_file must be IGNORED due to convention
	int image_file = open(image_name, O_RDONLY);
	if (image_file < 0)
	{
		fprintf(stderr, "error opening image '%s': %s\n", image_name,
				strerror(errno));
		free_file_names(header_name, image_name, junk);
		return NULL;
	}
	
	struct stat file_stats;
	if (fstat(image_file, &file_stats) != 0)
	{
		fprintf(stderr, "error reading image '%s': %s\n", image_name,
				strerror(errno));
		free_file_names(header_name, image_name, junk);
		close(image_file);
		return NULL;
	}
	
	int64_t file_size = file_stats.st_size;
	int64_t rounded_size = file_size;
	int64_t page_size = sysconf(_SC_PAGESIZE);
	if (rounded_size % page_size) rounded_size += page_size + (-rounded_size % page_size);
	
	if (file_size < file_content->size * file_content->bperpix)
	{
		fprintf(stderr, "file length is too short for specified dimensions\n");
		return NULL;
	}
	
	file_content->voxels = mmap(NULL, rounded_size, PROT_READ, MAP_SHARED, image_file, 0x00);
	if (!file_content->voxels)
	{
		fprintf(stderr, "error reading image '%s': %s\n", image_name,
				strerror(errno));
		free_common_format(file_content);
		close(image_file);
		return NULL;
	}
	file_content->is_mmaped = 1;//signifies not malloced
	file_content->mmap_offset = 0;
	
	for (i = 0; i < 3; ++i)
	{
		printf("common:\t%f\t%f\t%f\t%f\n", file_content->sform[i][0], file_content->sform[i][1], file_content->sform[i][2], file_content->sform[i][3]);
	}
	
	/* if they use a t4 file, redefine order by dominant direction of the indexing vectors */
	/* since the header is now sane, we can do this with no problems */
	/* if the t4 file isnt used, anything not on the major axis is 0, so no problems */
	free(image_name);
	//close(image_file);
	return file_content;
}

int save_4dfp(common_format *iImage, char *fFile, int argc, char* argv[], char control, int save_center)
{
	char *header_name = NULL, *image_name = NULL, *t4_name = NULL;
	int i, j;
	float tempf;
	
	char reswap = ((local_endian() && (control == 'l' || control == 'L')) || (!local_endian() && (control == 'b' || control == 'B')));
	
	get_file_names(fFile, &header_name, &image_name, &t4_name);
	
	to_lpi(iImage);//doesnt actually change any spacing or center, just finds what needs to be flipped/swapped and stores in orientation/order
	
	/* reverse lookup to reorder header fields for axis swapping */
	vindex revorder[4]; /* used for consistency, needed in nifti writes due to header written before voxels */
	for (i = 0; i < 4; ++i)
	{
		revorder[iImage->order[i]] = i;
	}
	
	iImage->orientation ^= (uint8_t) 01 << revorder[0];
	iImage->orientation ^= (uint8_t) 01 << revorder[1];
	/* satisfy 4dfp 1, -1, -1 spacing convention so slice number doesnt jump in some viewers */
	/* NOTE: orientation was already set up to change to 1, 1, 1 spacing, and quirks flip the x and z spacing again later, so x and y need flipping now */
	
	auto_orient_header(iImage);/* do any needed voxel flipping on the common header too, before adding quirks */

	float spacing[4] = {0.0f};/* because writeifhe needs float, so lose some precision */
	spacing[3] = fabs(iImage->timestep);
	double center[3] = {0.0};
	/* mmppix is magnitude of the spacial vectors that i, j, and k map to */
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			spacing[i] += iImage->sform[j][i] * iImage->sform[j][i];
		}
		spacing[i] = sqrt(spacing[i]);
	}
	spacing[0] = -spacing[0];/* keep the +, -, - convention in the .ifh, x and z are flipped later */
	spacing[1] = -spacing[1];/* we do this here to specify we want the flips to take place before the t4 transform is applied */
	
	float t4trans[4][4] = {{0.0f}}, *antiwarning;
	antiwarning = (float*)t4trans;
	double temp, temp2, det = 0.0;
	t4trans[3][3] = 1.0;
	
	/* separate spacing/center from sform, put remainder into a t4 file if needed */
	/* first, invert sform to get t4, and the center for the plumb image */
	for (i = 0; i < 3; ++i)
	{/* rescale by the indexing vectors */
		for (j = 0; j < 3; ++j)
		{
			iImage->sform[i][j] /= spacing[j];
		}
	}
	for (i = 0; i < 3; ++i)
	{/* determinant, can use diagonals trick for 3x3 */
		temp = 1.0;
		temp2 = 1.0;
		for (j = 0; j < 3; ++j)
		{
			temp  *= iImage->sform[j][(i + j) % 3];
			temp2 *= iImage->sform[j][(i - j + 3) % 3];
		}
		det += temp - temp2;
		temp = 1.0;
	}/* adjugate */
	int a, b, c, d;
	for (i = 0; i < 3; ++i)
	{
		a = (i + 1) % 3;/* wraparound avoids +/- alternation */
		b = (i + 2) % 3;
		for (j = 0; j < 3; ++j)
		{
			c = (j + 1) % 3;
			d = (j + 2) % 3;
			t4trans[j][i] = iImage->sform[a][c] * iImage->sform[b][d] - iImage->sform[a][d] * iImage->sform[b][c];
		}
	}/* divide by determinant */
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			t4trans[i][j] /= det;
		}
	}/* we want t4 translation to be 0, so calculate the center for center/mmppix */
	for (i = 0; i < 3; ++i)
	{/* plug sform offset (true coords of 0, 0, 0) into t4 to get the coordinates 0, 0, 0, should be at after only applying spacing/center */
		for (j = 0; j < 3; ++j)
		{
			center[i] += iImage->sform[j][3] * t4trans[i][j];
		}
	}
			
	/* center to corner adjustment */
	/* FIX: see parse_4dfp, interpolation and non-interpolated display conventions cancel each other */
	/*for (i = 0; i < 3; ++i)
	{
		center[i] -= spacing[i] / 2.0;
	}*/

	/* redo fortran 1-indexing AFTER changing to reverse center quirk */
	for (i = 0; i < 3; ++i)
	{
		center[i] = -center[i];
		center[i] += spacing[i];
	}
	
	/* vrtflip_ can be substituted for these lines (use orientation = 2) */
	/* center in 4dfp is not relative to the first written voxel, orientation 2 flips X and Z */
	center[0]	=  spacing[0] * (iImage->length[revorder[0]] + 1) - center[0];
	center[0]	= -center[0];
	spacing[0]	= -spacing[0];
	center[2]	=  spacing[2] * (iImage->length[revorder[2]] + 1) - center[2];
	center[2]	= -center[2];
	spacing[2]	= -spacing[2];
	
	if (save_center) for (i = 0; i < 3; ++i)
	{
		printf("out:\t%f\t%f\n", center[i], spacing[i]);
	}
	
	IFH header_data = {
		/*.interfile = ,*/
		.version_of_keys = "3.3",
		/*.conversion_program = ,*/
		/*.name_of_data_file = ,*/
		.number_format = "float",
		/*.imagedata_byte_order = ,*/
		.number_of_bytes_per_pixel = 4,
		.number_of_dimensions = 4,
		.matrix_size = {
			iImage->length[revorder[0]],
			iImage->length[revorder[1]],
			iImage->length[revorder[2]],
			iImage->length[revorder[3]]
		},
		.orientation = 2,
		.scaling_factor = {
			fabs(spacing[0]),
			fabs(spacing[1]),
			fabs(spacing[2]),
			fabs(spacing[3])
		},
		.mmppix = {
			spacing[0],
			spacing[1],
			spacing[2]
		},
		.center = {
			center[0],
			center[1],
			center[2]
		}
	};
	
	
	/*snprintf(header_data.interfile, sizeof header_data.interfile,
			 "%s", iImage->file);*/
	for (i = 0; i < (int)(sizeof header_data.interfile - 1) && iImage->file[i]; ++i)
	{
		header_data.interfile[i] = iImage->file[i];
	}
	header_data.interfile[i] = '\0';
	
	/*snprintf(header_data.conversion_program, sizeof header_data.conversion_program,
			 "%s", argv[0]);*/

	for (i = 0; i < (int)(sizeof header_data.conversion_program - 1) && argv[0][i]; ++i)
	{
		header_data.conversion_program[i] = argv[0][i];
	}
	header_data.conversion_program[i] = '\0';
	
	/*snprintf(header_data.name_of_data_file, sizeof header_data.name_of_data_file,
			 "%s", basename(image_name));*/

	char* tempstr = basename(image_name);
	for (i = 0; i < (int)(sizeof header_data.name_of_data_file - 1) && tempstr[i]; ++i)
	{
		header_data.name_of_data_file[i] = tempstr[i];
	}
	header_data.name_of_data_file[i] = '\0';
	
	if (iImage->haveCenter && save_center)
	{
		if (Writeifh(argv[0], header_name, &header_data, control) != 0)
		{
			fprintf(stderr, "error writing header '%s'\n", header_name);
			free_file_names(header_name, image_name, t4_name);
			return 0x00;
		}
	} else {/* use non-center header form if nifti had no center */
		if (writeifhe(argv[0], header_name, header_data.matrix_size, header_data.scaling_factor, 2, control) != 0)
		{
			fprintf(stderr, "error writing header '%s'\n", header_name);
			free_file_names(header_name, image_name, t4_name);
			return 0x00;
		}		
	}
	
	char write_t4 = 0;
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			if (fabs(t4trans[i][j]) != 0.0f && (t4trans[i][j] != 1.0f || i != j))
			{
				write_t4 = 1;
			}
		}
	}
	
	if (write_t4 && save_center)
	{
		printf("writing t4 file: %s\n", t4_name);
		for (i = 0; i < 3; ++i)
		{
			printf("t4: %f\t%f\t%f\t%f\n", t4trans[i][0], t4trans[i][1], t4trans[i][2], t4trans[i][3]);
		}
		for (i = 0; i < 3; ++i)
		{/* transpose so the columns are the second index */
			for (j = i + 1; j < 4; ++j)
			{
				tempf = t4trans[i][j];
				t4trans[i][j] = t4trans[j][i];
				t4trans[j][i] = tempf;
			}
		}
		FILE * t4_file = fopen(t4_name, "w");
		if (!t4_file)
		{
			fprintf(stderr, "error opening t4 file '%s': %s\n", t4_name, strerror(errno));
			free_file_names(header_name, image_name, t4_name);
			return 0x00;
		}
		fprintf(t4_file, "%s\n", rcsid);
		fprintf(t4_file, "%s", argv[0]);
		for (i = 1; i < argc; ++i) fprintf(t4_file, " %s", argv[i]);
		fprintf(t4_file, "\n");
		t4_write(t4_file, antiwarning, 0);/* scale has already been applied */
		fclose(t4_file);
	}
	
	/*int image_file = open(image_name, O_RDWR | O_CREAT | O_TRUNC, 0644);
	if (image_file < 0)
	{
		fprintf(stderr, "error opening image '%s': %s\n", image_name, strerror(errno));
		free_file_names(header_name, image_name, t4_name);
		return 0x00;
	}//*/
	FILE * image_stream = fopen(image_name, "wb");
	if (image_stream == NULL)
	{
		fprintf(stderr, "error opening image '%s'\n", image_name);
		free_file_names(header_name, image_name, t4_name);
		return 0x00;
	}
	//begin non-mmap writing code
	void* outmem = malloc(iImage->size * sizeof(float));//needs to be rewritten to use frame at a time, eventually
	if (!outmem)
	{
		fprintf(stderr, "error allocating memory: %s\n", strerror(errno));
		free_file_names(header_name, image_name, t4_name);
		return 0x00;
	}
	auto_orient(outmem, iImage, reswap);
	//ssize_t written = write(image_file, outmem, iImage->size * sizeof(float));
	size_t written = fwrite(outmem, sizeof(float), iImage->size, image_stream);
	free(outmem);
	if (written != (size_t)iImage->size)
	{
		fprintf(stderr, "could not write the entire file '%s'\n", image_name);
		free_file_names(header_name, image_name, t4_name);
		fclose(image_stream);
		return 0x00;
	}
	free_file_names(header_name, image_name, t4_name);
	fclose(image_stream);
	return 1;
	//old mmap code - problematic, and gave no advantage
	/*unsigned char buffer[4096] = {};
	int64_t total = 0;
	int64_t target = iImage->size * 4;
	int64_t page_size = sysconf(_SC_PAGESIZE);
	int64_t rounded_size = target;
	if (rounded_size % page_size) rounded_size += page_size + (-rounded_size % page_size);
	
	while (total < target)
	{
		int write_size = (total + (int)sizeof buffer > target)? (target - total) : (int)sizeof buffer;
		
		if ((int) write(image_file, buffer, write_size) != write_size)
		{
			fprintf(stderr, "error writing image '%s': %s\n", image_name,
					strerror(errno));
			free(header_name);
			close(image_file);
			return 0x00;
		}
		else total += write_size;
	}

	void *output_image = mmap(NULL, rounded_size, PROT_READ | PROT_WRITE, MAP_SHARED, image_file, 0x00);
	if (output_image == -1)
	{
		fprintf(stderr, "error writing image '%s': %s\n", image_name,
				strerror(errno));
		free(header_name);
		close(image_file);
		return 0x00;
	}
	
	auto_orient(output_image, iImage, reswap);
	msync(output_image, target, MS_SYNC);
	munmap(output_image, rounded_size);
	
	free_file_names(header_name, image_name, t4_name);
	close(image_file);
	return 0x01;//*/
}
