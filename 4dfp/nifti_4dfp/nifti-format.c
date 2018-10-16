/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/nifti-format.c,v 1.14 2012/04/23 18:05:06 coalsont Exp $*/
/*$Log: nifti-format.c,v $
 * Revision 1.14  2012/04/23  18:05:06  coalsont
 * fix output of pixdim to ensure it matches sform when oblique
 * check pixdim versus sform on reading
 *
 * Revision 1.13  2012/02/09  00:28:40  coalsont
 * changed write() to fwrite() to get around 2GB limitation
 *
 * Revision 1.12  2011/12/02  05:09:26  coalsont
 * changed file writing to not use mmap
 *
 * Revision 1.11  2010/11/24  03:29:33  coalsont
 * correct for rounding errors in qform to sform conversion, set NaNs to 0
 *
 * Revision 1.10  2010/08/16  18:18:44  coalsont
 * used 64bit integers to handle file sizes
 *
 * Revision 1.9  2010/07/02  16:00:11  coalsont
 * added int cast to remove warning on linux64
 *
 * Revision 1.8  2010/06/18  21:00:59  coalsont
 * moved definitions around to reduce dependencies
 *
 * Revision 1.7  2009/08/20  20:14:21  coalsont
 * handle malformed nifti more gracefully
 * code cleanup
 *
 * Revision 1.6  2009/08/17  22:57:44  coalsont
 * mmap_offset so munmap uses the original pointer
 *
 * Revision 1.5  2009/08/12  15:57:13  coalsont
 * fixed crash: set the mmaped flag in the common_format
 *
 * Revision 1.4  2009/08/07  22:11:43  coalsont
 * check file size to return null instead of crashing
 *
 * Revision 1.2  2009/07/29  20:02:28  coalsont
 * format string changes for t4 file
 *
 * Revision 1.1  2009/07/08  01:01:35  coalsont
 * Initial revision
 **/
//static char* rcsid = "$Id: nifti-format.c,v 1.14 2012/04/23 18:05:06 coalsont Exp $";
/* Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009 */
/* Kevin P. Barry [ta0kira@users.berlios.de], BMCLAB, 22 May 2009 */

#include "nifti-format.h"

#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <libgen.h>
#include <math.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <rec.h>

#include "common-format.h"
#include "transform.h"


#define HEADER_EXT ".nii"

static void get_file_name(const char *fFile, char **hHeader)
{
	size_t original_length = strlen(fFile);
	int extension_index = original_length;
	
	
	if (original_length > strlen(HEADER_EXT) &&
		strcmp(HEADER_EXT, fFile + original_length - strlen(HEADER_EXT)) == 0)
		extension_index = original_length - strlen(HEADER_EXT);
	
	
	*hHeader = calloc(1, extension_index + strlen(HEADER_EXT) + 1);
	strncpy(*hHeader, fFile, extension_index);
	strncpy(*hHeader + extension_index, HEADER_EXT, strlen(HEADER_EXT));
}

static void flip_header(struct nifti_1_header *hHeader)
{
	flip_endian(&hHeader->sizeof_hdr, sizeof(int32_t), 1);
	flip_endian(&hHeader->extents, sizeof(int32_t), 1);
	flip_endian(&hHeader->session_error, sizeof(int16_t), 1);
	flip_endian(&hHeader->dim, sizeof(int16_t), 8);
	flip_endian(&hHeader->intent_p1, sizeof(float), 1);
	flip_endian(&hHeader->intent_p2, sizeof(float), 1);
	flip_endian(&hHeader->intent_p3, sizeof(float), 1);
	flip_endian(&hHeader->intent_code, sizeof(int16_t), 1);
	flip_endian(&hHeader->datatype, sizeof(int16_t), 1);
	flip_endian(&hHeader->bitpix, sizeof(int16_t), 1);
	flip_endian(&hHeader->slice_start, sizeof(int16_t), 1);
	flip_endian(&hHeader->pixdim, sizeof(float), 8);
	flip_endian(&hHeader->vox_offset, sizeof(float), 1);
	flip_endian(&hHeader->scl_slope, sizeof(float), 1);
	flip_endian(&hHeader->scl_inter, sizeof(float), 1);
	flip_endian(&hHeader->slice_end, sizeof(int16_t), 1);
	flip_endian(&hHeader->cal_max, sizeof(float), 1);
	flip_endian(&hHeader->cal_min, sizeof(float), 1);
	flip_endian(&hHeader->slice_duration, sizeof(float), 1);
	flip_endian(&hHeader->toffset, sizeof(float), 1);
	flip_endian(&hHeader->glmax, sizeof(int32_t), 1);
	flip_endian(&hHeader->glmin, sizeof(int32_t), 1);
	flip_endian(&hHeader->qform_code, sizeof(int16_t), 1);
	flip_endian(&hHeader->sform_code, sizeof(int16_t), 1);
	flip_endian(&hHeader->quatern_b, sizeof(float), 1);
	flip_endian(&hHeader->quatern_c, sizeof(float), 1);
	flip_endian(&hHeader->quatern_d, sizeof(float), 1);
	flip_endian(&hHeader->qoffset_x, sizeof(float), 1);
	flip_endian(&hHeader->qoffset_y, sizeof(float), 1);
	flip_endian(&hHeader->qoffset_z, sizeof(float), 1);
	flip_endian(&hHeader->srow_x, sizeof(float), 4);
	flip_endian(&hHeader->srow_y, sizeof(float), 4);
	flip_endian(&hHeader->srow_z, sizeof(float), 4);
}

static int auto_flip(struct nifti_1_header *hHeader)
{
	if (hHeader->dim[0] <= 7 && hHeader->dim[0] > 0) return 0x00;
	flip_header(hHeader);
	return 0x01;
}


common_format *parse_nifti(const char *fFile, char rec_open)
{
	char *header_name = NULL;
	char mystr[1024];
	get_file_name(fFile, &header_name);
	
	int header_file = open(header_name, O_RDONLY);
	int i, j;
	if (header_file < 0)
	{
		fprintf(stderr, "error opening image '%s': %s\n", header_name,
				strerror(errno));
		free(header_name);
		return NULL;
	}
	
	struct stat file_stats;
	if (fstat(header_file, &file_stats) != 0)
	{
		fprintf(stderr, "error reading image '%s': %s\n", header_name,
				strerror(errno));
		free(header_name);
		close(header_file);
		return NULL;
	}
	
	int64_t file_size = file_stats.st_size;
	int64_t rounded_size = file_size;
	int64_t page_size = sysconf(_SC_PAGESIZE);
	if (rounded_size % page_size) rounded_size += page_size + (-rounded_size % page_size);
	
	
	if (file_size < (int)sizeof(header_block))
	{
		fprintf(stderr, "header '%s' is incomplete\n", header_name);
		free(header_name);
		close(header_file);
		return NULL;
	}
	
	
	void *image_data = mmap(NULL, rounded_size, PROT_READ, MAP_SHARED, header_file, 0x00);
	//close(header_file);
	if (!image_data)
	{
		fprintf(stderr, "error reading image '%s': %s\n", header_name,
				strerror(errno));
		free(header_name);
		return NULL;
	}
	
	
	struct nifti_1_header modified_header;
	memcpy(&modified_header, image_data, sizeof(struct nifti_1_header));
	
	
	int flipped = auto_flip(&modified_header);
	
	if (modified_header.sizeof_hdr != sizeof(struct nifti_1_header))
	{
		fprintf(stderr, "WARNING: incorrect header '%s', 'sizeof_hdr' is not %i\n", header_name, (int)sizeof(struct nifti_1_header));
	}
	if (modified_header.dim[0] < 3 || modified_header.dim[0] > 7)
	{
		fprintf(stderr, "ERROR: unusable header '%s', 'dim[0]' is not in the usable range 3-7\n", header_name);
		free(header_name);
		munmap(image_data, rounded_size);
		return NULL;
	}
	if (strcmp(modified_header.magic, "n+1") != 0 )
	{
		fprintf(stderr, "WARNING: incorrect header '%s', 'magic' is not \"n+1\\0\"\n", header_name);
	}
	
	
	/*if (modified_header.datatype != DT_FLOAT32 || modified_header.bitpix != 32)
	{
		fprintf(stderr, "%s: can only convert 32-bit float ('%s')\n", program_name, header_name);
		fprintf(stderr, "bitpix: %i\n", modified_header.bitpix);
		free(header_name);
		close(header_file);
		return NULL;
	}*/
	
	
	common_format *file_content = calloc(1, sizeof(common_format));
	if (!file_content)
	{
		fprintf(stderr, "allocation error: %s\n", strerror(errno));
		free(header_name);
		munmap(image_data, rounded_size);
		return NULL;
	}
	
	file_content->is_mmaped = 1;//signifies not malloced
	
	file_content->endian = local_endian() ? !flipped : flipped;
	
	file_content->f_offset = 0;//file position tracker for eread
	
	file_content->file = header_name;
	
	
	file_content->dims = 4;
	
	printf("datatype: %i\t bitpix: %i\n", modified_header.datatype, modified_header.bitpix);
	if (rec_open)
	{
		sprintf(mystr, "datatype: %i\t bitpix: %i\n", modified_header.datatype, modified_header.bitpix);
		printrec(mystr);
	}
	switch (modified_header.datatype)
	{
		case DT_UINT8:
			file_content->bperpix = 1;
			file_content->isFloat = 0;
			file_content->isSigned = 0;
			break;
		case DT_UINT16:
			file_content->bperpix = 2;
			file_content->isFloat = 0;
			file_content->isSigned = 0;
			break;
		case DT_UINT32:
			file_content->bperpix = 4;
			file_content->isFloat = 0;
			file_content->isSigned = 0;
			break;
		case DT_UINT64:
			file_content->bperpix = 8;
			file_content->isFloat = 0;
			file_content->isSigned = 0;
			break;
		case DT_INT8:
			file_content->bperpix = 1;
			file_content->isFloat = 0;
			file_content->isSigned = 1;
			break;
		case DT_INT16:
			file_content->bperpix = 2;
			file_content->isFloat = 0;
			file_content->isSigned = 1;
			break;
		case DT_INT32:
			file_content->bperpix = 4;
			file_content->isFloat = 0;
			file_content->isSigned = 1;
			break;
		case DT_INT64:
			file_content->bperpix = 8;
			file_content->isFloat = 0;
			file_content->isSigned = 1;
			break;
		case DT_FLOAT32:
			file_content->bperpix = 4;
			file_content->isFloat = 1;
			file_content->isSigned = 1;
			break;
		case DT_FLOAT64:
			file_content->bperpix = 8;
			file_content->isFloat = 1;
			file_content->isSigned = 1;
			break;
		case DT_FLOAT128:
			file_content->bperpix = 16;
			file_content->isFloat = 1;
			file_content->isSigned = 1;
			break;
		default:
			fprintf(stderr, "datatype unsupported ('%s'), value: %i\n", file_content->file, modified_header.datatype);
			free(header_name);
			munmap(image_data, rounded_size);
			free(file_content);
			return NULL;
	};
	if (modified_header.bitpix != 8 * file_content->bperpix)
	{
		fprintf(stderr, "bitpix and datatype do not match ('%s'), using dataype\n", file_content->file);
	}
	
	file_content->scale = modified_header.scl_slope;
	file_content->offset = modified_header.scl_inter;
	if (file_content->scale == 0.0)/* its in nifti1.h, wasn't my idea */
	{
		file_content->scale = 1.0;
		file_content->offset = 0.0;
	}
	
	if (rec_open)
	{
		sprintf(mystr, "scale: %f\t offset: %f\n", file_content->scale, file_content->offset);
		printrec(mystr);
	}
	
	file_content->length[0] = modified_header.dim[1 + 0];
	file_content->length[1] = modified_header.dim[1 + 1];
	file_content->length[2] = modified_header.dim[1 + 2];
	file_content->length[3] = modified_header.dim[1 + 3];
	if (modified_header.dim[0] >= 5 && modified_header.dim[1 + 4]) file_content->length[3] *= modified_header.dim[1 + 4];
	if (modified_header.dim[0] >= 6 && modified_header.dim[1 + 5]) file_content->length[3] *= modified_header.dim[1 + 5];
	if (modified_header.dim[0] >= 7 && modified_header.dim[1 + 6]) file_content->length[3] *= modified_header.dim[1 + 6];
	
	if (file_content->length[3] != modified_header.dim[1 + 3])
	{
		fprintf(stderr, "merging additional dimmensions with time ('%s')\n", file_content->file);
		if (rec_open)
		{
			printrec("additional dimensions merged with time\n");
		}
	}
	
	if (file_content->length[3] == 0) file_content->length[3] = 1;
	
	const vindex axes[4] = AXES_NII;
	
	memcpy(file_content->order, axes, sizeof file_content->order);
	
	file_content->size = 1;
	if (file_content->length[ file_content->order[0] ]) file_content->size *= file_content->length[ file_content->order[0] ];
	if (file_content->length[ file_content->order[1] ]) file_content->size *= file_content->length[ file_content->order[1] ];
	if (file_content->length[ file_content->order[2] ]) file_content->size *= file_content->length[ file_content->order[2] ];
	if (file_content->length[ file_content->order[3] ]) file_content->size *= file_content->length[ file_content->order[3] ];
	
	if ((int) modified_header.vox_offset < (int)sizeof(header_block))
	{
		fprintf(stderr, "WARNING: incorrect header '%s', 'vox_offset' is less than 352\nrecalculating by file size...\n", header_name);
		modified_header.vox_offset = file_size - (file_content->size) * (file_content->bperpix);
	}
	
	if (modified_header.vox_offset < 348 || file_size < file_content->size * file_content->bperpix + modified_header.vox_offset)
	{
		fprintf(stderr, "file length is too short for specified dimensions\n");
		free(header_name);
		munmap(image_data, rounded_size);
		free(file_content);
		return NULL;
	}
	
	file_content->voxels = image_data + (int) modified_header.vox_offset;
	
	file_content->mmap_offset = (int)modified_header.vox_offset;
	
	for (i = 0; i < 3; ++i)
		for (j = 0; j < 4; ++j)
			file_content->sform[i][j] = 0.0;
	
	file_content->haveCenter = 1;

	if ( modified_header.sform_code )
	{/* prefer sform, due to simplicity */
		for (i = 0; i < 4; ++i)
		{
			file_content->sform[0][i] = modified_header.srow_x[i];
			file_content->sform[1][i] = modified_header.srow_y[i];
			file_content->sform[2][i] = modified_header.srow_z[i];
		}
		/* sanity check pixdim against sform, because FSL can use both in flirt */
		const float TOLERANCE = 0.999f;
		for (i = 0; i < 3; ++i)
		{
			float calcspace = 0.0f;
			for (j = 0; j < 3; ++j)
			{
				calcspace += file_content->sform[j][i] * file_content->sform[j][i];
			}
			calcspace = sqrt(calcspace);
			float headerspace = modified_header.pixdim[i + 1];
			if (calcspace != headerspace && (calcspace == 0.0f || headerspace == 0.0f || calcspace / headerspace < TOLERANCE || headerspace / calcspace < TOLERANCE))
			{
				fprintf(stderr, "WARNING: nifti pixdims do not match sform, other software (FSL) may interpret this volume differently.\n");
				if (rec_open)
				{
					printrec("nifti pixdim and sform mismatch!");
				}
			}
		}
	} else {
		if (modified_header.qform_code)
		{/* quaternions may be compact for rotation, but conceptually obtuse, and ugly to convert, and requires separate scale (when constrained to magnitude 1) and offset */
			double b = modified_header.quatern_b, c = modified_header.quatern_c, d = modified_header.quatern_d;
			double a = sqrt(1.0 - (b * b + c * c + d * d));/* given implementation uses 10 floats instead of 12, so it has questionable value */
			if (a != a || a > 1.0 || a < 0.0) a = 0.0;/* make sure a isn't NaN or something equally bad */
			file_content->sform[0][0] = a * a + b * b - c * c - d * d;/* TODO: try to find a more elegant way to convert this */
			file_content->sform[0][1] = 2 * b * c - 2 * a * d;/* formulas copypasted from nifti1.h */
			file_content->sform[0][2] = 2 * b * d + 2 * a * c;/* reused products could be precalculated, but left alone for what remains of clarity */
			file_content->sform[1][0] = 2 * b * c + 2 * a * d;
			file_content->sform[1][1] = a * a + c * c - b * b - d * d;
			file_content->sform[1][2] = 2 * c * d - 2 * a * b;
			file_content->sform[2][0] = 2 * b * d - 2 * a * c;
			file_content->sform[2][1] = 2 * c * d + 2 * a * b;
			file_content->sform[2][2] = a * a + d * d - c * c - b * b;/* rotation matrix complete */
			if (rec_open)
			{
				sprintf(mystr, "qform: %12.6f\t%12.6f\t%12.6f\t%12.6f\n", a, b, c, d);
				printrec(mystr);
				sprintf(mystr, "qfac: %12.6f\n", modified_header.pixdim[0] < 0.0 ? -1.0f : 1.0f);
				printrec(mystr);
				sprintf(mystr, "pixdim: %12.6f\t%12.6f\t%12.6f\n", modified_header.pixdim[1], modified_header.pixdim[2], modified_header.pixdim[3]);
				printrec(mystr);
				sprintf(mystr, "qoffset: %12.6f\t%12.6f\t%12.6f\n", modified_header.qoffset_x, modified_header.qoffset_y, modified_header.qoffset_z);
				printrec(mystr);
			}
			/* quaternion does pixdim[] and qfac (pixdim[0]) first, then rotate (R), then qoffset */
			if (modified_header.pixdim[0] < 0.0) modified_header.pixdim[3] = -modified_header.pixdim[3];/* read nifti1.h:1005 for yourself, i can't make this stuff up */
			/* quaternions can't encode an axis handedness shift (maybe with a negative real component), so they use a float as a boolean and make it a special case */
			for (i = 0; i < 4; ++i)
			{
				file_content->sform[i][0] *= modified_header.pixdim[1];
				file_content->sform[i][1] *= modified_header.pixdim[2];
				file_content->sform[i][2] *= modified_header.pixdim[3];
			}
			file_content->sform[0][3] = modified_header.qoffset_x;
			file_content->sform[1][3] = modified_header.qoffset_y;
			file_content->sform[2][3] = modified_header.qoffset_z;
		} else {/* do it the originless, analyze 7.5 way */
			file_content->sform[0][0] = modified_header.pixdim[1];
			file_content->sform[1][1] = modified_header.pixdim[2];
			file_content->sform[2][2] = modified_header.pixdim[3];
			file_content->haveCenter = 0;
			printrec("no origin information, spacing only\n");
		}
	}
	
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			if (file_content->sform[i][j] != file_content->sform[i][j])
			{
				file_content->sform[i][j] = 0.0;/* set NaNs to zero */
			}
		}
	}
	
	/* rec file */
	for (i = 0; i < 3; ++i)
	{
		printf("common:\t%f\t%f\t%f\t%f\n", file_content->sform[i][0], file_content->sform[i][1], file_content->sform[i][2], file_content->sform[i][3]);
		if (rec_open)
		{
			sprintf(mystr, "common sform[%i]: %12.6f\t%12.6f\t%12.6f\t%12.6f\n", i, file_content->sform[i][0], file_content->sform[i][1], file_content->sform[i][2], file_content->sform[i][3]);
			printrec(mystr);
		}
	}
	
	return file_content;
}

int save_nifti(common_format *iImage, const char *fFile, char* program_name, char control)
{
	char *header_name = NULL;
	
	char reswap = ((local_endian() && (control == 'l' || control == 'L')) || (!local_endian() && (control == 'b' || control == 'B')));
	
	get_file_name(fFile, &header_name);
	int i, j;

	to_lpi(iImage);
	
	auto_orient_header(iImage);
	vindex revorder[4]; /* reverse lookup to reorder length field for axis swapping */
	for (i = 0; i < 4; ++i)
	{
		revorder[iImage->order[i]] = i;
	}
	for (i = 0; i < 3; ++i)
	{
		printf("out:\t%f\t%f\t%f\t%f\n", iImage->sform[i][0], iImage->sform[i][1], iImage->sform[i][2], iImage->sform[i][3]);
	}
	
	float spacing[3] = { 0.0f, 0.0f, 0.0f };//calculate the values to put in pixdim
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			spacing[i] += iImage->sform[j][i] * iImage->sform[j][i];
		}
		spacing[i] = sqrt(spacing[i]);
	}
	
	header_block image_header = {
		.header = {
			.sizeof_hdr = sizeof(struct nifti_1_header),
			
			/*.data_type = ,*/
			/*.db_name = ,*/
			/*.extents = ,*/
			/*.session_error = ,*/
			/*.regular = ,*/
			/*.dim_info = ,*/
			
			.dim = {
				iImage->dims,
				iImage->length[ revorder[0] ],
				iImage->length[ revorder[1] ],
				iImage->length[ revorder[2] ],
				iImage->length[ revorder[3] ]
			},
			
			/*.intent_p1 = ,*/
			/*.intent_p2 = ,*/
			/*.intent_p3 = ,*/
			
			/*.intent_code = ,*/
			.datatype = DT_FLOAT32,
			.bitpix = 32,
			/*.slice_start = ,*/
			.pixdim = {
				1.0,//this is NOT the dimension of anything
				spacing[0],
				spacing[1],
				spacing[2],
				iImage->timestep
			},
			.vox_offset = (float) sizeof(header_block),
			/*.scl_slope = ,*/
			/*.scl_inter = ,*/
			/*.slice_end = ,*/
			/*.slice_code = ,*/
			.xyzt_units= SPACE_TIME_TO_XYZT(NIFTI_UNITS_MM, NIFTI_UNITS_SEC),
			/*.cal_max = ,*/
			/*.cal_min = ,*/
			/*.slice_duration = ,*/
			/*.toffset = ,*/
			/*.glmax = ,*/
			/*.glmin = ,*/
			
			/*.descrip = ,*/
			/*.aux_file = ,*/
			
			/*.qform_code = ,*/
			.sform_code = NIFTI_XFORM_TALAIRACH,
			
			/*.quatern_b = ,*/
			/*.quatern_c = ,*/
			/*.quatern_d = ,*/
			/*.qoffset_x = ,*/
			/*.qoffset_y = ,*/
			/*.qoffset_z = ,*/
			
			.srow_x = {
				iImage->sform[0][0], iImage->sform[0][1], iImage->sform[0][2],
				iImage->sform[0][3]
			},
			
			.srow_y = {
				iImage->sform[1][0], iImage->sform[1][1], iImage->sform[1][2],
				iImage->sform[1][3]
			},
			
			.srow_z = {
				iImage->sform[2][0], iImage->sform[2][1], iImage->sform[2][2],
				iImage->sform[2][3]
			},
			
			/*.intent_name = ,*/
			
			.magic = "n+1"
		},
		
		.extender = {
			.extension = {}
		}
	};
	
	
	/*snprintf(image_header.header.descrip, sizeof image_header.header.descrip,
			 "%s converted with %s", basename(iImage->file), program_name);*/
	
	char* tempstr = basename(iImage->file);
	for (i = 0; i < (int)(sizeof image_header.header.descrip - 1) && tempstr[i]; ++i)
	{
		image_header.header.descrip[i] = tempstr[i];
	}
	tempstr = " converted with ";
	for (j = 0; i < (int)(sizeof image_header.header.descrip - 1) && tempstr[j]; ++i)
	{
		image_header.header.descrip[i] = tempstr[j];
		++j;
	}
	for (j = 0; i < (int)(sizeof image_header.header.descrip - 1) && program_name[j]; ++i)
	{
		image_header.header.descrip[i] = program_name[j];
		++j;
	}
	image_header.header.descrip[i] = '\0';
	
	if (reswap) flip_header(&(image_header.header));
	/*int header_file = open(header_name, O_RDWR | O_CREAT | O_TRUNC, 0644);
	if (header_file < 0)
	{
		fprintf(stderr, "error opening image '%s': %s\n", header_name,
				strerror(errno));
		free(header_name);
		return 0x00;
	}//*/
	FILE * header_stream = fopen(header_name, "wb");
	if (header_stream == NULL)
	{
		fprintf(stderr, "error opening image '%s'\n", header_name);
		free(header_name);
		return 0x00;
	}
	if (fwrite(&image_header, 1, sizeof(image_header), header_stream) != sizeof(image_header))
	{
		fprintf(stderr, "error writing header '%s'\n", header_name);
		free(header_name);
		fclose(header_stream);
		return 0x00;
	}
	//begin non-mmap writing code
	void* outmem = malloc(iImage->size * sizeof(float));//needs to be rewritten to use frame at a time, eventually
	if (!outmem)
	{
		fprintf(stderr, "error allocating memory: %s\n", strerror(errno));
		free(header_name);
		return 0x00;
	}
	auto_orient(outmem, iImage, reswap);
	//ssize_t written = write(header_file, outmem, iImage->size * sizeof(float));
	size_t written = fwrite(outmem, sizeof(float), iImage->size, header_stream);
	free(outmem);
	if (written != (size_t)iImage->size)
	{
		fprintf(stderr, "could not write the entire file '%s'\n", header_name);
		free(header_name);
		fclose(header_stream);
		return 0x00;
	}
	free(header_name);
	fclose(header_stream);
	return 1;	
	//old mmap code - problematic, and gave no advantage
	/*unsigned char buffer[4096] = {};
	int64_t total = 0;
	int64_t target = iImage->size * 4;
	int64_t page_size = sysconf(_SC_PAGESIZE);
	int64_t rounded_size = sizeof(image_header) + target;
	if (rounded_size % page_size) rounded_size += page_size + (-rounded_size % page_size);
	
	while (total < target)
	{
		int write_size = (total + (int)sizeof buffer > target)? (target - total) : (int)sizeof buffer;
		
		if ((int) write(header_file, buffer, write_size) != write_size)
		{
			fprintf(stderr, "error writing image '%s': %s\n", header_name,
					strerror(errno));
			free(header_name);
			close(header_file);
			return 0x00;
		}
		else total += write_size;
	}
	
	
	void *output_image = mmap(NULL, rounded_size, PROT_READ | PROT_WRITE, MAP_SHARED, header_file, 0x00);
	if (output_image == -1)
	{
		fprintf(stderr, "error writing image '%s': %s\n", header_name,
				strerror(errno));
		free(header_name);
		close(header_file);
		return 0x00;
	}
	
	auto_orient(output_image + sizeof(image_header), iImage, reswap);
	msync(output_image, target, MS_SYNC);
	munmap(output_image, rounded_size);
	
	free(header_name);
	close(header_file);
	return 0x01;//*/
}

FILE* fopen_nifti(char* fFile, int* dim, float* spacing, float* center, char* program_name, char control)
{
	char *header_name = NULL;
	
	char reswap = ((local_endian() && (control == 'l' || control == 'L')) || (!local_endian() && (control == 'b' || control == 'B')));
	
	get_file_name(fFile, &header_name);
	int i, j;
	FILE* header_file = fopen(header_name, "wb");
	if (header_file == 0)
	{
		fprintf(stderr, "error opening image '%s': %s\n", header_name,
				strerror(errno));
		free(header_name);
		return 0x00;
	}
	header_block image_header = {
		.header = {
			.sizeof_hdr = sizeof(struct nifti_1_header),
			
			/*.data_type = ,*/
			/*.db_name = ,*/
			/*.extents = ,*/
			/*.session_error = ,*/
			/*.regular = ,*/
			/*.dim_info = ,*/
			
			.dim = {
				4,
				dim[0],
				dim[1],
				dim[2],
				dim[3]
			},
			
			/*.intent_p1 = ,*/
			/*.intent_p2 = ,*/
			/*.intent_p3 = ,*/
			
			/*.intent_code = ,*/
			.datatype = DT_FLOAT32,
			.bitpix = 32,
			/*.slice_start = ,*/
			.pixdim = {
				1.0,//this is NOT the dimension of anything
				spacing[0],
				spacing[1],
				spacing[2],
				spacing[3]
			},
			.vox_offset = (float) sizeof(header_block),
			/*.scl_slope = ,*/
			/*.scl_inter = ,*/
			/*.slice_end = ,*/
			/*.slice_code = ,*/
			.xyzt_units= SPACE_TIME_TO_XYZT(NIFTI_UNITS_MM, NIFTI_UNITS_SEC),
			/*.cal_max = ,*/
			/*.cal_min = ,*/
			/*.slice_duration = ,*/
			/*.toffset = ,*/
			/*.glmax = ,*/
			/*.glmin = ,*/
			
			/*.descrip = ,*/
			/*.aux_file = ,*/
			
			/*.qform_code = ,*/
			.sform_code = NIFTI_XFORM_TALAIRACH,
			
			/*.quatern_b = ,*/
			/*.quatern_c = ,*/
			/*.quatern_d = ,*/
			/*.qoffset_x = ,*/
			/*.qoffset_y = ,*/
			/*.qoffset_z = ,*/
			
			.srow_x = {
				spacing[0], 0.0f, 0.0f, center[0]
			},
			
			.srow_y = {
				0.0f, spacing[1], 0.0f, center[1]
			},
			
			.srow_z = {
				0.0f, 0.0f, spacing[2], center[2]
			},
			
			/*.intent_name = ,*/
			
			.magic = "n+1"
		},
		
		.extender = {
			.extension = {}
		}
	};
	
	
	/*snprintf(image_header.header.descrip, sizeof image_header.header.descrip,
	 "%s converted with %s", basename(iImage->file), program_name);*/
	
	char* tempstr = basename(fFile);
	for (i = 0; i < (int)(sizeof image_header.header.descrip - 1) && tempstr[i]; ++i)
	{
		image_header.header.descrip[i] = tempstr[i];
	}
	tempstr = " converted with ";
	for (j = 0; i < (int)(sizeof image_header.header.descrip - 1) && tempstr[j]; ++i)
	{
		image_header.header.descrip[i] = tempstr[j];
		++j;
	}
	for (j = 0; i < (int)(sizeof image_header.header.descrip - 1) && program_name[j]; ++i)
	{
		image_header.header.descrip[i] = program_name[j];
		++j;
	}
	image_header.header.descrip[i] = '\0';
	
	if (reswap) flip_header(&(image_header.header));
	
	if (fwrite(&image_header, sizeof(header_block), 1, header_file) != sizeof(header_block)) return NULL;
	
	return header_file;
}

FILE* fopen_nifti_sform(char* fFile, int* dim, float sform[3][4], float timestep, char* program_name, char control)
{//sform is (float*)(float[3][4])
	char *header_name = NULL;
	
	char reswap = ((local_endian() && (control == 'l' || control == 'L')) || (!local_endian() && (control == 'b' || control == 'B')));
	
	get_file_name(fFile, &header_name);
	int i, j;
	FILE* header_file = fopen(header_name, "wb");
	if (header_file == 0)
	{
		fprintf(stderr, "error opening image '%s': %s\n", header_name,
				strerror(errno));
		free(header_name);
		return 0x00;
	}
	header_block image_header = {
		.header = {
			.sizeof_hdr = sizeof(struct nifti_1_header),
			
			/*.data_type = ,*/
			/*.db_name = ,*/
			/*.extents = ,*/
			/*.session_error = ,*/
			/*.regular = ,*/
			/*.dim_info = ,*/
			
			.dim = {
				4,
				dim[0],
				dim[1],
				dim[2],
				dim[3]
			},
			
			/*.intent_p1 = ,*/
			/*.intent_p2 = ,*/
			/*.intent_p3 = ,*/
			
			/*.intent_code = ,*/
			.datatype = DT_FLOAT32,
			.bitpix = 32,
			/*.slice_start = ,*/
			.pixdim = {
				1.0,//this is NOT the dimension of anything
				sform[0][0],
				sform[1][1],
				sform[2][2],
				timestep
			},
			.vox_offset = (float) sizeof(header_block),
			/*.scl_slope = ,*/
			/*.scl_inter = ,*/
			/*.slice_end = ,*/
			/*.slice_code = ,*/
			.xyzt_units = SPACE_TIME_TO_XYZT(NIFTI_UNITS_MM, NIFTI_UNITS_SEC),
			/*.cal_max = ,*/
			/*.cal_min = ,*/
			/*.slice_duration = ,*/
			/*.toffset = ,*/
			/*.glmax = ,*/
			/*.glmin = ,*/
			
			/*.descrip = ,*/
			/*.aux_file = ,*/
			
			/*.qform_code = ,*/
			.sform_code = NIFTI_XFORM_TALAIRACH,
			
			/*.quatern_b = ,*/
			/*.quatern_c = ,*/
			/*.quatern_d = ,*/
			/*.qoffset_x = ,*/
			/*.qoffset_y = ,*/
			/*.qoffset_z = ,*/
			
			.srow_x = {
				sform[0][0], sform[0][1], sform[0][2], sform[0][3]
			},
			
			.srow_y = {
				sform[1][0], sform[1][1], sform[1][2], sform[1][3]
			},
			
			.srow_z = {
				sform[2][0], sform[2][1], sform[2][2], sform[2][3]
			},
			
			/*.intent_name = ,*/
			
			.magic = "n+1"
		},
		
		.extender = {
			.extension = {}
		}
	};
	
	
	/*snprintf(image_header.header.descrip, sizeof image_header.header.descrip,
	 "%s converted with %s", basename(iImage->file), program_name);*/
	
	char* tempstr = basename(fFile);
	for (i = 0; i < (int)(sizeof image_header.header.descrip - 1) && tempstr[i]; ++i)
	{
		image_header.header.descrip[i] = tempstr[i];
	}
	tempstr = " converted with ";
	for (j = 0; i < (int)(sizeof image_header.header.descrip - 1) && tempstr[j]; ++i)
	{
		image_header.header.descrip[i] = tempstr[j];
		++j;
	}
	for (j = 0; i < (int)(sizeof image_header.header.descrip - 1) && program_name[j]; ++i)
	{
		image_header.header.descrip[i] = program_name[j];
		++j;
	}
	image_header.header.descrip[i] = '\0';
	
	if (reswap) flip_header(&(image_header.header));
	
	if (fwrite(&image_header, sizeof(header_block), 1, header_file) != sizeof(header_block)) return NULL;
	
	return header_file;
}
