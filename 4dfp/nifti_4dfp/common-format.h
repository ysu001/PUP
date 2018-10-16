/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/common-format.h,v 1.5 2010/08/16 18:20:29 coalsont Exp $*/
/*$Log: common-format.h,v $
 * Revision 1.5  2010/08/16  18:20:29  coalsont
 * use 64bit integer for file size
 *
 * Revision 1.4  2010/06/18  21:00:59  coalsont
 * moved things around to reduce dependencies
 *
 * Revision 1.3  2009/08/17  22:57:44  coalsont
 * common format convienience wrappers
 *
 * Revision 1.2  2009/08/07  20:42:04  coalsont
 * new methods, added field for "file offset" to emulate file pointer
 *
 * Revision 1.1  2009/07/08  00:54:37  coalsont
 * Initial revision
 **/
/* Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009 */
/* Kevin P. Barry [ta0kira@users.berlios.de], BMCLAB, 22 May 2009 */

#ifndef common_format_h
#define common_format_h

#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>

/*#include <stdint.h>*/

/* workaround for nonexistant stdint.h: if find/replace were done, incompatibilities with different OS type sizes would be
 * harder to resolve.  better fix: check for if stdint.h should exist, possibly within mystdint.h. */
#include "mystdint.h"

#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2
#define T_AXIS 3

#define AXES_4DFP { X_AXIS, Y_AXIS, Z_AXIS, T_AXIS }
#define AXES_NII { X_AXIS, Y_AXIS, Z_AXIS, T_AXIS }

typedef int32_t vindex;

typedef struct
{
	char *file;

	void *voxels;
	vindex dims, length[4];
	int64_t size;

	double sform[3][4];/* same as nifti's sform, from index to coordinates */
	double timestep;/* 4dfp doesn't support an "origin" for time, so forget it */
	double scale, offset;/* handle "store values as different values" stuff */
	char isFloat, isSigned, haveCenter, is_mmaped;
	int mmap_offset;//to get the root pointer of the mmap back when parsing nifti

	uint8_t bperpix, orientation;
	int endian; /*0: little, !0: big*/

	/*this specifies axis translations between formats*/
	vindex order[4];
	//keep track of "file position" for eread
	int64_t f_offset;
} common_format;

static inline void free_common_format(common_format *iImage)
{
	if (iImage)
	{
		if (iImage->file) free((void*) iImage->file);
		if (iImage->voxels)
		{
			if (iImage->is_mmaped)
			{
				int64_t page_size = sysconf(_SC_PAGESIZE);
				int64_t rounded_size = ((iImage->size * iImage->bperpix + page_size - 1) / page_size) * page_size;
				munmap(iImage->voxels - iImage->mmap_offset, rounded_size);
			} else {
				free(iImage->voxels);
			}
		}
		free(iImage);
	}
}

static inline void rewind_common_format(common_format* iImage)
{
	iImage->f_offset = 0;
}

static inline void fseek_common_format(common_format* iImage, long offset, int origin)
{
	switch (origin)
	{
		case SEEK_SET:
			iImage->f_offset = offset / 4;//index, NOT bytes
			break;
		case SEEK_CUR:
			iImage->f_offset += offset / 4;
			break;
		case SEEK_END:
			iImage->f_offset = iImage->size + offset / 4;
	}
}

common_format* malloc_common_format(int* dim, float* spacing, float* center, char control);
common_format* malloc_common_format_sform(int* dim, float* sform, float timestep, char control);
int eread_common(float *imgt, int64_t n, int isbig, common_format* myfile);
int ewrite_common(float *imgt, int64_t n, char control, common_format* myfile);//ONLY for malloced, currently

#endif /*common_format_h*/
