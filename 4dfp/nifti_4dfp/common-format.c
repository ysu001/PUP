/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/common-format.c,v 1.5 2012/09/08 01:04:45 avi Exp $*/
/*$Log: common-format.c,v $
 * Revision 1.5  2012/09/08  01:04:45  avi
 * eliminate in-line type defining and, hence, requirement for gcc -std=c99
 *
 * Revision 1.4  2010/08/16  18:20:29  coalsont
 * cast to 64bit integer in get_voxel
 *
 * Revision 1.3  2010/06/18  21:00:59  coalsont
 * moved things around to reduce dependencies
 *
 * Revision 1.2  2009/08/17  22:57:44  coalsont
 * common format convienience wrappers
 *
 * Revision 1.1  2009/08/07  20:42:04  coalsont
 * Initial revision
 **/
/*
 *  common-format.c
 *  
 *
 *  Created by Tim Coalson on 8/7/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


static char* rcsid = "$Id: common-format.c,v 1.5 2012/09/08 01:04:45 avi Exp $";

#include <stdio.h>
#include <string.h>

#include "common-format.h"
#include "transform.h"

common_format* malloc_common_format(int* dim, float* spacing, float* center, char control)
{//only supports 4 dimensions, for true nifti compatibility use the nifti1.1 clib
	//IMPORTANT: "center" is true coordinates at the sample point of voxel 0, simplest convention

	int		i;
	common_format*	ret;

	ret = (common_format*)calloc(1, sizeof(common_format));
	ret->size = 1;
	for (i = 0; i < 3; ++i)
	{
		ret->sform[i][i] = spacing[i];
		ret->sform[i][3] = center[i];
		ret->length[i] = dim[i];
		if (dim[i]) ret->size *= dim[i];
	}
	ret->order[3] = 3;
	ret->endian = local_endian();
	if (control == 'b' || control == 'B') ret->endian = 1;
	if (control == 'l' || control == 'L') ret->endian = 0;
	ret->timestep = spacing[3];
	ret->length[3] = dim[3];
	if (dim[3]) ret->size *= dim[3];
	ret->voxels = calloc(ret->size, 4);
	ret->is_mmaped = 0;
	ret->bperpix = 4;
	ret->isSigned = 1;
	ret->isFloat = 1;
	ret->haveCenter = 1;
	to_lpi(ret);//let user know how the data will actually be saved
	ret->dims = 4;
	ret->f_offset = 0;
	ret->file = NULL;
	ret->scale = 1.0;
	ret->offset = 0.0;
	return ret;
}

common_format* malloc_common_format_sform(int* dim, float* sform, float timestep, char control)
{//only supports 4 dimensions, for true nifti compatibility use the nifti1.1 clib
	//IMPORTANT: fastest index is column

	int		i, j;
	common_format*	ret;

	ret = (common_format*)calloc(1, sizeof(common_format));
	ret->size = 1;
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			ret->sform[i][j] = sform[i * 4 + j];
		}
		ret->length[i] = dim[i];
		if (dim[i]) ret->size *= dim[i];
	}
	ret->endian = local_endian();
	if (control == 'b' || control == 'B') ret->endian = 1;
	if (control == 'l' || control == 'L') ret->endian = 0;
	ret->timestep = timestep;
	ret->voxels = calloc(ret->size, 4);
	ret->is_mmaped = 0;
	ret->bperpix = 4;
	ret->isSigned = 1;
	ret->isFloat = 1;
	ret->haveCenter = 1;
	to_lpi(ret);//let user know how data will actually be saved
	ret->dims = 4;
	ret->f_offset = 0;
	ret->file = NULL;
	ret->scale = 1.0;
	ret->offset = 0.0;
	return ret;
}


int eread_common(float *imgt, int64_t n, int isbig, common_format* myfile)
{
	int64_t i, end;
	int swab_flag;
	
	swab_flag = (local_endian() != 0) != (isbig != 0);
	if (0) printf ("eread swab_flag=%d\n", swab_flag);
	end = myfile->f_offset + n;
	if (end > myfile->length[0] * myfile->length[1] * myfile->length[2] * myfile->length[3]) return -1;
	for (i = myfile->f_offset; i < end; i++)
	{
		imgt[i] = get_value_float(myfile->voxels + i * myfile->bperpix, myfile, swab_flag);
	}
	myfile->f_offset = end;
	return 0;
}

int ewrite_common(float *imgt, int64_t n, char control, common_format* myfile)
{
	int64_t		i, j, k, end, jend;
	int		swab_flag;

	if (myfile->bperpix != 4 || !myfile->isFloat) return -1;
	swab_flag = (local_endian() && (control == 'l' || control == 'L')) || (!local_endian() && (control == 'b' || control == 'B'));
	//maintain endianness in "memory" because it might be mmaped
	if (0) printf ("ewrite swab_flag=%d\n", swab_flag);
	end = myfile->f_offset + n;
	if (end > myfile->length[0] * myfile->length[1] * myfile->length[2] * myfile->length[3]) return -1;
	for (i = myfile->f_offset; i < end; i++)
	{
		if (swab_flag)
		{
			jend = (i + 1) * 4;
			k = jend - 1;
			for (j = i * 4; j < jend; ++j, --k)
			{
				((int8_t*)myfile->voxels)[j] = ((int8_t*)imgt)[k];
			}
		} else {
			((float*)myfile->voxels)[i] = imgt[i];
		}
	}
	myfile->f_offset = end;
	return 0;
}
