/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/transform.c,v 1.6 2012/09/08 01:08:17 avi Exp $*/
/*$Log: transform.c,v $
 * Revision 1.6  2012/09/08  01:08:17  avi
 * eliminate in-line type defining and, hence, requirement for gcc -std=c99
 *
 * Revision 1.5  2010/08/16  18:20:29  coalsont
 * cast to 64bit integer in get_voxel
 *
 * Revision 1.4  2010/06/18  21:00:59  coalsont
 * moved things around to reduce dependencies
 *
 * Revision 1.3  2009/08/17  22:57:44  coalsont
 * debug code
 *
 * Revision 1.2  2009/08/07  20:42:04  coalsont
 * restructuring, endian output control
 *
 * Revision 1.1  2009/07/08  01:01:35  coalsont
 * Initial revision
 **/
static char* rcsid = "$Id: transform.c,v 1.6 2012/09/08 01:08:17 avi Exp $";
/* Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009 */
/* Kevin P. Barry [ta0kira@users.berlios.de], BMCLAB, 22 May 2009 */

#include "transform.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mystdint.h"

typedef uint8_t  width1;
typedef uint16_t width2;
typedef uint32_t width4;
typedef double width8;


/*(borrowed from 'nifti1_io.h' from nifticlib-1.1.0)*/
int local_endian()
{
	union {
		unsigned char byte;
		short         word;
	} test;
	
	test.word = 01;
	
	return !test.byte;
}

void *get_voxel(void *dData, const vindex pPos[4], const vindex length[4], uint8_t bperpix)
{
	return dData + ( pPos[0] + length[0] *
	                (pPos[1] + length[1] *
					 (pPos[2] + (int64_t)length[2] *
					  pPos[3])) ) * bperpix;
}

void flip_endian(void *dData, int wWidth, int nNumber)
{
	int I, J;
	
	for (I = 0; I < nNumber; I++)
		for (J = 0; J < wWidth / 2; J++)
		{
			unsigned char temp = ((unsigned char*) dData)[I * wWidth + J];
			((unsigned char*) dData)[I * wWidth + J] = ((unsigned char*) dData)[(I + 1) * wWidth - J - 1];
			((unsigned char*) dData)[(I + 1) * wWidth - J - 1] = temp;
		}
}

float get_value_float(const void* data, common_format* iImage, int rReverse)
{
	static	uint8_t buffer[16];
	static	uint8_t wWidth;
	static	float ret = 0.0f;
	int	I;

	wWidth = iImage->bperpix;
	if (rReverse)
	{/* won't work with compound types */
		for (I = 0; I < wWidth; I++)
			buffer[I] = ((uint8_t*) data)[wWidth - I - 1];
	} else {
		memcpy(buffer, data, wWidth);
	}
	/* ugly (and possibly slow), but needed because 4dfp doesn't really support anything but float32 */
	if (iImage->isFloat)
	{
		switch (wWidth)
		{
			case 4:
				ret = *(float*)buffer;
				break;
			case 8:
				ret = *(double*)buffer;
				break;
			case 16:
				ret = *(long double*)buffer;
				break;
		};
	} else {
		if (iImage->isSigned)
		{
			switch (wWidth)
			{
				case 1:
					ret = *(int8_t*)buffer;
					break;
				case 2:
					ret = *(int16_t*)buffer;
					break;
				case 4:
					ret = *(int32_t*)buffer;
					break;
				case 8:
					ret = *(int64_t*)buffer;
					break;
			};
		} else {
			switch (wWidth)
			{
				case 1:
					ret = *(uint8_t*)buffer;
					break;
				case 2:
					ret = *(uint16_t*)buffer;
					break;
				case 4:
					ret = *(uint32_t*)buffer;
					break;
				case 8:
					ret = *(uint64_t*)buffer;
					break;
			};
		}
	}/* all other cases are caught on the parsing end */
	
	/* apply scale instead of encoding it in the header/t4 */
	if (iImage->scale != 1.0 || iImage->offset != 0.0)
	{
		ret = ret * iImage->scale + iImage->offset;
	}
	return ret;
}

int map_voxel(void *lLeft, const void *rRight, common_format *iImage, int rReverse, char reswap)
{
	int	i, j, ret;
	int8_t	buffer[4];

	*(float*)buffer = get_value_float(rRight, iImage, rReverse);
	ret = !(*(float*) buffer == *(float*) buffer); /*NaN check*/

	if (reswap)
	{
		j = 3;
		for (i = 0; i < 4; ++i, --j)
		{
			((int8_t*)lLeft)[i] = buffer[j];
		}
	} else {
		*(float*)lLeft = *(float*)buffer;
	}
	return ret;
}


void to_lpi(common_format *iImage)
{
	double	max;
	char	used = 0;
	int	i, j, k;

	for (i = 0; i < 2; ++i)
	{
		max = -1.0;
		k = -1;
		for (j = 0; j < 3; ++j)
		{
			if ((!(used & (1<<j))) && fabs(iImage->sform[j][i]) > max)
			{
				max = fabs(iImage->sform[j][i]);
				k = j;
			}
		}
		used |= (1<<k);
		iImage->order[i] = k;
	}
	switch (used)
	{
		case 3:
			iImage->order[2] = 2;
			break;
		case 5:
			iImage->order[2] = 1;
			break;
		case 6:
			iImage->order[2] = 0;
	};/* length[0] refers to fastest moving index */
	
	/* major axis order is now specified by ->order, sform first row is fastest index, ->order says what axis moves the most */
	iImage->orientation = 0;
	for (i = 0; i < 3; ++i)
	{
		if (iImage->sform[iImage->order[i]][i] < 0.0)
		{
			iImage->orientation ^= 01 << i;
		}
	}
}

void auto_orient_header(common_format *iImage)
{/* header needs to be written before calling auto_orient for nifti, so we fix the header first */
	double	sform[3][4];
	int	i, j;

	for (i = 0; i < 3; ++i)
	{/* flip axes */
		if (iImage->orientation & (01 << i))
		{
			for (j = 0; j < 3; ++j)
			{
				iImage->sform[j][3] = ((iImage->length[i] - 1) * iImage->sform[j][i] + iImage->sform[j][3]);
				iImage->sform[j][i] = -iImage->sform[j][i];
			}
		}
	}
	for (i = 0; i < 3; ++i)
	{/* reorder axes to x, y, z[, t] */
		for (j = 0; j < 3; ++j)
		{
			sform[i][iImage->order[j]] = iImage->sform[i][j];
		}
	}
	for (i = 0; i < 3; ++i)
	{/* load it back into the sform */
		for (j = 0; j < 3; ++j)
		{
			iImage->sform[i][j] = sform[i][j];
		}
	}/* NOTE: length is still from the original encoding */
	/* leave order/length as is, because auto_orient needs it */
}

void auto_orient(void *dData, common_format *iImage, char reswap)
{
	int	swap_required, nan_found = 0;
	int	i;
	int	val_flip[4];

	swap_required = !compare_endian(iImage->endian);
	vindex in_val[4] = {}, out_val[4], target_length[4]; /* encode both in nifti axes order, 4dfp orientation 2 allows us to do this with some flips */
	for (i = 0; i < 4; ++i)
	{
		target_length[ iImage->order[i] ] = iImage->length[i];
		val_flip[i] = iImage->orientation & (01 << i);
	}
	printf("flip\tin\tout\n");
	for (i = 0; i < 4; ++i)
	{
		printf("%s\t%i\t%i\n", val_flip[i] ? "yes" : "no", iImage->length[i], target_length[i]);
	}
	/* scan linearly through input file to allow readahead */
	/* could be read incrementally from a file pointer instead of mmap needing it to fit in address space */
	//char print = 1;//for debugging
	//char axes[] = "ijks";
	for (in_val[3] = 0; in_val[3] < iImage->length[3]; in_val[3]++)
	{
		out_val[ iImage->order[3] ] = val_flip[3]?
		(iImage->length[3] - in_val[3] - 1) : in_val[3];
		//if (print) printf("in(%c): %i\t->\tout(%c): %i\n", axes[3], in_val[3], axes[iImage->order[3]], out_val[iImage->order[3]]);
		/* recalculate corresponding output index only when it changes
		 * using iImage->order here means it doesn't need to be passed to get_voxel */
		for (in_val[2] = 0; in_val[2] < iImage->length[2]; in_val[2]++)
		{
			out_val[ iImage->order[2] ] = val_flip[2]?
			(iImage->length[2] - in_val[2] - 1) : in_val[2];
			//if (print) printf("in(%c): %i\t->\tout(%c): %i\n", axes[2], in_val[2], axes[iImage->order[2]], out_val[iImage->order[2]]);
			
			for (in_val[1] = 0; in_val[1] < iImage->length[1]; in_val[1]++)
			{
				out_val[ iImage->order[1] ] = val_flip[1]?
				(iImage->length[1] - in_val[1] - 1) : in_val[1];
				//if (print) printf("in(%c): %i\t->\tout(%c): %i\n", axes[1], in_val[1], axes[iImage->order[1]], out_val[iImage->order[1]]);
				
				for (in_val[0] = 0; in_val[0] < iImage->length[0]; in_val[0]++)
				{
					out_val[ iImage->order[0] ] = val_flip[0]?
					(iImage->length[0] - in_val[0] - 1) : in_val[0];
					//if (print) printf("in(%c): %i\t->\tout(%c): %i\n", axes[0], in_val[0], axes[iImage->order[0]], out_val[iImage->order[0]]);
					//print = 0;
					/* map_voxel(to, from), returns true if the voxel has a value of NaN */
					if (map_voxel(get_voxel(dData, out_val, target_length, 4),
								  get_voxel(iImage->voxels, in_val, iImage->length, iImage->bperpix),
								  iImage, swap_required, reswap))
					{
						if (!nan_found)
						{
							nan_found = 1;
							fprintf(stderr, "NaN value found, possible incorrect endianness in image '%s'\n", iImage->file);
						}
					}
				}
			}
		}
	}
}

void invert_affine(float affineIn[3][4], float affineOut[3][4])
{
	int	a, b, c, d, i, j;
	double	temp, temp2, det = 0.0;

	for (i = 0; i < 3; ++i)
	{/* determinant, can use diagonals trick for 3x3 */
		temp = 1.0;
		temp2 = 1.0;
		affineOut[i][3] = 0.0f;
		for (j = 0; j < 3; ++j)
		{
			temp *= affineIn[j][(i + j) % 3];
			temp2 *= affineIn[j][(i - j + 3) % 3];
		}
		det += temp - temp2;
	}/* adjugate */
	for (i = 0; i < 3; ++i)
	{
		a = (i + 1) % 3;/* wraparound avoids +/- alternation */
		b = (i + 2) % 3;
		for (j = 0; j < 3; ++j)
		{
			c = (j + 1) % 3;
			d = (j + 2) % 3;
			affineOut[j][i] = affineIn[a][c] * affineIn[b][d] - affineIn[a][d] * affineIn[b][c];
		}
	}/* divide by determinant */
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			affineOut[i][j] /= det;
		}
	}/* multiply offset by inverse t4 and negate */
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			affineOut[i][3] -= affineOut[i][j] * affineIn[j][3];
		}
	}
}

