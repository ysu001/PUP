/*$Header: /data/petsun4/data1/src_solaris/aff_conv/RCS/aff_conv.c,v 1.9 2012/04/23 19:27:07 coalsont Exp $*/
/*$Log: aff_conv.c,v $
 * Revision 1.9  2012/04/23  19:27:07  coalsont
 * revised help info to clarify world matrices
 *
 * Revision 1.8  2012/04/23  17:55:56  coalsont
 * use sqrt because sqrtf doesn't exist on petsun43
 *
 * Revision 1.7  2012/04/23  17:45:13  coalsont
 * fixed bad loop nesting in new flirt matrix logic
 *
 * Revision 1.6  2012/04/23  16:33:08  coalsont
 * more helpful help info
 * better math for flirt scale space
 *
 * Revision 1.5  2012/04/18  23:34:30  coalsont
 * minor edit to remove c99 syntax
 *
 * Revision 1.4  2012/04/18  21:31:46  coalsont
 * fixed uninitialized values in xr3d function
 *
 * Revision 1.3  2011/11/23  22:40:44  coalsont
 * fixes from testing, improved help info
 *
 * Revision 1.2  2011/11/23  21:18:55  coalsont
 * changed to allow separate original and desired spaces cleanly
 * not yet tested
 **/
#include "stdio.h"
#include "math.h"
#include "string.h"

#include "common-format.h"
#include "nifti-format.h"
#include "4dfp-format.h"
#include "transform.h"

#include "t4_io.h"

static char rcsid[] = "$Id: aff_conv.c,v 1.9 2012/04/23 19:27:07 coalsont Exp $";

void viewmat(float* mat[], int rows, int cols)
{
	int i, j;
	for (i = 0; i < rows; ++i)
	{
		for (j = 0; j < cols; ++j)
		{
			printf("%f\t", mat[i][j]);
		}
		printf("\n");
	}
}

void view4(float mat[4][4])
{
	int i, j;
	for (i = 0; i < 4; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			printf("%f\t", mat[i][j]);
		}
		printf("\n");
	}
}

void write44(FILE* output, float input[4][4])
{
	int i;
	for (i = 0; i < 4; ++i)
	{
		fprintf(output, "%f %f %f %f\n", input[i][0], input[i][1], input[i][2], input[i][3]);
	}
}

void read34(FILE* input, float output[4][4])
{
	int i;
	char* pattern = "%f %f %f %f";
	for (i = 0; i < 3; ++i)
	{
		if (fscanf(input, pattern, &output[i][0], &output[i][1], &output[i][2], &output[i][3]) != 4)
		{
			printf("error reading input affine matrix\n");
			exit(-1);
		}
		output[3][i] = 0.0f;
	}
	output[3][3] = 1.0f;
}

void mult4(float A[4][4], float B[4][4], float output[4][4])
{
	int i, j, k;
	for (i = 0; i < 4; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			output[i][j] = 0.0f;
			for (k = 0; k < 4; ++k)
			{
				output[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void calcrref(float* matrix[], int rows, int cols)
{//assumes more cols than rows
	int i, j, k, temp;
	float tempf, tempf2;
	for (i = 0; i < rows; ++i)
	{
		tempf = fabs(matrix[i][i]);//search for pivot
		temp = i;
		for (j = i + 1; j < rows; ++j)
		{
			tempf2 = fabs(matrix[j][i]);
			if (tempf2 > tempf)
			{
				tempf = tempf2;
				temp = j;
			}
		}
		if (i != temp)
		{
			for (j = i; j < cols; ++j)
			{//skip the waste that will end up 0's and 1's
				tempf = matrix[i][j];
				matrix[i][j] = matrix[temp][j];
				matrix[temp][j] = tempf;
			}
		}
		tempf = matrix[i][i];//pivot
		for (j = i + 1; j < cols; ++j)
		{//again, skip the 0's and 1's
			matrix[i][j] /= tempf;
			for (k = 0; k < i; ++k)
			{
				matrix[k][j] -= matrix[k][i] * matrix[i][j];
			}
			for (++k; k < rows; ++k)
			{
				matrix[k][j] -= matrix[k][i] * matrix[i][j];
			}
		}
	}//rref complete for all cols >= rows, just pretend rowsXrows is I
}

float det3(float input[4][4])
{//use diagonals trick
	float ret = input[0][0] * input[1][1] * input[2][2];
	ret += input[1][0] * input[2][1] * input[0][2];
	ret += input[2][0] * input[0][1] * input[1][2];
	ret -= input[0][2] * input[1][1] * input[2][0];
	ret -= input[1][2] * input[2][1] * input[0][0];
	ret -= input[2][2] * input[0][1] * input[1][0];
	return ret;
}

void inv4(float input[4][4], float output[4][4])
{//use gaussian elimination
	int i, j;
	float* array[4];
	for (i = 0; i < 4; ++i)
	{
		array[i] = (float*)malloc(8 * sizeof(float));
		for (j = 0; j < 4; ++j)
		{
			array[i][j] = input[i][j];
		}
		for (j = 4; j < 8; ++j)
		{
			array[i][j] = 0.0f;
		}
		array[i][i + 4] = 1.0f;
	}
	calcrref(array, 4, 8);//lazy rref
	for (i = 0; i < 4; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			output[i][j] = array[i][j + 4];
		}
		free(array[i]);
	}
}

void parse_t4(char* input, common_format* source, common_format* target, float world[4][4])
{
	float invmat[4][4], flat[16];
	t4_read(input, flat);
	int i, j;
	for (i = 0; i < 4; ++i)
	{
		for (j = 0; j < 4; ++j)
		{//order[3] should always be 3, no worries
			invmat[source->order[i]][target->order[j]] = flat[i + j * 4];
		}
	}
	inv4(invmat, world);
}

void parse_world(FILE* input, float world[4][4])
{
	read34(input, world);//dont require the implicit last row
}

void parse_flirt(FILE* input, common_format* source, common_format* target, float world[4][4])
{
	float flirtmat[4][4], srcscale[4][4], trgscale[4][4], sourcemat[4][4], targmat[4][4], temp[4][4], temp2[4][4], temp3[4][4];
	read34(input, flirtmat);//now, do stuff
	int i, j;
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 4; ++j)
		{//flirt doesn't understand origins, so use "scaled coordinates"
			sourcemat[i][j] = source->sform[i][j];
			targmat[i][j] = target->sform[i][j];
			srcscale[i][j] = 0.0f;
			trgscale[i][j] = 0.0f;
		}
		srcscale[3][i] = 0.0f;
		trgscale[3][i] = 0.0f;
		sourcemat[3][i] = 0.0f;
		targmat[3][i] = 0.0f;
	}
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			srcscale[i][i] += sourcemat[j][i] * sourcemat[j][i];
			trgscale[i][i] += targmat[j][i] * targmat[j][i];
		}
		srcscale[i][i] = sqrt(srcscale[i][i]);
		trgscale[i][i] = sqrt(trgscale[i][i]);
	}
	srcscale[3][3] = 1.0f;
	trgscale[3][3] = 1.0f;
	sourcemat[3][3] = 1.0f;
	targmat[3][3] = 1.0f;
	if (det3(sourcemat) > 0.0f)
	{
		srcscale[0][3] = (source->length[0] - 1) * srcscale[0][0];
		srcscale[0][0] = -srcscale[0][0];
	}
	if (det3(targmat) > 0.0f)
	{
		trgscale[0][3] = (target->length[0] - 1) * trgscale[0][0];
		trgscale[0][0] = -trgscale[0][0];
	}
	mult4(flirtmat, srcscale, temp2);
	inv4(trgscale, temp);
	mult4(temp, temp2, temp3);//now in voxel space
	inv4(sourcemat, temp);
	mult4(temp3, temp, temp2);
	mult4(targmat, temp2, world);
}

void parse_xr3d(FILE* input, common_format* source, common_format* target, float world[4][4])
{
	float xr3dmat[4][4], srcmat[4][4], trgmat[4][4], srccent[4][4], trgcent[4][4], temp[4][4], temp2[4][4], temp3[4][4];
	read34(input, xr3dmat);//dont require the implicit last row
	int i, j;
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 4; ++j)
		{//order[3] should always be 3, no worries
			srcmat[i][j] = source->sform[i][j];
			trgmat[i][j] = target->sform[i][j];
			srccent[i][j] = 0.0f;
			trgcent[i][j] = 0.0f;
		}
		srccent[i][i] = fabs(source->sform[source->order[i]][i]);
		srccent[i][3] = -srccent[i][i] * (source->length[i] / 2);//yes, integer division
		trgcent[i][i] = fabs(target->sform[target->order[i]][i]);
		trgcent[i][3] = -trgcent[i][i] * (target->length[i] / 2);
		srcmat[3][i] = 0.0f;
		trgmat[3][i] = 0.0f;
		srccent[3][i] = 0.0f;
		trgcent[3][i] = 0.0f;
	}
	srcmat[3][3] = 1.0f;
	trgmat[3][3] = 1.0f;
	srccent[3][3] = 1.0f;
	trgcent[3][3] = 1.0f;
	inv4(xr3dmat, temp);
	mult4(temp, srccent, temp2);
	inv4(trgcent, temp);
	mult4(temp, temp2, temp3);//now in voxel space
	inv4(srcmat, temp);
	mult4(temp3, temp, temp2);
	mult4(trgmat, temp2, world);
}

void save_t4(FILE* output, common_format* source, common_format* target, float world[4][4])
{
	float invmat[4][4], flat[16];
	inv4(world, invmat);
	int i, j;
	for (i = 0; i < 4; ++i)
	{
		for (j = 0; j < 4; ++j)
		{//order[3] should always be 3, no worries
			flat[i + j * 4] = invmat[source->order[i]][target->order[j]];
		}
	}
	t4_write(output, flat, 0.0f);//0 for scale should skip writing scale
}

void save_world(FILE* output, float world[4][4])
{
	write44(output, world);//write out the implicit last row
}

void save_flirt(FILE* output, common_format* source, common_format* target, float world[4][4])
{
	float flirtmat[4][4], srcscale[4][4], trgscale[4][4], sourcemat[4][4], targmat[4][4], temp[4][4], temp2[4][4], temp3[4][4];
	int i, j;
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 4; ++j)
		{//flirt doesn't understand origins, so use "scaled coordinates"
			sourcemat[i][j] = source->sform[i][j];
			targmat[i][j] = target->sform[i][j];
			srcscale[i][j] = 0.0f;
			trgscale[i][j] = 0.0f;
		}
		srcscale[3][i] = 0.0f;
		trgscale[3][i] = 0.0f;
		sourcemat[3][i] = 0.0f;
		targmat[3][i] = 0.0f;
	}
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			srcscale[i][i] += sourcemat[j][i] * sourcemat[j][i];
			trgscale[i][i] += targmat[j][i] * targmat[j][i];
		}
		srcscale[i][i] = sqrt(srcscale[i][i]);
		trgscale[i][i] = sqrt(trgscale[i][i]);
	}
	srcscale[3][3] = 1.0f;
	trgscale[3][3] = 1.0f;
	sourcemat[3][3] = 1.0f;
	targmat[3][3] = 1.0f;
	if (det3(sourcemat) > 0.0f)
	{
		srcscale[0][3] = (source->length[0] - 1) * srcscale[0][0];
		srcscale[0][0] = -srcscale[0][0];
	}
	if (det3(targmat) > 0.0f)
	{
		trgscale[0][3] = (target->length[0] - 1) * trgscale[0][0];
		trgscale[0][0] = -trgscale[0][0];
	}
	/*inv4(world, temp2);
	mult4(temp2, targmat, temp);
	inv4(sourcemat, temp2);
	mult4(temp2, temp, temp3);//temp3 should be fslvoxmat
	inv4(trgscale, temp);
	mult4(temp3, temp, temp2);
	mult4(srcscale, temp2, temp);
	inv4(temp, flirtmat);//*/
	mult4(world, sourcemat, temp);//some matrix algebra, less inversion, preserves a bit more at float32 precision
	inv4(targmat, temp2);//of course, it gets written with lower precision due to text format, so it doesn't really matter
	mult4(temp2, temp, temp3);//temp3 is in voxel space
	inv4(srcscale, temp);
	mult4(temp3, temp, temp2);
	mult4(trgscale, temp2, flirtmat);
	write44(output, flirtmat);
}

void save_xr3d(FILE* output, common_format* source, common_format* target, float world[4][4])
{
	float xr3dmat[4][4], srcmat[4][4], trgmat[4][4], srccent[4][4], trgcent[4][4], temp[4][4], temp2[4][4], temp3[4][4];
	int i, j;
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 4; ++j)
		{//order[3] should always be 3, no worries
			srcmat[i][j] = source->sform[i][j];
			trgmat[i][j] = target->sform[i][j];
			srccent[i][j] = 0.0f;
			trgcent[i][j] = 0.0f;
		}
		srccent[i][i] = fabs(source->sform[source->order[i]][i]);
		srccent[i][3] = -srccent[i][i] * (source->length[i] / 2);//yes, integer division
		trgcent[i][i] = fabs(target->sform[target->order[i]][i]);
		trgcent[i][3] = -trgcent[i][i] * (target->length[i] / 2);
		srcmat[3][i] = 0.0f;
		trgmat[3][i] = 0.0f;
		srccent[3][i] = 0.0f;
		trgcent[3][i] = 0.0f;
	}
	srcmat[3][3] = 1.0f;
	trgmat[3][3] = 1.0f;
	srccent[3][3] = 1.0f;
	trgcent[3][3] = 1.0f;
	mult4(world, srcmat, temp);
	inv4(trgmat, temp2);
	mult4(temp2, temp, temp3);
	inv4(srccent, temp);
	mult4(temp3, temp, temp2);
	mult4(trgcent, temp2, temp);
	inv4(temp, xr3dmat);
	write44(output, xr3dmat);
}

int main(int argc, char** argv)
{
	printf("%s\n", rcsid);
	if (argc < 6 || strlen(argv[1]) < 2)
	{	//linebreak guide:                                                   72+"\t"= 74|   80|
		printf("usage: %s <conversion string> <original source volume> <original target\n", argv[0]);
		printf("\tvolume> <input affine file> <desired source volume> <desired target\n");
		printf("\tvolume> <output affine file>\n\n");
		printf("\t4dfp volumes are required for t4 and xr3d, other formats (except world\n");
		printf("\tmatrices) require nifti.  For world matrices, volumes are not required,\n");
		printf("\tuse any string in place of the volume arguments.  World matrices are a\n");
		printf("\tsimple coordinate to coordinate transform that follows the nifti\n");
		printf("\tcoordinate convention (mm units, +x=right, +y=anterior, +z=superior),\n");
		printf("\tif you multiply a source coordinate (x,y,z) through the matrix (and add\n");
		printf("\tthe fourth column), you get the corresponding (x',y',z') coordinate in\n");
		printf("\tthe target space.\n\n");
		printf("\tThis utility has not been tested much on images that are not transverse\n");
		printf("\t(axis order left/right, posterior/anterior, inferior/superior), check\n");
		printf("\toutput for sanity.\n\n");
		printf("\tConversion string is a 2 character string, where the first character\n");
		printf("\tgives the current format, and second character gives the desired output\n");
		printf("\tformat.\n\n");
		printf("\tFormat characters are:\n");
		printf("\t4\t_t4 file from 4dfp utilities\n");
		printf("\tf\tflirt matrix from FSL\n");
		printf("\tw\tworld matrix, (x,y,z)->(x',y',z')\n");
		printf("\tx\textracted matrices from 4dfp xr3d.mat files\n\n");
		printf("Example: %s 4f T1.4dfp.ifh $REFDIR/711-2B.4dfp.ifh T1_to_711-2B_t4 T1.nii 711-2B.nii T1_to_711-2B.mat\n\n", argv[0]);
		printf("This example will convert the 'T1_to_711-2B_t4' t4 file, which warps\n");
		printf("'T1.4dfp.ifh' to the 711-2B space, to FSL's flirt matrix format such that it\n");
		printf("will perform the same warp when applied to the file 'T1.nii' with the reference\n");
		printf("volume '711-2B.nii'.\n\n");
		return 0;
	}
	common_format* origsource = NULL, *origtarget = NULL, *newsource = NULL, *newtarget = NULL;
	FILE* input = NULL, *output = NULL;
	char inopen = 1;
	switch (argv[1][0])//open input files
	{
	case 'x':
		input = fopen(argv[4], "r");
		if (!input)
		{
			printf("error: cannot access input affine file\n");
			return -1;
		}
		origsource = parse_4dfp(argv[2], NULL);
		origtarget = parse_4dfp(argv[3], NULL);
		break;
	case '4':
		inopen = 0;
		origsource = parse_4dfp(argv[2], NULL);
		origtarget = parse_4dfp(argv[3], NULL);
		break;
	case 'f':
		input = fopen(argv[4], "r");
		if (!input)
		{
			printf("error: cannot access input affine file\n");
			return -1;
		}
		origsource = parse_nifti(argv[2], 0);
		origtarget = parse_nifti(argv[3], 0);
		break;
	case 'w':
		input = fopen(argv[4], "r");
		break;
	default:
		printf("error: unknown affine format specifier\n");
		return -1;
	};
	if ((argv[1][0] != 'w') && (!origsource || !origtarget))//sanity check volumes
	{
		printf("error: cannot open original volume files\n");
		return -1;
	}
	if (inopen && !input)//sanity check affine
	{
		printf("error: cannot access input affine file\n");
		return -1;
	}
	switch (argv[1][1])//open volumes for output space
	{
	case '4':
	case 'x':
		newsource = parse_4dfp(argv[5], NULL);
		newtarget = parse_4dfp(argv[6], NULL);
		break;
	case 'f':
		newsource = parse_nifti(argv[5], 0);
		newtarget = parse_nifti(argv[6], 0);
		break;
	case 'w':
		break;//world to world doesn't need the volumes
	default:
		printf("error: unknown affine format specifier\n");
		return -1;
	};
	if ((argv[1][1] != 'w') && (!newsource || !newtarget))//sanity check volumes
	{
		printf("error: cannot open original volume files\n");
		return -1;
	}
	if (origsource && origtarget)
	{
		to_lpi(origsource);//merely fills in the order and orientation fields
		to_lpi(origtarget);
	}
	if (newsource && newtarget)
	{
		to_lpi(newsource);//ditto
		to_lpi(newtarget);
	}
	float world[4][4];
	switch (argv[1][0])//parse it
	{
	case '4':
		parse_t4(argv[4], origsource, origtarget, world);
		break;
	case 'w':
		parse_world(input, world);
		break;
	case 'f':
		parse_flirt(input, origsource, origtarget, world);
		break;
	case 'x':
		parse_xr3d(input, origsource, origtarget, world);
		break;
	default:
		printf("unknown format character: '%c'\n", argv[1][0]);
	};
	if (inopen) fclose(input);
	output = fopen(argv[7], "w");
	if (!output)
	{
		printf("error: cannot access output affine file\n");
		return -1;
	}
	switch(argv[1][1])
	{
	case '4':
		save_t4(output, newsource, newtarget, world);
		break;
	case 'w':
		save_world(output, world);
		break;
	case 'f':
		save_flirt(output, newsource, newtarget, world);
		break;
	case 'x':
		save_xr3d(output, newsource, newtarget, world);
		break;
	default:
		printf("unknown format character: '%c'\n", argv[1][1]);
	};
	free_common_format(origsource);
	free_common_format(origtarget);
	free_common_format(newsource);
	free_common_format(newtarget);
	fclose(output);
	return 0;
}
