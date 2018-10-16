/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/nifti_4dfp.c,v 1.9 2011/09/13 03:40:37 avi Exp $*/
/*$Log: nifti_4dfp.c,v $
 * Revision 1.9  2011/09/13  03:40:37  avi
 * option -N (optionally suppress saving mmppic and center fields in output 4dfp ifh)
 * option -@ to preserve compatibility with other 4dfp tools (@[b|l] still works)
 * neater usage
 *
 * Revision 1.8  2010/07/08  15:27:21  coalsont
 * used getroot to strip extensions
 *
 * Revision 1.7  2010/06/18  21:00:59  coalsont
 * added ifh2hdr system call on -4
 *
 * Revision 1.6  2009/08/20  20:14:21  coalsont
 * print rcsid
 *
 * Revision 1.5  2009/08/17  22:57:44  coalsont
 * superficial
 *
 * Revision 1.4  2009/08/07  20:42:04  coalsont
 * added endian control for output
 *
 * Revision 1.3  2009/07/22  01:01:14  coalsont
 * fixed uninitialized warnings
 *
 * Revision 1.2  2009/07/08  01:01:35  coalsont
 * program for 4dfp/nifti interconversion
 **/
static char* rcsid = "$Id: nifti_4dfp.c,v 1.9 2011/09/13 03:40:37 avi Exp $";
/* Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009 */
/* Kevin P. Barry [ta0kira@users.berlios.de], BMCLAB, 22 May 2009 */

#include <stdio.h>
#include <string.h>

#include "common-format.h"
#include "4dfp-format.h"
#include "nifti-format.h"
#include "transform.h"
#include "rec.h"
#include "endianio.h"

#define MAXL 256

#define HEADER_EXT	".4dfp.ifh"
#define IMAGE_EXT	".4dfp.img"

static void get_image_name(char *fFile, char **iImage)
{/* used to get the base name for the .rec file */
	size_t original_length = strlen(fFile);
	/*int extension_index = original_length;
	
	if (original_length > strlen(HEADER_EXT) &&
		strcmp(HEADER_EXT, fFile + original_length - strlen(HEADER_EXT)) == 0)
		extension_index = original_length - strlen(HEADER_EXT);//*/
	char* temp = malloc(original_length + 1);
	getroot(fFile, temp);
	int extension_index = strlen(temp);
	free(temp);
	
	*iImage = calloc(1, extension_index + strlen(IMAGE_EXT) + 1);
	strncpy(*iImage, fFile, extension_index);
	strncpy(*iImage + extension_index, IMAGE_EXT, strlen(IMAGE_EXT));
}

int main(int argc, char *argv[])
{
	printf("%s\n", rcsid);
	char *program_name = NULL;
	char haveDir = 0, toNii = 0, haveT4 = 0, control = '\0', save_center = 1;
	char* inname = NULL, *outname = NULL, *t4name = NULL, *ptr, c, *imgname;
	int k, i;
	program_name = argv[0];
	
	/******************************/
	/* get command line arguments */
	/******************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			ptr = argv[i];
			while ((c = *(++ptr))) switch (c) {
				case '4': haveDir++;	toNii = 0;		break;
				case 'N': save_center = 0;			break;
				case 'n': haveDir++;	toNii = 1;		break;
				case 'T': haveT4++; t4name = argv[i + 1]; ++i;	break;
				case '@': control = *(++ptr);			break;
			}
		} else if (*argv[i] == '@') {
			control = argv[i][1];
		} else switch (k) {
		 	case 0: inname  = argv[i];	k++; break;
		 	case 1: outname = argv[i];	k++; break;
			default: k++;
		}	
	}
	if (k != 2 || haveDir != 1 || haveT4 > 1) {
		printf ("Usage: %s <-4 or -n> <infile> <outfile> [options]\n", program_name);
		printf (" e.g.: %s -n time_BOXzstat_333_t88.4dfp.ifh time_BOXzstat_333_t88.nii\n", program_name);
		printf ("\toptions\n");
		printf ("\t-T <t4 file>\tspecify a t4 file to use converting TO NIfTI from 4dfp\n");
		printf ("\t-n\tconvert TO NIfTI from 4dfp\n");
		printf ("\t-4\tconvert TO 4dfp from NIfTI\n");
		printf ("\t-N\tsuppress saving of mmppix and center fields in output ifh\n");
		printf ("\t-@<val>\tspecify endianness for output, b or B for big, l or L for little\n");
		printf ("N.B.:\texactly one of -4 or -n must be specified\n");
		printf ("N.B.:\t\".4dfp.ifh\" or \".nii\" are appended to filenames specified without extension\n");
		printf ("N.B.:\toption -N has effect only on converting nii->4dfp\n");
		printf ("N.B.:\toption -T has effect only on converting 4dfp->nii\n");
		return 1;
	}
	common_format *neutral_image = NULL;
	int save_status = 0;
	
	if (toNii)
	{/*convert from 4dfp to nifti*/
		if ((neutral_image = parse_4dfp(inname, t4name)))
			save_status = save_nifti(neutral_image, outname, program_name, control);
	}	
	else
	{/*convert from nifti to 4dfp*/
		get_image_name(outname, &imgname);
		startrece(imgname, argc, argv, rcsid, '\0');/* let startrece figure out local endian, even though we have a function for it */
		if ((neutral_image = parse_nifti(inname, 1)))
			save_status = save_4dfp(neutral_image, outname, argc, argv, control, save_center);
		free(imgname);
		endrec();
#ifndef STANDALONE
		char command[256];
		sprintf(command, "ifh2hdr %s", outname);
		system(command);
#endif
	}

	if (!neutral_image || !save_status)
	{
		fprintf(stderr, neutral_image ?
				"%s: error saving file '%s'\n" : "%s: error parsing file '%s'\n",
				program_name, neutral_image ? outname : inname);
		free_common_format(neutral_image);
		return 1;
	}
	
	free_common_format(neutral_image);
	return 0;
}
