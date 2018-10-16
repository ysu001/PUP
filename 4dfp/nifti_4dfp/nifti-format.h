/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/nifti-format.h,v 1.4 2010/06/18 21:00:59 coalsont Exp $*/
/*$Log: nifti-format.h,v $
 * Revision 1.4  2010/06/18  21:00:59  coalsont
 * added define checks to avoid __attribute__ on non-gcc compilers
 *
 * Revision 1.3  2009/08/17  22:57:44  coalsont
 * fopen_nifti
 *
 * Revision 1.2  2009/08/07  20:42:04  coalsont
 * endian control
 *
 * Revision 1.1  2009/07/08  00:54:37  coalsont
 * Initial revision
 **/
/* Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009 */
/* Kevin P. Barry [ta0kira@users.berlios.de], BMCLAB, 22 May 2009 */

#ifndef nifti_format_h
#define nifti_format_h

#include "nifti1.h"
#include "common-format.h"
#include "stdio.h"

typedef struct
{
	struct nifti_1_header  header;
	struct nifti1_extender extender;
#ifndef NOT_GCC
} __attribute__ ((packed)) header_block;
#else
} header_block;
#endif

common_format *parse_nifti(const char* fFile, char rec_open);//ONLY for reading, currently
int save_nifti(common_format* iImage, const char* fFile, char* program_name, char control);
FILE* fopen_nifti(char* fFile, int* dim, float* spacing, float* center, char* program_name, char control);//ONLY for writing
FILE* fopen_nifti_sform(char* fFile, int* dim, float sform[3][4], float timestep, char* program_name, char control);//ditto

#endif /*nifti_format_h*/
