/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/transform.h,v 1.4 2010/08/16 18:20:29 coalsont Exp $*/
/*$Log: transform.h,v $
 * Revision 1.4  2010/08/16  18:20:29  coalsont
 * cast to 64bit integer in get_voxel
 *
 * Revision 1.3  2010/06/18  21:00:59  coalsont
 * added functions to .h so they could be used elsewhere
 *
 * Revision 1.2  2009/08/07  20:42:04  coalsont
 * restructuring, endian control
 *
 * Revision 1.1  2009/07/08  00:54:37  coalsont
 * Initial revision
 **/
/* Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009 */
/* Kevin P. Barry [ta0kira@users.berlios.de], BMCLAB, 22 May 2009 */

#ifndef transform_h
#define transform_h

#include "common-format.h"


int local_endian();
void flip_endian(void*, int width, int number);
void to_lpi(common_format* iImage);
void auto_orient_header(common_format*);
float get_value_float(const void* data, common_format* iImage, int rReverse);
void auto_orient(void*, common_format*, char);
void invert_affine(float affineIn[3][4], float affineOut[3][4]);
int map_voxel(void *lLeft, const void *rRight, common_format *iImage, int rReverse, char reswap);

void *get_voxel(void *dData, const vindex pPos[4], const vindex length[4], uint8_t bperpix);

static inline int compare_endian(int eEndian)
{ return local_endian() ^ (eEndian == 0); }


#endif /*transform_h*/
