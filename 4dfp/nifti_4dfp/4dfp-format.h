/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/4dfp-format.h,v 1.4 2011/09/13 03:43:25 avi Exp $*/
/*$Log: 4dfp-format.h,v $
 * Revision 1.4  2011/09/13  03:43:25  avi
 * save_4dfp() now has additional argument (save_center)
 *
 * Revision 1.3  2010/07/08  15:27:21  coalsont
 * removed const to allow getroot without warnings
 *
 * Revision 1.2  2009/08/07  20:42:04  coalsont
 * endian control
 *
 * Revision 1.1  2009/07/08  00:54:37  coalsont
 * Initial revision
 **/
/* Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009 */
/* Kevin P. Barry [ta0kira@users.berlios.de], BMCLAB, 22 May 2009 */

#ifndef fdfp_format_h
#define fdfp_format_h

#include "common-format.h"

common_format *parse_4dfp(char*, char* t4name);
int save_4dfp(common_format*, char*, int argc, char* argv[], char control, int save_center);

#endif /*fdfp_format_h*/
