/*$Header: /data/petsun4/data1/src_solaris/ecat/RCS/ecat_header.h,v 1.1 2010/09/01 17:36:51 larsc Exp larsc $*/
/*$Log: ecat_header.h,v $
 * Revision 1.1  2010/09/01  17:36:51  larsc
 * Initial revision
 **/
#include "matrix.h"

void show_main_header	(Main_header *mh);
void show_subheader	(MatrixFile *mptr, MatrixData *matrix);
void show_scan_subheader	(Scan_subheader   *sh);
void show_Scan3D_subheader	(Scan3D_subheader *sh);
void show_image_subheader	(MatrixData       *matrix);
void show_attn_subheader	(Attn_subheader   *ah);
void show_norm_subheader	(Norm_subheader   *nh);
void show_norm3d_subheader	(Norm3D_subheader *nh);
static char* storage_order (int idx);
void day_time (char *str, int day, int month, int year, int hour, int minute, int sec);
