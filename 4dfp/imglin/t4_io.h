/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4_io.h,v 1.2 2010/03/03 01:25:34 avi Exp $*/
/*$Log: t4_io.h,v $
 * Revision 1.2  2010/03/03  01:25:34  avi
 * t4_write returns void
 *
 * Revision 1.1  2009/02/23  06:14:49  avi
 * Initial revision
 **/

int	t4_read1 (FILE *fp, float *t4);
void 	t4_init (float *t4);
int	t4_read (char *t4file, float *t4);
int	iget_t4_scale (char *t4file, float *scale);
void	t4_list (FILE *fp, float *t4);
void	t4_write (FILE *fp, float *t4, float scale);
