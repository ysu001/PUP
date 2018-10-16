/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4_sub.h,v 1.1 2009/09/09 05:53:40 avi Exp $*/
/*$Log: t4_sub.h,v $
 * Revision 1.1  2009/09/09  05:53:40  avi
 * Initial revision
 **/

void 	t4_init_	(float *t4);
int	t4_read_	(char *t4file, float *t4);
int	iget_t4_scale_	(char *t4file, float *scale);
void	t4_list_	(float *t4);
int	t4_write_	(float *t4, char *t4file);
