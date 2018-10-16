#include "parse_common.h"
#include "nifti-format.h"
#include "4dfp-format.h"

common_format* parse_common_format(char* root, char rec_open)
{
	common_format* ret = parse_nifti(root, rec_open);
	if (ret == NULL) ret = parse_4dfp(root, NULL);//second argument is t4 file
	return ret;
}

