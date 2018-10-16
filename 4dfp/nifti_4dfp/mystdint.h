/*$Header: /data/petsun4/data1/src_solaris/nifti_4dfp/RCS/mystdint.h,v 1.2 2010/07/02 15:52:12 coalsont Exp $*/
/*$Log: mystdint.h,v $
 * Revision 1.2  2010/07/02  15:52:12  coalsont
 * updated to not redefine types, removes warning with fortran link
 *
 * Revision 1.1  2009/07/29  20:53:23  coalsont
 * Initial revision
 **/
/* Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009 */
/*
 *  mystdint.h
 *  
 *
 *  Created by NRG on 6/12/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MYSTDINT_H
#define MYSTDINT_H
	#ifdef HAVE_STDINT
		#include "stdint.h"
	#else
	#ifndef int8_t
		typedef char int8_t;
	#endif
	#ifndef int16_t
		typedef short int16_t;
	#endif
	#ifndef int32_t
		typedef int int32_t;
	#endif
	#ifndef int64_t
		typedef long long int64_t;
	#endif
	#ifndef uint8_t
		typedef unsigned char uint8_t;
	#endif
	#ifndef uint16_t
		typedef unsigned short uint16_t;
	#endif
	#ifndef uint32_t
		typedef unsigned int uint32_t;
	#endif
	#ifndef uint64_t
		typedef unsigned long long uint64_t;
	#endif
	#endif
#endif
