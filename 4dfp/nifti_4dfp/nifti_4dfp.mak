#$Log: nifti_4dfp.mak,v $
# Revision 1.13  2018/10/16  09:50:55  ysu001
# moved "-lm" from LINK to release target
#
# Revision 1.12  2012/09/08  01:14:55  avi
# eliiminate -sdt=c99
#
# Revision 1.11  2011/11/23  21:03:09  coalsont
# added missing source file
#
# Revision 1.10  2009/08/07  20:19:22  coalsont
# link flags for optimizing
#
# Revision 1.7  2009/07/29  20:55:12  coalsont
# mystdint.h added to HEADERS
#
# Revision 1.6  2009/07/29  20:50:57  coalsont
# conditional for linux, now needs gmake on solaris
#
# Revision 1.5  2009/07/22  01:24:01  coalsont
# conditional removed, make incompatibilities
# previous version: also changed to optimized flags for gcc
#
# Revision 1.4  2009/07/22  01:01:46  coalsont
# conditional to define HAVE_STDINT
#
# Revision 1.3  2009/07/08  01:00:41  coalsont
# added csrcs to checkin/checkout targets
#
# Tim Coalson [tsc5yc@mst.edu], NRG, 29 June 2009
# Kevin P. Barry [ta0kira@users.berlios.de], BMCLAB, 22 May 2009

PROG	= nifti_4dfp
TRX	= ${NILSRC}/TRX
IMGLIN	= ${NILSRC}/imglin
LOBJS	= ${TRX}/endianio.o ${TRX}/Getifh.o ${TRX}/rec.o ${IMGLIN}/t4_io.o
SUFFIXES= .o .c
ifeq ($(shell uname), SunOS)
	COMPILE = -W -Wall -O2
else
	COMPILE = -W -Wall -O2 -DHAVE_STDINT
endif
INCLUDE	= -isystem -I. -I${TRX} -I${IMGLIN}
LINK	= -O2
HEADERS	= common-format.h 4dfp-format.h nifti-format.h transform.h nifti1.h mystdint.h parse_common.h
IHEADERS= ${TRX}/ANALYZE.h ${TRX}/endianio.h ${TRX}/Getifh.h ${TRX}/ifh.h ${TRX}/rec.h ${IMGLIN}t4_io.h
CSRCS	= nifti_4dfp.c 4dfp-format.c nifti-format.c split.c transform.c common-format.c parse_common.c
OBJECTS	= nifti_4dfp.o 4dfp-format.o nifti-format.o split.o transform.o common-format.o parse_common.o

.PHONY: all clean

.SUFFIXES: ${SUFFIXES}

all: $(PROG)

.c.o: $(HEADERS) ${IHEADERS} nifti_4dfp.mak
	gcc $(COMPILE) $(INCLUDE) -c $< -o $@

$(PROG): $(OBJECTS) nifti_4dfp.mak
	gcc $(LINK) $(OBJECTS) ${LOBJS} -o $(PROG) -lm

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm -f $(PROG) $(OBJECTS)

checkout:
	co ${CSRCS} ${HEADERS}

checkin:
	ci ${CSRCS} ${HEADERS}

