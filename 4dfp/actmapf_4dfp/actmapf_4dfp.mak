# $Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/actmapf_4dfp.mak,v 1.9 2007/11/20 02:46:57 avi Exp $
# $Log: actmapf_4dfp.mak,v $
# Revision 1.9  2007/11/20  02:46:57  avi
# linux compliant
#
# Revision 1.8  2006/09/23  23:03:04  avi
# ${NILSRC} ${RELEASE}
#
# Revision 1.7  2006/09/21  19:56:59  avi
# use endian invariant subroutines
#
# Revision 1.6  2005/07/22  05:08:56  avi
# add conc.c and eliminate all references to libmri and librms
#
# Revision 1.5  2004/05/23  04:05:03  avi
# link with fimg_mode.o
#
# Revision 1.4  2004/05/16  01:55:57  avi
# chmod to program on release
#
# Revision 1.3  1999/03/11  09:17:43  avi
# Getifh.o
#
# Revision 1.2  1999/01/04  08:21:01  avi
# new release
#
# Revision 1.1  1998/10/08  23:51:41  avi
# Initial revision
#
PROG	= actmapf_4dfp
CSRCS	= ${PROG}.c expandf.c conc.c
OBJS	= ${CSRCS:.c=.o}
TRX	= ${NILSRC}/TRX
HST	= ${NILSRC}/img_hist_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${HST}/fimg_mode.o

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -I. -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm -lgfortran
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

${PROG}: ${OBJS} 
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
