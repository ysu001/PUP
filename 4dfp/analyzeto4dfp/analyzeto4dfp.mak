#$Header: /data/petsun4/data1/src_solaris/analyzeto4dfp/RCS/analyzeto4dfp.mak,v 1.6 2007/05/01 03:12:04 avi Exp $
#$Log: analyzeto4dfp.mak,v $
# Revision 1.6  2007/05/01  03:12:04  avi
# linux gcc compliant
#
# Revision 1.5  2006/10/04  05:02:45  avi
# ${PROG} ${RELEASE}
#
# Revision 1.4  2004/09/09  21:48:21  avi
# eliminate lmri and all FORTRAN dependencies
#
# Revision 1.3  2004/09/01  06:00:33  avi
# replace FORTRAN fflip.o with cflip.o
#

PROG	= analyzeto4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
FLP	= ${NILSRC}/flip_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/endianio.o ${TRX}/Getifh.o ${FLP}/cflip.o

CFLAGS	= -O -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc ${CFLAGS}
endif
LIBS	= -lm

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} -lm

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS}

checkin:
	ci ${CSRCS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
