#$Header: /data/petsun4/data1/src_solaris/mpetto4dfp/RCS/mpetto4dfp.mak,v 1.2 2007/04/28 00:34:50 avi Exp $
#$Log: mpetto4dfp.mak,v $
# Revision 1.2  2007/04/28  00:34:50  avi
# endian compliant
#
# Revision 1.1  2006/03/06  08:11:15  avi
# Initial revision
#

PROG	= mpetto4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
FLP	= ${NILSRC}/flip_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${FLP}/cflip.o

CFLAGS	= -O -I${TRX}
FC	= f77 -O -e -I4
CC	= cc ${CFLAGS}

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
