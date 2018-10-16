#$Header: /data/petsun4/data1/src_solaris/imato4dfp2/RCS/unpack_4dfp.mak,v 1.4 2008/07/24 01:03:32 avi Exp $
#$Log: unpack_4dfp.mak,v $
# Revision 1.4  2008/07/24  01:03:32  avi
# add ${FLP}/cflip.o to ${LOBJS}
#
# Revision 1.3  2007/09/19  02:38:22  avi
# linux compliant
#
# Revision 1.2  2007/05/03  19:55:17  avi
# Solaris 10; endian compliant
# eliminate -lrms
#
# Revision 1.1  2004/11/15  21:05:32  rsachs
# Initial revision
#

PROG	= unpack_4dfp
CSRCS	= ${PROG}.c
FSRCS	= unpack.f
TRX	= ${NILSRC}/TRX
FLP	= ${NILSRC}/flip_4dfp
OBJS	= ${FSRCS:.f=.o} ${CSRCS:.c=.o}
LOBJS	= ${TRX}/Getifh.o ${TRX}/endianio.o ${TRX}/rec.o ${FLP}/cflip.o

CFLAGS	= -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm -lgfortran
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

${PROG}: ${OBJS}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

