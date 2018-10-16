#$Header: /data/petsun4/data1/src_solaris/mprtir/RCS/2Dhist.mak,v 1.2 2007/08/07 00:44:29 avi Exp $
#$Log: 2Dhist.mak,v $
# Revision 1.2  2007/08/07  00:44:29  avi
# Solaris 10 and linux compliant
#

PROG	= 2Dhist
CSRCS	= ${PROG}.c
FSRCS	=
TRX	= ${NILSRC}/TRX
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LOBJS	= ${TRX}/rec.o ${TRX}/endianio.o ${TRX}/Getifh.o

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
