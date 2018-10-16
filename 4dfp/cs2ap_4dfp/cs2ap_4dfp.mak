#$Header: /data/petsun4/data1/src_solaris/cs2ap_4dfp/RCS/cs2ap_4dfp.mak,v 1.1 2010/08/29 02:57:18 avi Exp $
#$Log: cs2ap_4dfp.mak,v $
# Revision 1.1  2010/08/29  02:57:18  avi
# Initial revision
#

PROG	= cs2ap_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
BLR	= ${NILSRC}/imgblur_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${BLR}/fimgblur.o

CFLAGS	= -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore
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

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS}

checkin:
	ci ${CSRCS}

release: ${PROG}
	chmod 751 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
