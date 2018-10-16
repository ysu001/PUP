#$Header: /data/petsun4/data1/src_solaris/crop_4dfp/RCS/zero_slice_4dfp.mak,v 1.2 2010/08/20 04:46:32 avi Exp $Log: zero_slice_4dfp.mak,v $
# Revision 1.1  2010/08/20  04:44:12  avi
# Initial revision
#
# Revision 1.2  2006/03/13  07:17:58  avi

PROG	= zero_slice_4dfp
CSRCS	= ${PROG}.c
FSRCS	= fzero_slice.f
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
FLP	= ${NILSRC}/flip_4dfp
LOBJS   = ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${FLP}/cflip.o

CFLAGS	= -I${TRX} -I${FLP} -O
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

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

release: ${PROG}
	chmod 751 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
