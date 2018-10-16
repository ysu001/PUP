#$Header: /data/petsun4/data1/src_solaris/gauss_4dfp/RCS/hsphere_4dfp.mak,v 1.3 2007/04/17 05:29:01 avi Exp $
#$Log: hsphere_4dfp.mak,v $
# Revision 1.3  2007/04/17  05:29:01  avi
# gcc compliant
# fftsun.o -> fftsol.o
#
# Revision 1.2  2006/09/25  17:42:38  avi
# ${PROG} ${RELEASE}
#
PROG	= hsphere_4dfp
CSRCS	= ${PROG}.c
FSRCS	= hsphere3d.f
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${RMS}/fftsol.o ${RMS}/imgpad.o

CFLAGS	= -O -I${TRX} -I${RMS}
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
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}
