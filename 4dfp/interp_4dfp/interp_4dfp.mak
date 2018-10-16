#$Header: /data/petsun4/data1/src_solaris/interp_4dfp/RCS/interp_4dfp.mak,v 1.4 2007/11/20 07:07:25 avi Exp $
#$Log: interp_4dfp.mak,v $
# Revision 1.4  2007/11/20  07:07:25  avi
# linux compatible (fftsun.o -> fftsol.o)
#
# Revision 1.3  2007/05/03  17:52:12  avi
# Solaris 10; $NILSRC $RELEASE; eliminate -lrms
#
# Revision 1.2  2005/01/29  05:51:27  avi

PROG	= interp_4dfp
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
CSRCS	= ${PROG}.c
FSRCS	= spline.f
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${RMS}/imgpad.o ${RMS}/fftsol.o
LIBS	= -lm

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -O -I${RMS} -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm -lgfortran
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

${PROG}: $(OBJS)
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

clean:
	rm ${OBJS}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

release: ${PROG}
	chmod 751 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
