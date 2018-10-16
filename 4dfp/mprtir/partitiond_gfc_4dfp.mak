#$Header: /data/petsun4/data1/src_solaris/mprtir/RCS/partitiond_gfc_4dfp.mak,v 1.2 2008/01/02 01:21:16 avi Exp $
#$Log: partitiond_gfc_4dfp.mak,v $
# Revision 1.2  2008/01/02  01:21:16  avi
# linux gcc v4 compatible
#
# Revision 1.1  2007/10/15  21:52:25  avi
# Initial revision
#
PROG	= partitiond_gfc_4dfp
CSRCS	= ${PROG}.c
FSRCS	= fitgain3d.f powsub.f fnegdef.f fitgain3dd.f
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
BLUR	= ${NILSRC}/imgblur_4dfp
LOBJS   = ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${RMS}/matopr.o ${RMS}/eigen.o ${BLUR}/fimgblur.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= -lm

CFLAGS	= -O -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore -fcray-pointer
	LIBS	= -lm -lgfortran
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

${PROG}: $(OBJS)
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<
clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

release: ${PROG}
	chmod 711 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
