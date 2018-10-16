#$Header: /data/petsun4/data1/src_solaris/nonlin/RCS/morphc_4dfp.mak,v 1.4 2007/07/04 03:25:55 avi Exp $
#$Log: morphc_4dfp.mak,v $
# Revision 1.4  2007/07/04  03:25:55  avi
# gcc compliant
#
# Revision 1.3  2007/07/04  02:32:18  avi
# eliminate -lrms
#
# Revision 1.2  2007/07/03  05:13:45  avi
# endian compliant i/o
#
# Revision 1.1  2007/07/03  03:54:10  avi
# Initial revision
#
PROG	= morphc_4dfp
CSRCS	= ${PROG}.c piic.c iic.c lud.c
FSRCS	= fvfield.f util.f b2v4m.f inc.f imgval0.f cmplxfft3d.f
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
PAR     = ${NILSRC}/parzen_4dfp
RMS	= ${NILSRC}/librms
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${PAR}/parzen.o ${PAR}/imgpac.o ${PAR}/mskpac.o \
	  ${RMS}/fftsol.o ${RMS}/imgpad.o ${RMS}/matopr.o

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
	/bin/rm ${OBJS} ${PROG}
