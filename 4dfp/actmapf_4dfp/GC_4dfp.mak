#$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/GC_4dfp.mak,v 1.2 2008/09/14 02:50:43 avi Exp $
#$Log: GC_4dfp.mak,v $
# Revision 1.2  2008/09/14  02:50:43  avi
# linux compliant
#
# Revision 1.1  2007/07/08  03:50:35  avi
# Initial revision
#

PROG	= GC_4dfp
CSRCS	= ${PROG}.c expandf.c conc.c
FSRCS	= fGC.f
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o \
	  ${RMS}/matopr.o

OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= -lm

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -I. -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm -lgfortran
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

${PROG}: ${OBJS} 
	${FC} -o ${PROG} ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
