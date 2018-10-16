#$Header: /data/petsun4/data1/src_solaris/peak_4dfp/RCS/read_4dfp.mak,v 1.2 2007/09/23 03:01:50 avi Exp ${NILSRC}/peak_4dfp/RCS/read_4dfp.mak,v 1.1 2006/03/03 03:10:06 avi Exp avi $
#$Log: read_4dfp.mak,v $
# Revision 1.2  2007/09/23  03:01:50  avi
# Solaris10 and linux compliant
#
# Revision 1.1  2006/03/03  03:10:06  avi
# Initial revision
#

PROG	= read_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
RMS	= ${NILSRC}/librms
T4I	= ${NILSRC}/t4imgs_4dfp
REG	= ${NILSRC}/imgreg_4dfp
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${REG}/imgvalx.o ${T4I}/cvrtflip.o \
	  ${RMS}/polfit.o ${RMS}/matopr.o

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
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
clean:
	/bin/rm ${OBJS} ${PROG}
