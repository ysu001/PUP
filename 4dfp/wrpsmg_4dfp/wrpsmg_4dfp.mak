#$Header: /data/petsun4/data1/src_solaris/wrpsmg_4dfp/RCS/wrpsmg_4dfp.mak,v 1.3 2010/01/15 01:57:04 avi Exp $
#$Log: wrpsmg_4dfp.mak,v $
# Revision 1.3  2010/01/15  01:57:04  avi
# linux compliant
#
# Revision 1.2  2010/01/15  01:37:07  avi
# Solaris 10 compliant
#
# Revision 1.1  2006/03/06  05:28:03  avi
# Initial revision
#

PROG	= wrpsmg_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
LIN	= ${NILSRC}/imglin
FLP	= ${NILSRC}/flip_4dfp
REG	= ${NILSRC}/imgreg_4dfp
T4I	= ${NILSRC}/t4imgs_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${FLP}/cflip.o ${REG}/imgvalx.o \
	  ${LIN}/t4_sub.o ${T4I}/ft4imgo.o ${T4I}/set_rnan.o ${LIN}/t4scale.o \
	  ${RMS}/polfit.o ${RMS}/matopr.o ${RMS}/param6opr.o

CFLAGS	= -I${TRX} -I${LIN} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore
	Q	= $(wildcard /usr/lib*/libgfortran.so.1)
	ifeq ($(Q), "")
		LIBS	= -lm -lg2c
	else
		LIBS	= -lm -lgfortran
	endif
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
