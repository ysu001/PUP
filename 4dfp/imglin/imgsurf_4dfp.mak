#$Header: /data/petsun4/data1/src_solaris/imglin/RCS/imgsurf_4dfp.mak,v 1.2 2009/02/25 21:29:26 avi Exp $
#$Log: imgsurf_4dfp.mak,v $
# Revision 1.2  2009/02/25  21:29:26  avi
# accommodate 64bit architecture
#
# Revision 1.1  2007/09/05  04:05:03  avi
# Initial revision
#

PROG	= imgsurf_4dfp
CSRCS	= ${PROG}.c
FSRCS	= fimgsurf.f grad_opt.f
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
T4I	= ${NILSRC}/t4imgs_4dfp
REG	= ${NILSRC}/imgreg_4dfp
LOBJS	= ${TRX}/Getifh.o ${TRX}/endianio.o ${REG}/imgvalx.o \
	  ${RMS}/param6opr.o ${RMS}/polfit.o ${RMS}/matopr.o ${T4I}/ft4imgo.o ${T4I}/set_rnan.o

.c.o:
	${CC} -c $<
.f.o:
	${FC} -c $<

CFLAGS	= -I${TRX} -O
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

${PROG}: ${OBJS}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} 

checkin:
	ci ${CSRCS} 
