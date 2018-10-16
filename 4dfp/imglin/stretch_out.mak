#$Header: /data/petsun4/data1/src_solaris/imglin/RCS/stretch_out.mak,v 1.1 2010/03/03 01:31:41 avi Exp $
#$Log: stretch_out.mak,v $
# Revision 1.1  2010/03/03  01:31:41  avi
# Initial revision
#

PROG	= stretch_out
CSRCS	= ${PROG}.c t4_io.c
FSRCS	= stretchout.f
OBJS 	= ${FSRCS:.f=.o} ${CSRCS:.c=.o}
RMS	= ${NILSRC}/librms
LOBJS	= ${RMS}/matopr.o ${RMS}/eigen.o

CFLAGS	= -O -I. -I${RMS}
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
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS}

checkout:
	co ${CSRCS} ${FSRCS} 

checkin:
	ci ${CSRCS} ${FSRCS} 
