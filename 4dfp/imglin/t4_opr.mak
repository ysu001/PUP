#$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4_opr.mak,v 1.1 2010/02/19 02:06:17 avi Exp $
#$Log: t4_opr.mak,v $
# Revision 1.1  2010/02/19  02:06:17  avi
# Initial revision
#

PROG	= t4_opr
CSRCS	= ${PROG}.c t4_io.c
FSRCS	= 
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
