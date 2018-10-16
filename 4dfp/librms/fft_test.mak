#$Header: /data/petsun4/data1/src_solaris/librms/RCS/fft_test.mak,v 1.1 2007/03/27 05:00:56 avi Exp avi $
#$Log: fft_test.mak,v $
# Revision 1.1  2007/03/27  05:00:56  avi
# Initial revision
#

PROG	= fft_test
CSRCS	= ${PROG}.c
FSRCS	= fftsol.f npad.f
LIN	= ${NILSRC}/imglin
LOBJS	= ${LIN}/dnormal.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

.c.o:
	${CC} -c $<
.f.o:
	${FC} -c $<

CFLAGS	= -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore
	LIBS	= -lm -lgfortran
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

${PROG}: ${OBJS}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

clean:
	/bin/rm ${OBJS} ${PROG}

test:
ifeq (${OSTYPE}, linux)
	echo "linux"
else
	echo "not linux"
endif
