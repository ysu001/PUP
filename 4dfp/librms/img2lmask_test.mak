#$Header$
#$Log$

PROG	= img2lmask_test
CSRCS	= ${PROG}.c img2lmask.c
FSRCS	= gauss3d.f
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
LOBJS	= ${TRX}/Getifh.o ${TRX}/endianio.o ${RMS}/imgpad.o ${RMS}/fftsol.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

CFLAGS	= -O -I${TRX} -I${RMS}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -lcray-pointer
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

${PROG}: $(OBJS)
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

release: ${PROG}
	chmod 711 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
