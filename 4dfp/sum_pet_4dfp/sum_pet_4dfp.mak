#$Id:  $

PROG	= sum_pet_4dfp
CSRCS	= ${PROG}.c

TRX	= ${NILSRC}/TRX
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o 

CFLAGS	= -I${TRX}
ifeq (${OSTYPE}, linux)
CC	= gcc ${CFLAGS}
else
CC	= cc ${CFLAGS}
endif
LIBS	= -lm

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}
