#$Header: /data/petsun4/data1/src_solaris/algebra_4dfp/RCS/RFX_4dfp.mak,v 1.1 2010/12/31 07:26:24 avi Exp $
#$Log: RFX_4dfp.mak,v $
# Revision 1.1  2010/12/31  07:26:24  avi
# Initial revision
#

PROG	= RFX_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

CFLAGS	= -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc -std=c99 ${CFLAGS}	# -std=c99 enables use of isnormal()
	LIBS	= -lm
else
	CC	= cc ${CFLAGS}
	LIBS	= -lm -lsunmath
endif

.c.o:
	${CC} -c $<

${PROG}: ${PROG}.o 
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
