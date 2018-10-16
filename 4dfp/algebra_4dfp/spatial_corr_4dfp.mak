#$Header: /data/petsun4/data1/src_solaris/algebra_4dfp/RCS/spatial_corr_4dfp.mak,v 1.1 2010/03/17 20:31:37 avi Exp $
#$Log: spatial_corr_4dfp.mak,v $
# Revision 1.1  2010/03/17  20:31:37  avi
# Initial revision
#

PROG	= spatial_corr_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

CFLAGS	= -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc -std=c99 ${CFLAGS}	# -std=c99 enables use of isnormal()
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm -lsunmath
endif

.c.o:
	${CC} -c $<
.f.o:
	${FC} -c $<

${PROG}: spatial_corr_4dfp.o
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
