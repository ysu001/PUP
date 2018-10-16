#$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/LDA.mak,v 1.2 2011/08/29 19:40:49 avi Exp $
#$Log: LDA.mak,v $
# Revision 1.2  2011/08/29  19:40:49  avi
# remoave -std=c99 to anable use of M_PI on linux
#
# Revision 1.1  2011/08/22  03:36:13  avi
# Initial revision
#
PROG	= LDA
CSRCS	= ${PROG}.c
RMS	= ${NILSRC}/librms
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LOBJS	= ${RMS}/deigen.o ${RMS}/dgeigen.o ${RMS}/dmatopr.o ${RMS}/dmatinv.o

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -I${RMS} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}		# -std=c99 preclude use of M_PI
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm -lgfortran
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm -lsunmath
endif

${PROG}: ${OBJS} 
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
