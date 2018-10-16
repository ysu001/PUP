#$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/mat_algebra.mak,v 1.2 2010/02/15 05:36:33 avi Exp $
#$Log: mat_algebra.mak,v $
# Revision 1.2  2010/02/15  05:36:33  avi
# include dmatinv.o
#
# Revision 1.1  2009/03/05  03:48:35  avi
# Initial revision
#
PROG	= mat_algebra
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
	CC	= gcc -std=c99 ${CFLAGS}	# -std=c99 enables use of isnormal()
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
