#$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/covariance.mak,v 1.4 2011/02/21 00:35:57 avi Exp $
#$Log: covariance.mak,v $
# Revision 1.4  2011/02/21  00:35:57  avi
# include dmatinv.o
#
# Revision 1.3  2008/12/02  06:14:44  avi
# eigen -> deigen
#
# Revision 1.2  2008/02/15  22:32:21  avi
# linux compliant
#
# Revision 1.1  2006/11/26  00:06:34  avi
# Initial revision
#
PROG	= covariance
CSRCS	= ${PROG}.c expandf.c
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LOBJS	= ${TRX}/endianio.o ${TRX}/Getifh.o ${RMS}/deigen.o ${RMS}/matopr.o ${RMS}/dmatinv.o

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -I${TRX} -I${RMS} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm -lgfortran
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
	/bin/rm ${OBJS} ${PROG}
