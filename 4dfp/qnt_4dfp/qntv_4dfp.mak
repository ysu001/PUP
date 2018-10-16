#$Header: /data/petsun4/data1/src_solaris/qnt_4dfp/RCS/qntv_4dfp.mak,v 1.1 2011/03/04 09:01:36 avi Exp $
#$Log: qntv_4dfp.mak,v $
# Revision 1.1  2011/03/04  09:01:36  avi
# Initial revision
#
PROG	= qntv_4dfp
CSRCS	= ${PROG}.c
TRX     = ${NILSRC}/TRX
ACT     = ${NILSRC}/actmapf_4dfp
RMS     = ${NILSRC}/librms
LOBJS   = ${TRX}/Getifh.o ${TRX}/endianio.o ${TRX}/rec.o ${ACT}/conc.o ${ACT}/expandf.o ${RMS}/dsvdcmp0.o 
OBJS	= ${CSRCS:.c=.o} 
LIBS	= -lm 

.c.o:
	${CC} -c $<

CFLAGS	= -O -I${TRX} -I${ACT}
ifeq (${OSTYPE}, linux)
	CC	= gcc -std=c99 ${CFLAGS}	# -std=c99 enables use of isnormal()
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm -lsunmath
endif

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS} 

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS}

checkin:
	ci ${CSRCS}
