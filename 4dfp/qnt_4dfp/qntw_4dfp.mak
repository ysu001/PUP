#$Header: /data/petsun4/data1/src_solaris/qnt_4dfp/RCS/qntw_4dfp.mak,v 1.1 2011/04/18 07:25:35 avi Exp $
#$Log: qntw_4dfp.mak,v $
# Revision 1.1  2011/04/18  07:25:35  avi
# Initial revision
#
PROG	= qntw_4dfp
CSRCS	= ${PROG}.c
TRX     = ${NILSRC}/TRX
ACT     = ${NILSRC}/actmapf_4dfp
LOBJS   = ${TRX}/Getifh.o ${TRX}/endianio.o ${TRX}/rec.o ${ACT}/conc.o
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
