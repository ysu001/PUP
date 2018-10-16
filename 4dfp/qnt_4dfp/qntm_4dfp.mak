#$Header: /data/petsun4/data1/src_solaris/qnt_4dfp/RCS/qntm_4dfp.mak,v 1.1 2010/05/19 02:37:14 avi Exp $
#$Log: qntm_4dfp.mak,v $
# Revision 1.1  2010/05/19  02:37:14  avi
# Initial revision
#
PROG	= qntm_4dfp
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
