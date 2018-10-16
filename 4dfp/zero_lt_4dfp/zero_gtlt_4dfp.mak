#$Header: /data/petsun4/data1/src_solaris/zero_lt_4dfp/RCS/zero_gtlt_4dfp.mak,v 1.2 2007/08/19 00:10:05 avi Exp $
#$Log: zero_gtlt_4dfp.mak,v $
# Revision 1.2  2007/08/19  00:10:05  avi
# Solaris 10 and linux compliant
#
# Revision 1.1  2004/11/08  20:08:04  rsachs
# Initial revision
#

PROG	= zero_gtlt_4dfp
CSRCS	= ${PROG}.c
OBJS	= ${CSRCS:.c=.o}
TRX	= ${NILSRC}/TRX/
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o
LIBS	= -lm 

.c.o:
	${CC} -c $<

CFLAGS	= -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc  ${CFLAGS}
endif

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS} 

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS}

checkout:
	co ${CSRCS}

checkin:
	ci ${CSRCS}
