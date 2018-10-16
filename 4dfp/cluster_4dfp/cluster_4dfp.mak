#$Header: /data/petsun4/data1/src_solaris/cluster_4dfp/RCS/cluster_4dfp.mak,v 1.3 2009/07/02 05:50:57 avi Exp $
#$Log: cluster_4dfp.mak,v $
# Revision 1.3  2009/07/02  05:50:57  avi
# include cflip.o in ${LOBJS}
#
# Revision 1.2  2008/04/28  05:15:42  avi
# linux compliant
#
# Revision 1.1  2006/11/04  06:24:09  avi
# Initial revision
#
PROG	= cluster_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX     = ${NILSRC}/TRX
FLP     = ${NILSRC}/flip_4dfp
LOBJS   = ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${FLP}/cflip.o
LIBS	= -lm

.c.o:
	${CC} -c $<

CFLAGS	= -O -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc ${CFLAGS}
endif
LIBS	= -lm

${PROG}: $(OBJS)
	${CC} -o $@ $(OBJS) $(LOBJS) ${LIBS}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

release: ${PROG}
	chmod 711 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
