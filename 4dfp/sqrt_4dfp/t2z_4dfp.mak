#$Header: /data/petsun4/data1/src_solaris/sqrt_4dfp/RCS/t2z_4dfp.mak,v 1.3 2008/03/14 03:49:45 avi Exp $
#$Log: t2z_4dfp.mak,v $
# Revision 1.3  2008/03/14  03:49:45  avi
# linux compliant
#
# Revision 1.2  2006/09/25  00:38:04  avi
# ${PROG} ${RELEASE}
#
# Revision 1.1  2006/08/07  02:42:52  avi
# Initial revision
#

PROG	= t2z_4dfp
CSRCS	= ${PROG}.c bertulani_betai.c statsub.c
OBJS	= ${CSRCS:.c=.o}
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

.c.o:
	${CC} -c $<

CFLAGS	= -O -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	LIBS	= -lm
else
	CC	= cc ${CFLAGS}
	LIBS	= -lm
endif

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co $(CSRCS) 

checkin:
	ci $(CSRCS) 
