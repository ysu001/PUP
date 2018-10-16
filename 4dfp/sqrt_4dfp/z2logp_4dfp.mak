#$Header: /data/petsun4/data1/src_solaris/sqrt_4dfp/RCS/z2logp_4dfp.mak,v 1.3 2008/03/14 04:00:25 avi Exp $
#$Log: z2logp_4dfp.mak,v $
# Revision 1.3  2008/03/14  04:00:25  avi
# linux compliant
#
# Revision 1.2  2006/09/25  01:21:54  avi
# ${PROG} ${RELEASE}
#
# Revision 1.1  2006/08/07  02:46:44  avi
# Initial revision
#

PROG	= z2logp_4dfp
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
