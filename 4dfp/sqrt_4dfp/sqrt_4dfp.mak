#$Header: /data/petsun4/data1/src_solaris/sqrt_4dfp/RCS/sqrt_4dfp.mak,v 1.4 2008/03/14 02:23:10 avi Exp $
#$Log: sqrt_4dfp.mak,v $
# Revision 1.4  2008/03/14  02:23:10  avi
# linux compliant
#
# Revision 1.3  2006/09/24  05:40:27  avi
# ${PROG} ${RELEASE}
#
# Revision 1.2  2006/08/07  02:38:06  avi
# new ${TRX}
#
# Revision 1.1  2005/01/22  23:06:23  avi
# Initial revision
#

PROG	= sqrt_4dfp
CSRCS	= ${PROG}.c 
OBJS	= ${CSRCS:.c=.o}
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

.c.o:
	${CC} -c $<

CFLAGS	= -O -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc -std=c99 ${CFLAGS}	# -std=c99 enables use of isnormal()
	LIBS	= -lm
else
	CC	= cc ${CFLAGS}
	LIBS	= -lm -lsunmath
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
