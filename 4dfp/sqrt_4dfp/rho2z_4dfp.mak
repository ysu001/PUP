#$Header: /data/petsun4/data1/src_solaris/sqrt_4dfp/RCS/rho2z_4dfp.mak,v 1.4 2007/11/20 22:38:13 avi Exp $
#$Log: rho2z_4dfp.mak,v $
# Revision 1.4  2007/11/20  22:38:13  avi
# linux compliant (eliminate -lsunmath and all references to FORTRAN libs)
#
# Revision 1.3  2006/09/24  23:36:53  avi
# ${PROG} ${RELEASE}
#
# Revision 1.2  2006/08/07  02:34:46  avi
# new ${TRX}
#
# Revision 1.1  2005/09/13  03:22:06  avi
# Initial revision
#

PROG	= rho2z_4dfp
CSRCS	= ${PROG}.c 
OBJS	= ${CSRCS:.c=.o}
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

.c.o:
	${CC} -c $<

CFLAGS	= -O -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc ${CFLAGS}
endif
LIBS	= -lm

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
