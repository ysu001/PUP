#$Id: 4dfptoanalyze.mak,v 1.5 2007/09/21 04:17:48 avi Exp $
#$Log: 4dfptoanalyze.mak,v $
# Revision 1.5  2007/09/21  04:17:48  avi
# complete clean
#
# Revision 1.4  2007/05/01  02:53:27  avi
# linux gcc compliant
#
# Revision 1.3  2007/02/28  07:10:11  avi
# Solaris 10
#
# Revision 1.2  2004/09/10  05:10:04  avi
# eliminate all FORTRAN references
#
# Revision 1.1  1998/12/30  05:32:50  avi
# Initial revision
#

PROG	= 4dfptoanalyze
CSRCS   = ${PROG}.c
FSRCS   =
TRX     = ${NILSRC}/TRX
FLP	= ${NILSRC}/flip_4dfp
LOBJS   = ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${FLP}/cflip.o ${TRX}/Inithdr.o
OBJS    = ${FSRCS:.f=.o} ${CSRCS:.c=.o}
LIBS    = -lm

CFLAGS = -O -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc ${CFLAGS}
endif
LIBS	= -lm

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

clean:
	rm ${OBJS} ${PROG}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

checkout:
	co $(CSRCS) $(FSRCS) 

checkin:
	ci $(CSRCS) $(FSRCS) 
