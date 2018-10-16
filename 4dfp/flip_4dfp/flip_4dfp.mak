#$Id: flip_4dfp.mak,v 1.4 2007/05/02 01:27:30 avi Exp $
#$Log: flip_4dfp.mak,v $
# Revision 1.4  2007/05/02  01:27:30  avi
# linux gcc v3 compliant; $NILSRC and $RELEASE
#
# Revision 1.3  2004/11/16  22:27:32  rsachs
# Removed 'libmri'. Installed 'get_4d_images2.o'.
#
# Revision 1.2  2004/02/19  01:05:02  avi
# eliminate FORTRAN dependence
#
# Revision 1.1  1999/01/24  07:07:16  avi
# Initial revision
#

PROG	= flip_4dfp
CSRCS	= ${PROG}.c cflip.c
FSRCS	=
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= -lm

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

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
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

