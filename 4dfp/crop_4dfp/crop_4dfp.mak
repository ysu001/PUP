#$Id: crop_4dfp.mak,v 1.5 2010/08/20 04:59:10 avi Exp $
#$Log: crop_4dfp.mak,v $
# Revision 1.5  2010/08/20  04:59:10  avi
# look for flip_4dfp.h
#
# Revision 1.4  2008/05/05  04:14:48  avi
# Solaris 10 and linux compliant
#
# Revision 1.3  2004/11/18  21:39:34  rsachs
# Removed 'libmri'. Installed 'writeifh.c' & 'get_4d_images2.o'.
#
# Revision 1.2  2004/02/19  01:14:29  avi
# eliminate FORTRAN dependence (fflip.o)
#
# Revision 1.1  2002/10/19  23:05:47  avi
# Initial revision
#

PROG 	= crop_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
TRX     = ${NILSRC}/TRX
FLP	= ${NILSRC}/flip_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${FLP}/cflip.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

.c.o:
	${CC} -c $<

CFLAGS	= -I${TRX} -I${FLP} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc  ${CFLAGS}
endif
LIBS	= -lm

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}
