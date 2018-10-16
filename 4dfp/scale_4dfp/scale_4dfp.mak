#$Id: scale_4dfp.mak,v 1.5 2007/05/02 05:32:11 avi Exp $
#$Log: scale_4dfp.mak,v $
# Revision 1.5  2007/05/02  05:32:11  avi
# linux compliant
#
# Revision 1.4  2006/09/29  23:37:58  avi
# ${PROG} ${RELEASE}
#
# Revision 1.3  2004/09/23  20:59:28  rsachs
# Removed 'libmri'. Added 'Getifh.o', 'get_4d_images2.o'.
#
# Revision 1.2  2004/05/10  02:05:58  avi
# compile cc and set group to program on release
#
# Revision 1.1  1998/10/02  06:31:10  avi
# Initial revision
#

PROG	= scale_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
TRX     = ${NILSRC}/TRX
OBJS 	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LOBJS   = ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o 

CFLAGS	= -I${TRX}
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

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}
