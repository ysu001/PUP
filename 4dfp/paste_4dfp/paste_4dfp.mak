#$Id: paste_4dfp.mak,v 1.3 2007/07/31 06:37:47 avi Exp $
#$Log: paste_4dfp.mak,v $
# Revision 1.3  2007/07/31  06:37:47  avi
# Solaris 10/Linux and endian compliant
#
# Revision 1.2  2004/10/08  22:16:22  rsachs
# Removed 'libmri'. Added 'get_4d_images2.o'.
#
# Revision 1.1  1998/09/25  04:45:08  avi
# Initial revision
#

PROG  	= paste_4dfp
CSRCS 	= ${PROG}.c
FSRCS 	=
OBJS  	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

LIBS 	= -lm 
CFLAGS	= -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
endif

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
	rm ${OBJS}
