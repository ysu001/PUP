#$Id: reindex_4dfp.mak,v 1.4 2007/08/15 23:54:01 avi Exp $
#$Log: reindex_4dfp.mak,v $
# Revision 1.4  2007/08/15  23:54:01  avi
# Solaris10, endian and linux compliant
#
# Revision 1.3  2004/09/23  20:29:54  rsachs
# Removed 'libmri'. Added 'get_4d_images2.o'.
#
# Revision 1.2  2004/05/02  02:04:58  avi
# release with group set to program
#
# Revision 1.1  2000/09/24  22:32:33  avi
# Initial revision
#

PROG	= reindex_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS 	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

LIBS	= -lm 

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -I. -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm -lgfortran
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release:${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}
