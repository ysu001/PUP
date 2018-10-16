#$Header: /data/petsun4/data1/src_solaris/collate_slice_4dfp/RCS/collate_slice_4dfp.mak,v 1.3 2007/09/19 02:10:20 avi Exp $
#$Log: collate_slice_4dfp.mak,v $
# Revision 1.3  2007/09/19  02:10:20  avi
#
# Revision 1.2  2007/09/19  02:06:21  avi
# linux, Solaris 10 compliant
#
# Revision 1.1  2004/11/09  21:56:53  rsachs
# Initial revision
#

PROG	= collate_slice_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/Getifh.o ${TRX}/endianio.o ${TRX}/rec.o

CFLAGS	= -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc  ${CFLAGS}
endif

.c.o:
	${CC} -c $<
.f.o:
	${FC} -c $<

${PROG}: ${OBJS} 
	${CC} -o $@ ${OBJS} ${LOBJS} -lm

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
