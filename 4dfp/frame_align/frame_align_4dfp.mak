#$Header: /data/petsun4/data1/src_solaris/frame_align/RCS/frame_align_4dfp.mak,v 1.3 2011/07/26 04:51:10 avi Exp $
#$Log: frame_align_4dfp.mak,v $
# Revision 1.3  2011/07/26  04:51:10  avi
# add slice_sequence.c
#
# Revision 1.2  2007/09/19  00:56:21  avi
# linux compliant
#
# Revision 1.1  2005/01/28  06:04:55  avi
# Initial revision
#

PROG	= frame_align_4dfp
CSRCS	= ${PROG}.c rfft1d.c slice_sequence.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o
LIBS	= -lm

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
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

clean:
	rm ${OBJS} ${PROG}

release: ${PROG}
	chmod 751 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
