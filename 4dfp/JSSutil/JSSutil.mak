#$Id: JSSutil.mak,v 1.2 2008/12/30 03:28:58 avi Exp $
#$Log: JSSutil.mak,v $
# Revision 1.2  2008/12/30  03:28:58  avi
# include make -f JSSstatistics.c
#
# Revision 1.1  2007/08/28  03:43:44  avi
# Initial revision
#

CSRCS	= JSSnrutil.c lin_algebra.c random.c JSSstatistics.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -O -I.
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc ${CFLAGS}
endif

compile: ${OBJS}

clean:
	rm ${OBJS}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}
