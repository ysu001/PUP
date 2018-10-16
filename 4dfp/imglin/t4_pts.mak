#$Id: t4_pts.mak,v 1.4 2009/10/27 05:57:11 avi Exp $
#$Log: t4_pts.mak,v $
# Revision 1.4  2009/10/27  05:57:11  avi
# correct ${LIBS}
#
# Revision 1.3  2009/10/27  05:44:02  avi
# linux compliant C main
#
# Revision 1.2  2008/01/26  23:15:24  avi
# Solaris 10 compliant
#
# Revision 1.1  2008/01/26  23:07:42  avi
# Initial revision
#

PROG	= t4_pts
CSRCS	= ${PROG}.c t4_io.c 
FSRCS	= param12opr.f
OBJS	= ${FSRCS:.f=.o} ${CSRCS:.c=.o}
RMS	= ${NILSRC}/librms
LOBJS	= ${RMS}/param6opr.o ${RMS}/matopr.o ${RMS}/eigen.o

t4_pts: ${OBJS}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

CFLAGS	= -I. -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore
	LIBS	= -lm -lgfortran
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

release: ${PROG}
	chmod 751 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS} 

checkin:
	ci ${CSRCS} ${FSRCS}
