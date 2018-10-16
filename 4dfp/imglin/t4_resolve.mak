# $Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4_resolve.mak,v 1.2 2009/09/09 05:52:14 avi Exp $
# $Log: t4_resolve.mak,v $
# Revision 1.2  2009/09/09  05:52:14  avi
# Solaris10 and linux compliant
#
# Revision 1.1  2009/09/09  03:29:49  avi
# Initial revision
#

PROG	= t4_resolve
CSRCS	= ${PROG}.c
FSRCS	= CKsub.f t4_sub.f ft4_resolve.f getjxy.f
OBJS 	= ${FSRCS:.f=.o} ${CSRCS:.c=.o}
RMS	= ${NILSRC}/librms
LOBJS	= ${RMS}/matopr.o ${RMS}/param6opr.o ${RMS}/eigen.o

CFLAGS	= -O -I.
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore -fcray-pointer
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

${PROG}: ${OBJS}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS}

checkout:
	co ${CSRCS} ${FSRCS} 

checkin:
	ci ${CSRCS} ${FSRCS} 
