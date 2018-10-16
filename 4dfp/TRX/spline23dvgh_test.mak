#$Header: /data/petsun4/data1/src_solaris/TRX/RCS/spline23dvgh_test.mak,v 1.6 2011/07/07 01:08:23 avi Exp $
#$Log: spline23dvgh_test.mak,v $
# Revision 1.6  2011/07/07  01:08:23  avi
# restore use of cray-pointer
#
# Revision 1.5  2009/02/25  21:29:25  avi
# accommodate 64bit architecture
#
# Revision 1.4  2007/09/21  06:27:52  avi
# gcc v3 compliant
#
# Revision 1.3  2007/04/22  03:37:52  avi
# set flags for gcc v3
#
# Revision 1.2  2007/04/01  05:37:53  avi
# f90/gcc compatible
#
# Revision 1.1  2005/12/16  04:22:33  avi
# Initial revision
#

PROG	= spline23dvgh_test
CSRCS	= ${PROG}.c
FSRCS	= spline2dvgh.f spline3dvgh.f
RMS	= ${NILSRC}/librms
LOBJS	= ${RMS}/fftsol.o ${RMS}/matopr.o ${RMS}/npad.o

CFLAGS	= -O -I${NILSRC}/TRX
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore -fcray-pointer
	Q	= $(wildcard /usr/lib*/libgfortran.so.1)
	ifeq ($(Q), "")
		LIBS	= -lm -lg2c
	else
		LIBS	= -lm -lgfortran
	endif
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

${PROG}: ${OBJS}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

clean:
	/bin/rm ${OBJS} ${PROG}
