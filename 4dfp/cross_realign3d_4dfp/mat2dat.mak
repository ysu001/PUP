#$Id: mat2dat.mak,v 1.5 2010/02/10 07:00:33 avi Exp $
#$Log: mat2dat.mak,v $
# Revision 1.5  2010/02/10  07:00:33  avi
# include expandf.o in LOBJS
#
# Revision 1.4  2007/08/08  02:39:38  avi
# linux gcc and Solaris 10 compliant
#
# Revision 1.3  2006/09/28  22:00:06  avi
# ${PROG} ${RELEASE}
#
# Revision 1.2  2004/08/02  05:56:42  avi
# routine update
#
# Revision 1.1  1999/01/20  03:34:41  avi
# Initial revision
#

PROG	= mat2dat
CSRCS	= ${PROG}.c
OBJS 	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
RMS	= ${NILSRC}/librms
ACT	= ${NILSRC}/actmapf_4dfp
LOBJS	= ${RMS}/frmsmri.o ${RMS}/matopr.o ${RMS}/param6opr.o ${RMS}/polfit.o ${ACT}/expandf.o

.c.o:
	${CC} -c $<

.f.o:
	f77 -e -c -O $<

CFLAGS	= -O -I${RMS}
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
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS}
