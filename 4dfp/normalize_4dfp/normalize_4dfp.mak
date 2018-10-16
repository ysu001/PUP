# $Header: /data/petsun4/data1/src_solaris/normalize_4dfp/RCS/normalize_4dfp.mak,v 1.4 2007/09/19 03:38:48 avi Exp $
# $Log: normalize_4dfp.mak,v $
# Revision 1.4  2007/09/19  03:38:48  avi
# linux compliant
# eliminate -lrms
#
# Revision 1.3  2006/09/29  17:06:51  avi
# ${PROG} ${RELEASE}
#
# Revision 1.2  2005/12/19  02:20:22  avi
#  remove links to libmri
#
# Revision 1.1  1999/03/24  00:28:31  avi
# Initial revision
#
PROG	= normalize_4dfp
CSRCS	= ${PROG}.c img2hist.c
FSRCS	= fnormalize.f
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o \
	  ${RMS}/imag2mask.o ${RMS}/imgpad.o ${RMS}/gauss3d.o ${RMS}/fftsol.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

CFLAGS	= -I${TRX} -I${RMS} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
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

${PROG}: $(OBJS)
	${FC} -o $@ $(OBJS) $(LOBJS) ${LIBS}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
