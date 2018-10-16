#$Header: /data/petsun4/data1/src_solaris/peak_4dfp/RCS/peak_4dfp.mak,v 1.5 2007/09/23 03:54:57 avi Exp $
#$Log: peak_4dfp.mak,v $
# Revision 1.5  2007/09/23  03:54:57  avi
# linux compliant
#
# Revision 1.4  2007/06/04  04:16:03  avi
# include ${SRC}/t4imgs_4dfp/set_rnan.o called by ${SRC}/t4imgs_4dfp/ft4imgo.o
#
# Revision 1.3  2006/09/26  01:16:10  avi
# ${PROG} ${RELEASE}
#
# Revision 1.2  2006/09/07  04:58:05  avi
# TRX moved to /data/petsun4/data1/src_solaris
#
# Revision 1.1  2002/10/21  17:40:35  avi
# Initial revision
#

PROG	= peak_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
SRC	= ${NILSRC}
RMS	= ${NILSRC}/librms
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o \
	  ${SRC}/imgreg_4dfp/imgvalx.o ${SRC}/t4imgs_4dfp/ft4imgo.o ${SRC}/t4imgs_4dfp/set_rnan.o \
	  ${RMS}/polfit.o ${RMS}/matopr.o ${RMS}/param6opr.o

CFLAGS	= -I${TRX} -O
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

${PROG}: ${OBJS} 
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
