#$Header: /data/petsun4/data1/src_solaris/peak_4dfp/RCS/burn_sphere_4dfp.mak,v 1.5 2007/09/23 03:15:33 avi Exp $
#$Log: burn_sphere_4dfp.mak,v $
# Revision 1.5  2007/09/23  03:15:33  avi
# linux compliant; eliminate -lrms
#
# Revision 1.4  2007/06/06  22:05:43  avi
# add set_rnan.o
#
# Revision 1.3  2006/09/26  02:32:47  avi
# ${PROG} ${RELEASE}
#
# Revision 1.2  2006/08/05  01:50:46  avi
# update ${TRX}
#
# Revision 1.1  2006/03/03  03:44:52  avi
# Initial revision
#

PROG	= burn_sphere_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
SRC	= ${NILSRC}
TRX	= ${SRC}/TRX
RMS	= ${SRC}/librms
REG	= ${SRC}/imgreg_4dfp
T4I	= ${SRC}/t4imgs_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o \
	  ${REG}/imgvalx.o ${T4I}/set_rnan.o ${T4I}/ft4imgo.o \
	  ${RMS}/polfit.o ${RMS}/matopr.o ${RMS}/param6opr.o

CFLAGS	= -I${TRX} -O
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

${PROG}: ${OBJS} 
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
