#$Header: /data/petsun4/data1/src_solaris/peak_4dfp/RCS/index2atl.mak,v 1.6 2009/02/25 21:29:27 avi Exp $
#$Log: index2atl.mak,v $
# Revision 1.6  2009/02/25  21:29:27  avi
# accommodate 64bit architecture
#
# Revision 1.5  2007/09/23  03:35:36  avi
# gcc v3 and v4 compliant
#
# Revision 1.4  2007/05/01  20:18:58  avi
# typo
#
# Revision 1.3  2007/05/01  20:13:36  avi
# linux gcc v3 compliant (eliminate -lrms)
#
# Revision 1.2  2006/09/26  03:12:40  avi
# ${PROG} ${RELEASE}
#
# Revision 1.1  2002/10/21  17:40:47  avi
# Initial revision
#

PROG	= index2atl
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
RMS	= ${NILSRC}/librms
FT4	= ${NILSRC}/t4imgs_4dfp
REG	= ${NILSRC}/imgreg_4dfp
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/Getifh.o ${TRX}/endianio.o ${REG}/imgvalx.o ${FT4}/ft4imgo.o  ${FT4}/set_rnan.o \
	  ${RMS}/param6opr.o ${RMS}/matopr.o ${RMS}/polfit.o

CFLAGS	= -O -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore
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

${PROG}: ${OBJS} 
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
