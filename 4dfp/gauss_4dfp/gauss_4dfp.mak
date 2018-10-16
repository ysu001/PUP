#$Header: /data/petsun4/data1/src_solaris/gauss_4dfp/RCS/gauss_4dfp.mak,v 1.7 2009/02/25 21:18:47 avi Exp $
#$Log: gauss_4dfp.mak,v $
# Revision 1.7  2009/02/25  21:18:47  avi
# accommodate 64bit $ARHC
#
# Revision 1.6  2007/07/04  01:25:34  avi
# test for existence of gcc fortran library
#
# Revision 1.5  2007/05/04  01:56:52  avi
# include -I${RMS}
#
# Revision 1.4  2007/04/23  04:32:06  avi
# gcc v3 compliant (gauss3d.f and gauss3dd.f -> cgauss3d.c cgauss3dd.c)
# eliminate -lrms in link
#
# Revision 1.3  2006/09/25  18:34:12  avi
# ${PROG} ${RELEASE}
#
# Revision 1.2  2005/12/06  06:36:50  avi
# include conc.o in LOBJS

PROG	= gauss_4dfp
CSRCS	= ${PROG}.c cgauss3dd.c
FSRCS	=
ACT	= ${NILSRC}/actmapf_4dfp
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${ACT}/conc.o \
	  ${RMS}/cgauss3d.o ${RMS}/fftsol.o ${RMS}/imgpad.o

CFLAGS	= -O -I${TRX} -I${ACT} -I${RMS}
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
	@echo Q=${Q}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS} 

release: ${PROG} 
	chmod 751 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

.PHONY : clean test
clean:
	rm ${OBJS} ${PROG} 

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

test:
	@echo ${Q}

