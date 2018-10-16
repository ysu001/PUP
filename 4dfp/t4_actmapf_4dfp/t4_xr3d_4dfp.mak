#$Header: /data/petsun4/data1/src_solaris/t4_actmapf_4dfp/RCS/t4_xr3d_4dfp.mak,v 1.4 2007/09/20 03:15:28 avi Exp $
#$Log: t4_xr3d_4dfp.mak,v $
# Revision 1.4  2007/09/20  03:15:28  avi
# linux compliant (fftsun -> fftsol)
# eliminate -lrms
#
# Revision 1.3  2006/09/29  04:55:09  avi
# ${PROG} ${RELEASE}
#
# Revision 1.2  2006/09/29  02:48:43  avi
# rsachs modernization
#
# Revision 1.1  1999/01/12  06:38:02  avi
# Initial revision
#

PROG	= t4_xr3d_4dfp
CSRCS	= ${PROG}.c
FSRCS	= ft4ixyz.f
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
FLP	= ${NILSRC}/flip_4dfp
T4I	= ${NILSRC}/t4imgs_4dfp
LOBJS	= ${TRX}/rec.o ${FLP}/cflip.o ${TRX}/Getifh.o ${TRX}/endianio.o ${NILSRC}/imglin/t4_sub.o \
	  ${T4I}/ft4imgo.o ${T4I}/set_rnan.o ${TRX}/spline3dvgh.o ${NILSRC}/imgreg_4dfp/imgvalx.o \
	  ${RMS}/param6opr.o ${RMS}/npad.o ${RMS}/matopr.o ${RMS}/polfit.o ${RMS}/fftsol.o

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
	chmod 771 ${PROG}
	mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
