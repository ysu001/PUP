#$Header: /data/petsun4/data1/src_solaris/cross_realign3d_4dfp/RCS/cross_realign3d_4dfp.mak,v 1.8 2013/02/14 03:23:29 avi Exp $
#$Log: cross_realign3d_4dfp.mak,v $
# Revision 1.8  2013/02/14  03:23:29  avi
# remove redundant .c.o .c.f
#
# Revision 1.7  2007/08/08  00:32:15  avi
# linux compliant
#
# Revision 1.6  2007/04/16  04:34:39  avi
# Solaris 10
#
# Revision 1.5  2006/03/27  04:24:06  avi
# endian invariant i/o subroutines
#
# Revision 1.4  2005/08/28  00:01:49  avi
# eliminate -lmri
#
# Revision 1.3  2003/08/17  02:16:44  avi
# eliminate references to cross_realign3d_4dfp.h
#
# Revision 1.2  1999/03/08  02:27:11  avi
# expandf ()
#
# Revision 1.1  1999/01/20  03:31:43  avi
# Initial revision

PROG	= cross_realign3d_4dfp
CSRCS	= ${PROG}.c gauss3diz.c gauss2diz.c
FSRCS	= frmsfmri.f splinexyz.f alignmrixyz.f imgvalsixyz.f	# frmsfmri.f replaces frmsmri.f in librms
TRX	= ${NILSRC}/TRX
ACT	= ${NILSRC}/actmapf_4dfp
RMS	= ${NILSRC}/librms
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${ACT}/expandf.o \
	  ${RMS}/matopr.o ${RMS}/imgpad.o ${RMS}/fftsol.o ${RMS}/img2lmask.o ${RMS}/polfit.o ${RMS}/gauss3d.o ${RMS}/param6opr.o 
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= lm

CFLAGS	= -O -I${TRX} -I${RMS}
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
	chmod 771 ${PROG}
	mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}
