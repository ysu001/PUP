# $Header: /data/petsun4/data1/src_solaris/deband_4dfp/RCS/deband_4dfp.mak,v 1.6 2009/02/22 02:15:35 avi Exp ${LINSRC}/deband_4dfp/RCS/deband_4dfp.mak,v 1.4 2006/02/01 03:58:41 avi Exp avi $
# $Log: deband_4dfp.mak,v $
# Revision 1.6  2009/02/22  02:15:35  avi
# linux compliant
# imag2mask (in former librms library) -> img2lmask.o in librms subdirectory
#
# Revision 1.5  2006/09/27  22:17:46  avi
# ${PROG} ${RELEASE}
#
# Revision 1.4  2006/02/01  03:58:41  avi
# TRX moved
#
# Revision 1.3  2004/10/12  20:12:20  rsachs
# Added 'writeifh.o'.
#
# Revision 1.2  2004/10/11  22:03:48  rsachs
# Removed 'libmri'. Added 'Getifh.o','get_4d_images2.o'.
#
# Revision 1.1  1999/03/05  05:44:06  avi
# Initial revision
#

PROG	= deband_4dfp
CSRCS	= ${PROG}.c
FSRCS	= fdeband.f debandg.f
TRX	= ${NILSRC}/TRX
ACT	= ${NILSRC}/actmapf_4dfp
RMS	= ${NILSRC}/librms
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${ACT}/expandf.o \
	  ${RMS}/img2lmask.o ${RMS}/imgpad.o ${RMS}/fftsol.o ${RMS}/gauss3d.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

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

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

${PROG}: $(OBJS)
	${FC} -o $@ ${OBJS} ${LOBJS} ${FOBJS} ${LIBS}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

release: ${PROG}
	chmod 711 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
