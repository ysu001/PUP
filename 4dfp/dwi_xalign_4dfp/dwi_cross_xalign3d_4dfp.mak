####################################################################################
# Copyright 2003, 2004, 2005, 2006, 2007, 2008
# Washington University, Mallinckrodt Institute of Radiology.
# All Rights Reserved.
# This software may not be reproduced, copied, or distributed without written
# permission of Washington University. For further information contact A. Z. Snyder.
####################################################################################
#$Id: dwi_cross_xalign3d_4dfp.mak,v 1.4 2008/08/15 04:17:14 avi Exp $
#$Log: dwi_cross_xalign3d_4dfp.mak,v $
# Revision 1.4  2008/08/15  04:17:14  avi
# linux compliant
#
# Revision 1.3  2007/07/22  02:58:57  avi
# Solaris 10
#
# Revision 1.2  2005/12/30  00:23:52  avi
# get rec.h from imglin
# include fimgetae.o in LOBJS
#
# Revision 1.1  2005/11/29  08:08:20  avi
# Initial revision
#

PROG	= dwi_cross_xalign3d_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
LOBJS	= ${TRX}/fimgeta.o ${TRX}/fimgetae.o ${TRX}/etadeta.o ${TRX}/fomega.o ${TRX}/twoecal.o \
	  ${TRX}/Getifh.o ${TRX}/endianio.o ${TRX}/spline2dvgh.o ${TRX}/spline3dvgh.o ${TRX}/rec.o \
	  ${RMS}/fftsol.o ${RMS}/imgpad.o ${RMS}/matopr.o ${RMS}/gauss3d.o ${RMS}/eigen.o ${RMS}/determ12.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= -lm

CFLAGS	= -I${TRX} -I${RMS} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore
	LIBS	= -lm -lgfortran
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

${PROG}: ${OBJS}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}
