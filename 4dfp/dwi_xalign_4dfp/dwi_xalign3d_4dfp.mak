#########################################################################################
# Copyright 2003, 2004, 2005, 2006, 2007						#
# Washington University, Mallinckrodt Institute of Radiology.				#
# All Rights Reserved.									#
# This software may not be reproduced, copied, or distributed without written		#
# permission of Washington University. For further information contact A. Z. Snyder.	#
#########################################################################################
#$Id: dwi_xalign3d_4dfp.mak,v 1.5 2009/02/25 21:29:25 avi Exp $
#$Log: dwi_xalign3d_4dfp.mak,v $
# Revision 1.5  2009/02/25  21:29:25  avi
# accommodate 64bit architecture
#
# Revision 1.4  2008/08/15  03:50:04  avi
# linux compliant
#
# Revision 1.3  2007/06/21  01:03:05  avi
# eliminate -lrms
#
# Revision 1.2  2007/03/12  02:46:11  avi
# Solaris 10
#
# Revision 1.1  2005/12/29  07:37:20  avi
# Initial revision
#

PROG	= dwi_xalign3d_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
LOBJS	= ${TRX}/spline2dvgh.o ${TRX}/spline3dvgh.o ${TRX}/rec.o ${TRX}/endianio.o ${TRX}/Getifh.o \
	  ${TRX}/fimgeta.o ${TRX}/fimgetae.o ${TRX}/etadeta.o ${TRX}/fomega.o ${TRX}/twoecal.o \
	  ${RMS}/eigen.o ${RMS}/matopr.o ${RMS}/determ12.o ${RMS}/imgpad.o ${RMS}/fftsol.o ${RMS}/gauss3d.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= -lm

CFLAGS	= -I${TRX} -I${RMS} -O
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

${PROG}: ${OBJS}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

release: ${PROG}
	chmod 775 ${PROG}
	mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}
