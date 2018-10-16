#$Header: /home/usr/shimonyj/diff4dfp/RCS/diffRGB_4dfp.mak,v 1.7 2008/09/16 03:16:50 avi Exp $
#$Log: diffRGB_4dfp.mak,v $
# Revision 1.7  2008/09/16  03:16:50  avi
# linux compliant
#
# Revision 1.6  2008/04/14  01:50:38  adrian
# add pointer to JSSutil.h
#
# Revision 1.5  2006/09/04  06:27:29  avi
# new TRX and endianio modules
#
# Revision 1.4  2006/02/21  05:07:33  adrian
# install calls to t4tolin() and affine_DTrot() to deal with affine transformed DWI data
#
# Revision 1.3  2004/06/10  03:39:38  avi
# eliminate -lmri and dependency on FORTRAN
#
# Revision 1.2  2003/09/17  22:23:14  avi
# typos
#
# Revision 1.1  2003/09/17  22:21:10  avi
# Initial revision
#

PROG	= diffRGB_4dfp 
CSRCS	= diffRGB_4dfp.c dtensor.c t4tolin.c affine_DTrot.c get_dti_params.c
FSRCS	= 
JSS 	= ${NILSRC}/JSSutil
HST	= ${NILSRC}/img_hist_4dfp
TRX	= ${NILSRC}/TRX
FLP	= ${NILSRC}/flip_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${FLP}/cflip.o ${TRX}/Inithdr.o \
	  ${TRX}/endianio.o ${HST}/fimg_mode.o \
	  ${JSS}/JSSnrutil.o ${JSS}/lin_algebra.o 
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= -lm 

CFLAGS	= -O -I. -I${TRX} -I${JSS}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc  ${CFLAGS}
endif

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

checkin:
	ci ${CSRCS} ${FSRCS} 

checkout:
	co ${CSRCS} ${FSRCS} 

clean:
	/bin/rm ${OBJS} ${PROG}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
