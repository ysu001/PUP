#$Header: /home/usr/shimonyj/diff4dfp/RCS/diff_4dfp.mak,v 1.6 2010/01/25 05:48:43 avi Exp $
#$Log: diff_4dfp.mak,v $
# Revision 1.6  2010/01/25  05:48:43  avi
# include "knee" subroutines
#
# Revision 1.5  2008/09/16  03:14:14  avi
# ${JSS}/JSSnrutil.o ${JSS}/lin_algebra.o ${JSS}/random.o
#
# Revision 1.4  2008/06/11  03:29:18  avi
# remove sources bayes_mcmc.c hist.c
#
# Revision 1.3  2008/04/14  00:13:37  adrian
# add pointer to JSSutil.h
#
# Revision 1.2  2007/10/15  21:20:43  avi
# endianio.c Getifh.c compliant
#
# Revision 1.1  2006/03/18  07:39:12  avi
# Initial revision
#

PROG	= diff_4dfp
CSRCS	= ${PROG}.c get_dti_params.c dtensor.c nonlinear.c feat.c morph.c fuzzy.c contour.c
TRX	= ${NILSRC}/TRX
JSS 	= ${NILSRC}/JSSutil
HST	= ${NILSRC}/img_hist_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/endianio.o ${TRX}/Getifh.o \
	  ${JSS}/JSSnrutil.o ${JSS}/lin_algebra.o ${JSS}/random.o \
	  ${HST}/fimg_mode.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= -lm

CFLAGS	= -O -I${TRX} -I. -I${JSS}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc  ${CFLAGS}
endif

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

clean:
	/bin/rm ${OBJS} ${PROG}

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
