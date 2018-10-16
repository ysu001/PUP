#$Header: /home/usr/shimonyj/diff4dfp/RCS/whisker_4dfp.mak,v 1.3 2008/09/16 02:56:14 avi Exp $
#$Log: whisker_4dfp.mak,v $
# Revision 1.3  2008/09/16  02:56:14  avi
# linux compliant
#
# Revision 1.2  2005/12/08  03:10:44  avi
# add code for dealing with affine transformed DWI data (t4tolin.c and affine_DTrot.c)
#
# Revision 1.1  2005/08/06  01:08:32  avi
# Initial revision
#

PROG	= whisker_4dfp
CSRCS	= ${PROG}.c dtensor.c t4tolin.c affine_DTrot.c get_dti_params.c
TRX	= ${NILSRC}/TRX
JSS 	= ${NILSRC}/JSSutil
HST	= ${NILSRC}/img_hist_4dfp
FLP	= ${NILSRC}/flip_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/endianio.o ${TRX}/Getifh.o ${HST}/fimg_mode.o ${FLP}/cflip.o \
	  ${JSS}/JSSnrutil.o ${JSS}/lin_algebra.o
OBJS	= ${CSRCS:.c=.o}
LIBS	= -lm 

CFLAGS	= -O -I. -I${TRX} -I${JSS}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc  ${CFLAGS}
endif

.c.o:
	${CC} -c $<


${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

clean:
	/bin/rm ${OBJS} ${PROG}

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
