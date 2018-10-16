#$Id: rmspike_4dfp.mak,v 1.8 2007/04/17 06:04:47 avi Exp $
#$Log: rmspike_4dfp.mak,v $
# Revision 1.8  2007/04/17  06:04:47  avi
# gcc compliant
#
# Revision 1.7  2007/04/11  01:55:27  avi
# linux ready
#
# Revision 1.6  2006/09/27  20:36:10  avi
# rmspike_4dfp.mak
#
# Revision 1.5  2005/07/07  00:42:19  avi
# add ${NON}/writeifh.o to LOBJS
#
# Revision 1.4  2004/10/11  22:06:33  rsachs
# Removed 'libmri'. Installed 'Getifh.o','get_4d_images2.o'.
#
# Revision 1.3  2001/05/17  07:06:03  avi
#
# Revision 1.2  2001/05/17  07:03:39  avi
# expandf.o
#
# Revision 1.1  2000/11/14  07:35:46  avi
# Initial revision
#

PROG	= rmspike_4dfp
CSRCS	= ${PROG}.c
FSRCS	= frmspike.f
OBJS 	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
ACT     = ${NILSRC}/actmapf_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${ACT}/expandf.o ${TRX}/endianio.o

ifeq (${OSTYPE}, linux)
	CC	= gcc -O -I${TRX}
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm -lgfortran
else
	CC	= cc  -O -I${TRX}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

${PROG}: ${OBJS}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release:${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}
