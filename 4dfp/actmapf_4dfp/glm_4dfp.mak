#$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/glm_4dfp.mak,v 1.6 2011/02/06 03:38:33 avi Exp $
#$Log: glm_4dfp.mak,v $
# Revision 1.6  2011/02/06  03:38:33  avi
# eliminate dglm_4dfp.f
#
# Revision 1.5  2007/07/18  02:32:47  avi
# linux compliant
#
# Revision 1.4  2006/09/24  01:08:57  avi
# ${PROG} ${RELEASE}
#
# Revision 1.3  2006/05/04  05:26:14  avi
# new TRX
#
# Revision 1.2  2005/09/05  00:53:21  avi
# double precision GLM inversion using dglm_4dfp.f
#
# Revision 1.1  2004/11/27  05:45:24  avi
# Initial revision
#

PROG	= glm_4dfp
CSRCS	= ${PROG}.c expandf.c conc.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
HST	= ${NILSRC}/img_hist_4dfp
RMS	= ${NILSRC}/librms
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${HST}/fimg_mode.o ${RMS}/deigen.o ${RMS}/dmatinv.o
LIBS	= -lm

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -I. -I${TRX} -O
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
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
