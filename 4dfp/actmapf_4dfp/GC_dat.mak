#$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/GC_dat.mak,v 1.3 2008/09/14 02:56:15 avi Exp $
#$Log: GC_dat.mak,v $
# Revision 1.3  2008/09/14  02:56:15  avi
# linux compliant
#
# Revision 1.2  2007/07/15  03:04:30  avi
# ${NILSRC} ${RELEASE} no -lrms
#
# Revision 1.1  2006/07/26  02:25:00  avi
# Initial revision
#

PROG	= GC_dat
CSRCS	= ${PROG}.c expandf.c
FSRCS	= fGC.f
RMS	= ${NILSRC}/librms
LOBJS	= ${RMS}/matopr.o
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -O
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
