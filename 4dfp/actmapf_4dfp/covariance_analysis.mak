#$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/covariance_analysis.mak,v 1.3 2008/09/20 00:32:13 avi Exp $
#$Log: covariance_analysis.mak,v $
# Revision 1.3  2008/09/20  00:32:13  avi
# lose all mention of FORTRAN
#
# Revision 1.2  2006/11/25  23:51:08  avi
# typo
#
# Revision 1.1  2006/11/25  23:49:30  avi
# Initial revision
#
PROG	= covariance_analysis
CSRCS	= ${PROG}.c
OBJS	= ${CSRCS:.c=.o}
LIBS	= -lm

.c.o:
	${CC} -c $<

CFLAGS	= -O
CC	= cc ${CFLAGS}

${PROG}: ${OBJS} 
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
