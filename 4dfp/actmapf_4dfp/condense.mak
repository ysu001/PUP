#$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/condense.mak,v 1.3 2008/02/01 04:45:52 avi Exp $
#$Log: condense.mak,v $
# Revision 1.3  2008/02/01  04:45:52  avi
# unix compliant
#
# Revision 1.2  2006/10/02  01:55:32  avi
# ${PROG} ${RELEASE}
#
# Revision 1.1  2005/08/30  03:34:11  avi
# Initial revision
#
PROG	= condense
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= -lm

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -O
CC	= cc ${CFLAGS}

${PROG}: ${OBJS} 
	${CC} -o $@ ${OBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
