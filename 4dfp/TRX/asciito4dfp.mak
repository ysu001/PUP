#$Header: /data/petsun4/data1/src_solaris/TRX/RCS/asciito4dfp.mak,v 1.1 2007/09/08 21:00:51 avi Exp $
#$Log: asciito4dfp.mak,v $
# Revision 1.1  2007/09/08  21:00:51  avi
# Initial revision
#
PROG	= asciito4dfp
CSRCS	= ${PROG}.c Getifh.c endianio.c rec.c
FSRCS	= 
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

CC	= cc -O -I.

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LIBS}

clean:
	/bin/rm ${PROG} ${OBJS}

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
