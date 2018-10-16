#$Header: /data/petsun4/data1/src_solaris/analyzeto4dfp/RCS/hdr2txt.mak,v 1.2 2007/09/26 02:29:26 avi Exp $
#$Log: hdr2txt.mak,v $
# Revision 1.2  2007/09/26  02:29:26  avi
# set group program on release
#
# Revision 1.1  2007/09/21  23:36:19  avi
# Initial revision
#
PROG	= hdr2txt
CSRCS	= ${PROG}.c
FSRCS	=
TRX	= ${NILSRC}/TRX
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LOBJS	= ${TRX}/endianio.o ${TRX}/Getifh.o
LIBS	= -lm

CFLAGS	= -O -I${TRX}
CC	= cc ${CFLAGS}

.c.o:
	${CC} -c $<
.f.o:
	${FC} -c $<

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS} 

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
