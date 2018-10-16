#$Header: /data/petsun4/data1/src_solaris/TRX/RCS/endian_4dfp.mak,v 1.4 2007/05/03 21:44:11 avi Exp $
#$Log: endian_4dfp.mak,v $
# Revision 1.4  2007/05/03  21:44:11  avi
# gcc v3 FC compliant options
#
# Revision 1.3  2007/04/17  05:54:37  avi
# more general gcc compliant strategy
#
# Revision 1.2  2007/04/03  03:45:59  avi
# Linux/Solaris competent
#
# Revision 1.1  2006/03/26  01:33:07  avi
# Initial revision
#
PROG	= endian_4dfp
CSRCS	= ${PROG}.c Getifh.c endianio.c rec.c
FSRCS	= 
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

ifeq (${OSTYPE}, linux)
	CC	= gcc -O -I.
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore
	LIBS	= -lm
else
	CC	= cc -O -I.
	FC	= f77 -O -I4 -e
	LIBS	= -lm -lsunmath
endif

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LIBS}

endian_4dfp.o:
ifeq (${OSTYPE}, linux)
	${CC} -std=c99 -c endian_4dfp.c
else
	${CC}          -c endian_4dfp.c
endif

clean:
	/bin/rm ${PROG} ${OBJS}

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
