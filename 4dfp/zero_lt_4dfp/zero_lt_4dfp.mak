#$Header: /data/petsun4/data1/src_solaris/zero_lt_4dfp/RCS/zero_lt_4dfp.mak,v 1.7 2007/08/18 23:46:39 avi Exp $
#$Log: zero_lt_4dfp.mak,v $
# Revision 1.7  2007/08/18  23:46:39  avi
# remove redundant LIBS  = -lm
#
# Revision 1.6  2007/05/01  04:40:22  avi
# Solaris10 and linux gcc v3 compliant
#
# Revision 1.5  2004/11/05  22:26:12  rsachs
# Removed 'libmri'. Installed 'get_4d_images2.o', 'Getifh.o'.
#
# Revision 1.4  2001/08/02  01:13:18  avi
# ${PROG} made a prerequisite to release:
#
# Revision 1.3  2001/08/02  00:39:11  avi
#

PROG	= zero_lt_4dfp
CSRCS	= ${PROG}.c
OBJS	= ${CSRCS:.c=.o}
TRX     = ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o
LIBS	= -lm  

CFLAGS = -O -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc ${CFLAGS}
endif

.c.o:
	${CC} -c $<

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS} 

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS}

checkout:
	co ${CSRCS}

checkin:
	ci ${CSRCS}
