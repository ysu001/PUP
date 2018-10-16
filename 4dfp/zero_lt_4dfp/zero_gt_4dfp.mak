#$Header: /data/petsun4/data1/src_solaris/zero_lt_4dfp/RCS/zero_gt_4dfp.mak,v 1.3 2007/08/18 22:59:45 avi Exp $
#$Log: zero_gt_4dfp.mak,v $
# Revision 1.3  2007/08/18  22:59:45  avi
# Solaris10 and linux compliant
#
# Revision 1.2  2004/11/05  21:30:36  rsachs
# Removed 'libmri'. Installed 'get_4d_images2.o'.
#
# Revision 1.1  2001/08/02  00:50:56  avi
# Initial revision
#
PROG	= zero_gt_4dfp
CSRCS	= ${PROG}.c
OBJS	= ${CSRCS:.c=.o}
TRX	= ${NILSRC}/TRX/
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o
LIBS	= -lm 

.c.o:
	${CC} -c $<

CFLAGS	= -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
else
	CC	= cc  ${CFLAGS}
endif

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
