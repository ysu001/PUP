#$Header: /data/petsun4/data1/src_solaris/TRX/RCS/ifh2hdr.mak,v 1.7 2007/05/03 21:43:11 avi Exp $
#$Log: ifh2hdr.mak,v $
# Revision 1.7  2007/05/03  21:43:11  avi
# gcc v3 compliant FC options
#
# Revision 1.6  2007/04/17  05:57:50  avi
# more general gcc compliant strategy
#
# Revision 1.5  2007/04/03  04:36:56  avi
# Linux/Solaris competent
#
# Revision 1.4  2006/03/26  01:36:14  avi
# release directory now environmental variable
#
# Revision 1.3  2006/03/25  04:30:28  avi
# use endianio.c
#
# Revision 1.2  2006/03/16  05:32:36  avi
# cc -I option points to $cwd
#
# Revision 1.1  2005/12/16  02:24:06  avi
# Initial revision
#
PROG	= ifh2hdr
CSRCS	= ${PROG}.c Inithdr.c Getifh.c endianio.c
FSRCS	= 
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= -lm

ifeq (${OSTYPE}, linux)
	CC	= gcc -O -I.
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore
else
	CC	= cc  -O -I.
	FC	= f77 -O -I4 -e
endif

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
