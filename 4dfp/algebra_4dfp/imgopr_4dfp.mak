#$Header: /data/petsun4/data1/src_solaris/algebra_4dfp/RCS/imgopr_4dfp.mak,v 1.5 2007/07/05 04:42:10 avi Exp $
#$Log: imgopr_4dfp.mak,v $
# Revision 1.5  2007/07/05  04:42:10  avi
# Linux gcc compliant
#
# Revision 1.4  2006/09/24  02:55:35  avi
# ${PROG} ${RELEASE}
#
# Revision 1.3  2006/09/13  02:08:43  avi
# rec in /data/petsun4/data1/src_solaris/TRX
#
# Revision 1.2  2006/03/16  05:44:23  avi
# all #include now local
#
# Revision 1.1  2004/03/11  05:42:23  avi
# Initial revision
#

PROG	= imgopr_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

CFLAGS	= -I${TRX} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc -std=c99 ${CFLAGS}	# -std=c99 enables use of isnormal()
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm -lsunmath
endif

.c.o:
	${CC} -c $<
.f.o:
	${FC} -c $<

${PROG}: imgopr_4dfp.o 
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 755 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
