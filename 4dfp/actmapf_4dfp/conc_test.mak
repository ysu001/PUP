#$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/conc_test.mak,v 1.4 2006/09/23 04:36:45 avi Exp $
#$Log: conc_test.mak,v $
# Revision 1.4  2006/09/23  04:36:45  avi
# ${NILSRC} and ${RELEASE}
#
# Revision 1.3  2006/09/23  04:25:33  avi
# include endianio.o in LOBJS
#
# Revision 1.2  2006/08/02  23:59:46  avi
# rec.o now in /data/petsun4/data1/src_solaris/TRX
#
# Revision 1.1  2006/03/16  05:53:15  avi
# Initial revision
#

PROG	= conc_test
CSRCS	= ${PROG}.c conc.c
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
ACT	= ${NILSRC}/actmapf_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o
LIBS	= -lm

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

CFLAGS	= -I${TRX} -I${ACT} -O
CC	= cc ${CFLAGS}
FC 	= f77 -e -I4 -O

${PROG}: ${OBJS} 
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
