#$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/compute_defined_4dfp.mak,v 1.3 2008/06/01 04:06:37 avi Exp $
#$Log: compute_defined_4dfp.mak,v $
# Revision 1.3  2008/06/01  04:06:37  avi
# linux compliant
#
# Revision 1.2  2006/10/01  01:24:52  avi
# ${PROG} ${RELEASE}
#
# Revision 1.1  2006/08/10  03:34:41  avi
# Initial revision
#

PROG	= compute_defined_4dfp
CSRCS	= ${PROG}.c conc.c
FSRCS	=
LOBJS	=
TRX	= ${NILSRC}/TRX
ACT	= ${NILSRC}/actmapf_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

.c.o:
	${CC} -c $<

CFLAGS	= -I${TRX} -I${ACT} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc -std=c99 ${CFLAGS}	# -std=c99 enables use of isnormal()
	LIBS	= -lm
else
	CC	= cc ${CFLAGS}
	LIBS	= -lm -lsunmath
endif

${PROG}: ${OBJS} 
	${CC} -o ${PROG} ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}
