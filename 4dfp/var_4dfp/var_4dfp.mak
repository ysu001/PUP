#$Header: /data/petsun4/data1/src_solaris/var_4dfp/RCS/var_4dfp.mak,v 1.8 2007/08/05 01:22:03 avi Exp $
#$Log: var_4dfp.mak,v $
# Revision 1.8  2007/08/05  01:22:03  avi
# linux compliant
#
# Revision 1.7  2006/09/24  02:21:30  avi
# ${PROG} ${RELEASE}
#

PROG	= var_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${FSRCS:.f=.o} ${CSRCS:.c=.o}
TRX	= ${NILSRC}/TRX
ACT	= ${NILSRC}/actmapf_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${ACT}/expandf.o ${ACT}/conc.o

.c.o:
	$(CC) -c $<

.f.o:
	$(FC) -c $<

CFLAGS	= -O -I${ACT} -I${TRX}
ifeq (${OSTYPE}, linux)
	CC	= gcc -std=c99 ${CFLAGS}	# -std=c99 enables use of isnormal()
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
	LIBS	= -lm
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

${PROG}: $(OBJS)
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS}

checkout:
	co $(CSRCS) $(FSRCS) 

checkin:
	ci $(CSRCS) $(FSRCS) 
