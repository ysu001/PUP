#$Header: /data/petsun4/data1/src_solaris/var_4dfp/RCS/dvar_4dfp.mak,v 1.2 2008/02/01 04:16:22 avi Exp $
#$Log: dvar_4dfp.mak,v $
# Revision 1.2  2008/02/01  04:16:22  avi
# linux and Solaris 10 compliant
#
PROG	= dvar_4dfp
CSRCS	= ${PROG}.c
FSRCS	=
OBJS	= ${FSRCS:.f=.o} ${CSRCS:.c=.o}
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

.c.o:
	$(CC) -c $<

.f.o:
	$(FC) -c $<

CFLAGS	= -O -I${TRX}
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
