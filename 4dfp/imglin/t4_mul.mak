#$Id: t4_mul.mak,v 1.4 2009/02/25 21:29:26 avi Exp $
#$Log: t4_mul.mak,v $
# Revision 1.4  2009/02/25  21:29:26  avi
# accommodate 64bit architecture
#
# Revision 1.3  2007/09/21  21:43:24  avi
# gcc v3 and v4 compliant
#
# Revision 1.2  2007/04/30  17:54:34  avi
#

PROG	= t4_mul
CSRCS	= ${PROG}.c
FSRCS	=
RMS	= ${NILSRC}/librms
LOBJS	= ${RMS}/matopr.o
OBJS	= ${FSRCS:.f=.o} ${CSRCS:.c=.o}

CFLAGS	= -O -I${RMS}
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fno-second-underscore
	Q	= $(wildcard /usr/lib*/libgfortran.so.1)
	ifeq ($(Q), "")
		LIBS	= -lm -lg2c
	else
		LIBS	= -lm -lgfortran
	endif
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
	LIBS	= -lm
endif

${PROG}: ${OBJS}
	$(FC) -o $@ ${OBJS} ${LOBJS} ${LIBS}

.c.o:
	$(CC) -c $<

.f.o:
	$(FC) -c $<

clean:
	/bin/rm ${OBJS} ${PROG}

release: ${PROG} 
	chmod 771 ${PROG}
	mv ${PROG} ${RELEASE}

checkout:
	co $(CSRCS) $(FSRCS) 

checkin:
	ci $(CSRCS) $(FSRCS) 
