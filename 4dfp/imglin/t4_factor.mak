#$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4_factor.mak,v 1.5 2009/02/25 21:29:26 avi Exp $
#$Log: t4_factor.mak,v $
# Revision 1.5  2009/02/25  21:29:26  avi
# accommodate 64bit architecture
#
# Revision 1.4  2009/02/24  01:01:11  avi
# correct confusing variable name (TRX -> RMS)
#
# Revision 1.3  2009/02/23  06:28:13  avi
# correct ${LIBS} usage
#
# Revision 1.2  2009/02/23  06:18:31  avi
# linux compliant with C instead of FORTRAN main
#
# Revision 1.1  2006/09/29  00:23:37  avi
# Initial revision

PROG	= t4_factor
CSRCS	= ${PROG}.c t4_io.c
FSRCS	= t4_sub.f param12opr.f
OBJS 	= ${FSRCS:.f=.o} ${CSRCS:.c=.o}
RMS	= ${NILSRC}/librms
LOBJS	= ${RMS}/matopr.o ${RMS}/eigen.o ${RMS}/param6opr.o

CFLAGS	= -O -I.
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
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

.c.o:
	$(CC) -c $<

.f.o:
	$(FC) -c $<

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co $(CSRCS) $(FSRCS) 

checkin:
	ci $(CSRCS) $(FSRCS) 

