#$Id: T2C_4dfp.mak,v 1.5 2009/02/25 21:29:25 avi Exp $
#$Log: T2C_4dfp.mak,v $
# Revision 1.5  2009/02/25  21:29:25  avi
# accommodate 64bit architecture
#
# Revision 1.4  2007/08/31  05:11:06  avi
# linx and X86 Solaris compliant
#
# Revision 1.3  2007/01/18  03:51:09  avi
# new TRX
#
# Revision 1.2  2004/11/24  20:48:53  rsachs
# Added 'writeifh.c' to C source modules.
#
# Revision 1.1  2004/03/11  07:53:22  avi
# Initial revision
#

PROG	= T2C_4dfp
CSRCS	= ${PROG}.c
FSRCS	= reorder.f
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX	= ${NILSRC}/TRX
FLP	= ${NILSRC}/flip_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${FLP}/cflip.o

CFLAGS	= -I${TRX} -O
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

${PROG}: $(OBJS)
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

.c.o:
	${CC} -c $<
.f.o:
	${FC} -c $<

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co $(CSRCS) 

checkin:
	ci $(CSRCS) 
