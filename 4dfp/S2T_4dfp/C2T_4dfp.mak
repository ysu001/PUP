#$Id: C2T_4dfp.mak,v 1.6 2009/02/25 21:29:24 avi Exp $
#$Log: C2T_4dfp.mak,v $
# Revision 1.6  2009/02/25  21:29:24  avi
# accommodate 64bit architecture
#
# Revision 1.5  2007/08/31  04:48:58  avi
# linux and X86 Solaris compliant
#
# Revision 1.4  2007/01/18  03:25:05  avi
# new TRX pointers
#
# Revision 1.3  2004/12/03  22:06:08  rsachs
# Added ${PROG} to the 'release' fragment.
#
# Revision 1.2  2004/11/23  20:52:14  rsachs
# Added 'writeifh.o' to list of object modules.
#
# Revision 1.1  2004/03/11  05:45:46  avi
# Initial revision
#

PROG	= C2T_4dfp
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
	$(CC) -c $<
.f.o:
	$(FC) -c $<

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co $(CSRCS) 

checkin:
	ci $(CSRCS) 
