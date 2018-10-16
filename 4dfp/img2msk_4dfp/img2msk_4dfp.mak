#$Id: img2msk_4dfp.mak,v 1.6 2009/02/25 21:29:26 avi Exp $
#$Log: img2msk_4dfp.mak,v $
# Revision 1.6  2009/02/25  21:29:26  avi
# accommodate 64bit architecture
#
# Revision 1.5  2007/09/21  22:22:59  avi
# linux v3 and v4 compliant
#
# Revision 1.4  2007/03/02  07:38:10  avi
# Solaris 10
#
# Revision 1.3  2004/09/16  21:11:40  rsachs
# Added 'get_4d_images2.c' & 'Getifh.o'.
#
# Revision 1.2  1998/12/21  23:03:47  avi
#

PROG 	= img2msk_4dfp
TRX	= ${NILSRC}/TRX
CSRCS	= ${PROG}.c
FSRCS	= fimg2msk.f
OBJS 	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

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

${PROG}: ${OBJS}
	${FC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS}

checkout:
	co $(CSRCS) $(FSRCS) 

checkin:
	ci $(CSRCS) $(FSRCS) 
