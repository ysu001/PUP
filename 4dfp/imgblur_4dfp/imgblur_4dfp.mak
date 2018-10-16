#$Id: imgblur_4dfp.mak,v 1.6 2009/07/18 01:30:16 avi Exp $
#$Log: imgblur_4dfp.mak,v $
# Revision 1.6  2009/07/18  01:30:16  avi
# Solaris10/Linux compliant
#
# Revision 1.5  2004/09/20  21:38:40  rsachs
# Removed 'libmri'. Added 'get_4d_images2.o' & 'Getifh.o'.
#
# Revision 1.4  1999/01/29  06:30:20  avi
# LOBJS
#

PROG	= imgblur_4dfp
CSRCS	= ${PROG}.c
FSRCS	= fimgblur.f
TRX	= ${NILSRC}/TRX
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o

OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
LIBS	= -lm

CFLAGS	= -O -I${TRX}
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

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

${PROG}: $(OBJS) 
	${FC} -o $@ ${LOBJS} ${OBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG} 

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}



