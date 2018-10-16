#$Id: imgreg_4dfp.mak,v 1.10 2009/02/25 21:29:27 avi Exp $
#$Log: imgreg_4dfp.mak,v $
# Revision 1.10  2009/02/25  21:29:27  avi
# accommodate 64bit architecture
#
# Revision 1.9  2007/08/10  03:28:14  avi
# conditional gcc FORTRAN library finding
#
# Revision 1.8  2007/04/25  05:13:54  avi
# linux distribution compliant
# eliminate -lrms
#
# Revision 1.7  2006/09/26  17:56:15  avi
# ${PROG} ${RELEASE}
#
# Revision 1.6  2004/09/24  19:34:07  rsachs
# Removed 'libmri'. Added
#
# Revision 1.5  1998/12/26  04:37:19  avi
# t4_sub.o and aviparse.o taken from /data/petsun4/data1/src_solaris/imglin
#
# Revision 1.4  1998/12/24  07:29:28  avi
# LOBJS
#
# Revision 1.3  1998/12/14  07:50:55  avi
# use objects in TRX
#
# Revision 1.2  1998/10/12  20:08:50  mcavoy
#

PROG	= imgreg_4dfp
CSRCS	= ${PROG}.c
FSRCS	= fimgreg.f imgvalm.f imgvalx.f ffind.f
TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
LIN	= ${NILSRC}/imglin
FLP	= ${NILSRC}/flip_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${FLP}/cflip.o ${LIN}/t4_sub.o \
	  ${RMS}/polfit.o ${RMS}/matopr.o ${RMS}/eigen.o ${RMS}/param6opr.o
OBJS	= ${FSRCS:.f=.o} ${CSRCS:.c=.o}
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
