#$Id: t4imgs_4dfp.mak,v 1.15 2009/02/25 21:29:27 avi Exp $
#$Log: t4imgs_4dfp.mak,v $
# Revision 1.15  2009/02/25  21:29:27  avi
# accommodate 64bit architecture
#
# Revision 1.14  2007/08/10  02:57:44  avi
# conditionally compile on linus using either -lgfortran or -lg2c
#
# Revision 1.13  2007/04/29  05:20:22  avi
# gcc v3 compliant: eliminate -librms; add set_rnan.f
#
# Revision 1.12  2006/09/26  23:32:42  avi
# correct cc options
#
# Revision 1.11  2006/09/26  23:26:31  avi
# ${PROG} ${RELEASE}
#
# Revision 1.10  2006/02/13  07:33:33  avi
# eliminate dependency on TEC and libmri
#
# Revision 1.9  2004/09/03  19:54:26  rsachs
#
# Revision 1.8  2004/08/23  23:58:23  rsachs
# Added ft4imgn.f
#
# Revision 1.7  2001/01/07  05:03:38  avi
# link with object modules needed for 3D cubic spline interpolation
#
# Revision 1.6  1999/01/10  04:08:45  avi
# t4scale included in LOBJS
#
# Revision 1.5  1998/12/29  07:51:23  avi
# use ${NILSRC}/imglin/t4_sub.o instead of t4_read.o
#
# Revision 1.4  1998/12/25  08:42:19  avi
# Revision 1.3  1998/06/28  22:06:32  avi
# rec.c
#
# Revision 1.2  1997/10/11  07:04:58  avi
# include to_711-2B.f
#
# Revision 1.1  1997/10/11  00:18:46  avi
# Initial revision
#
#	Makefile:	t4imgs_4dfp (solaris)
#	Authors:	Avi Snyder
#	Date:		11-Aug-97
#
PROG	= t4imgs_4dfp
CSRCS	= ${PROG}.c
FSRCS	= ft4imgo.f to_711-2B.f ft4imgn.f set_rnan.f
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

TRX	= ${NILSRC}/TRX
RMS	= ${NILSRC}/librms
T4A	= ${NILSRC}/t4_actmapf_4dfp
LIN	= ${NILSRC}/imglin
REG	= ${NILSRC}/imgreg_4dfp
FLP	= ${NILSRC}/flip_4dfp
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o ${FLP}/cflip.o ${TRX}/spline3dvgh.o \
	  ${RMS}/fftsol.o ${RMS}/imgpad.o ${RMS}/polfit.o ${RMS}/matopr.o ${RMS}/param6opr.o \
	  ${REG}/imgvalx.o ${LIN}/t4_sub.o ${LIN}/t4scale.o ${T4A}/ft4ixyz.o
LIBS	= -lm

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

.c.o:
	${CC} -c $<

.f.o:
	${FC} -c $<

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS} ${FSRCS}

checkin:
	ci ${CSRCS} ${FSRCS}

