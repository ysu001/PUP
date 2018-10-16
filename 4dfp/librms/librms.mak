#$Header: /data/petsun4/data1/src_solaris/librms/RCS/librms.mak,v 1.5 2011/04/28 21:16:40 avi Exp $
#$Log: librms.mak,v $
# Revision 1.5  2011/04/28  21:16:40  avi
# dsvdcmp0.
#
# Revision 1.4  2009/01/20  01:00:22  avi
# dmatopr.f and dgeigen.c
#
# Revision 1.3  2007/09/21  07:22:02  avi
# -I${NILSRC}/TRX in ${CFLAGS}
#
# Revision 1.2  2007/09/21  06:13:10  avi
# deigen.f dmatinv.f
#
# Revision 1.1  2007/09/18  23:24:46  avi
# Initial revision
#

CSRCS	= cgauss3d.c imag2mask.c img2lmask.c dgeigen.c dsvdcmp0.c
FSRCS	= deconv.f deigen.f determ12.f dmatinv.f eigen.f fftsol.f frmsmri.f \
	  gauss3d.f hipass3d.f imgpad.f matopr.f dmatopr.f npad.f param6opr.f polfit.f deigen.f dmatinv.f
OBJS  = ${FSRCS:.f=.o} ${CSRCS:.c=.o} 

CFLAGS	= -O -I. -I${NILSRC}/TRX
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	FC	= gcc -O -ffixed-line-length-132 -fcray-pointer
else
	CC	= cc ${CFLAGS}
	FC	= f77 -O -I4 -e
endif

.c.o:
	$(CC) -c $<

.f.o:
	$(FC) -c $<

RMS: ${OBJS}
	echo done

clean: 
	rm $(OBJS)

checkout: 
	co $(CSRCS) $(FSRCS)
