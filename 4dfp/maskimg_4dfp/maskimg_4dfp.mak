#$Header: /data/petsun4/data1/src_solaris/maskimg_4dfp/RCS/maskimg_4dfp.mak,v 1.10 2006/09/25 02:32:47 avi Exp $
#$Log: maskimg_4dfp.mak,v $
# Revision 1.10  2006/09/25  02:32:47  avi
# ${PROG} ${RELEASE}
#
# Revision 1.9  2006/03/23  04:59:26  avi
# variaable endian subroutines in new TRX
#
# Revision 1.8  2006/03/16  07:53:05  avi
# all #include local
#
# Revision 1.7  2006/02/22  04:30:05  avi
# point to local TRX
# properly define LOBJS
#
# Revision 1.6  2004/09/17  20:18:04  rsachs
# Moved 'get_4d_images2' to the list of objects.
#
# Revision 1.5  2004/09/17  20:08:25  rsachs
# Removed the 'mri' library.
#
# Revision 1.4  2004/09/16  20:17:49  rsachs
# Added 'get_4dfp_images2.c' & 'Getifh.o'.
#
# Revision 1.3  2003/10/05  21:14:59  avi
# Revision 1.2  1998/12/19  07:04:58  avi
# new release
#
# Revision 1.1  1998/10/12  20:36:20  mcavoy
# Initial revision
#
PROG	= maskimg_4dfp
TRX     = ${NILSRC}/TRX
CSRCS	= ${PROG}.c 
OBJS	= ${CSRCS:.c=.o}
LOBJS	= ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o
LIBS	= -lm

CCFLAGS = -O -I${TRX}
CC	= cc ${CCFLAGS}

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

.c.o:
	${CC} -c $<

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

checkout:
	co $(CSRCS)

checkin:
	ci $(CSRCS)
