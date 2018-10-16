# $Id: img_hist_4dfp.mak,v 1.9 2006/09/23 21:53:44 avi Exp $
# $Log: img_hist_4dfp.mak,v $
# Revision 1.9  2006/09/23  21:53:44  avi
# correct disasterous rscahs OBJS legacy
#
# Revision 1.8  2006/09/23  21:25:09  avi
# ${NILSRC} ${RELEASE}
#
# Revision 1.7  2006/06/28  01:01:22  avi
# new TRX and endianio.o links
#
# Revision 1.6  2004/09/17  21:16:04  rsachs
# Removed -lmri. Installed 'Getifh.o', 'get_4dfp_dimo'.
#
# Revision 1.5  2003/05/19  23:40:47  avi
#
# Revision 1.4  1999/07/10  02:21:41  avi
# cleaned
#
# Revision 1.3  1999/07/10  02:13:38  avi
# Revision 1.2  1998/10/09  17:57:16  mcavoy

PROG	= img_hist_4dfp
CSRCS	= ${PROG}.c 
FSRCS	=
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}
TRX     = ${NILSRC}/TRX
LOBJS	= ${TRX}/Getifh.o ${TRX}/endianio.o
LIBS	= -lm

CCFLAGS	= -O -I${TRX}
CC	= cc ${CCFLAGS}

.c.o:
	$(CC) -c $<

${PROG}: $(OBJS)
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS} 

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE} 

clean:
	rm ${OBJS} ${PROG}

checkout:
	co $(CSRCS) 

checkin:
	ci $(CSRCS) 
