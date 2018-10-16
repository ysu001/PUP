#$Header: /data/petsun4/data1/src_solaris/imgmax_4dfp/RCS/imgmax_4dfp.mak,v 1.6 2006/09/24 01:30:10 avi Exp $
#$Log: imgmax_4dfp.mak,v $
#Revision 1.6  2006/09/24 01:30:10  avi
#${PROG} ${RELEASE}
#

PROG 	= imgmax_4dfp
CSRCS	= ${PROG}.c 
TRX	= ${NILSRC}/TRX
OBJS	= ${CSRCS:.c=.o}
LOBJS	= ${TRX}/Getifh.o ${TRX}/endianio.o
LIBS	= -lm

CCFLAGS	= -O -I${TRX}
CC	= cc ${CCFLAGS}

${PROG}: $(OBJS)
	$(CC) -o $@ $(OBJS) ${LOBJS} ${LIBS} 

.c.o:
	$(CC) -c $<

release: ${PROG}
	chmod 771 ${PROG} 
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS}

checkout:
	co $(CSRCS) 

checkin:
	ci $(CSRCS) 
