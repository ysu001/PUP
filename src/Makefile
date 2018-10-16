PROG	= calrsfmat PreprocessRSFtest2 rsfroieval2 roieval2 roieval3 fs2brain loganREFROI
CSRCS	= nrutil.c Getifh.c endianio.c cgauss3d.c RSF.c loganREF.c
FSRCs	= fftsol.f
OBJS	= ${CSRCS:.c=.o} fftsol.o
CFLAGS	= -O -I. -std=c99
CC	= gcc ${CFLAGS}
FC	= gcc  -O -ffixed-line-length-132 -fno-second-underscore

LIBS	= -lm -lgfortran


.c.o:
	${CC} -c $<
.f.o:
	${FC} -c $<
	
calrsfmat: calrsfmat.c ${OBJS}
	${CC} -o $@ calrsfmat.c ${OBJS} ${LIBS}

PreprocessRSFtest2: PreprocessRSFtest2.c ${OBJS}
	${CC} -o $@ PreprocessRSFtest2.c ${OBJS} ${LIBS}		

rsfroieval2: rsfroieval2.c ${OBJS}
	${CC} -o $@ rsfroieval2.c ${OBJS} ${LIBS}

roieval2: roieval2.c ${OBJS}
	${CC} -o $@ roieval2.c ${OBJS} ${LIBS}
	
roieval3: roieval3.c ${OBJS}
	${CC} -o $@ roieval3.c ${OBJS} ${LIBS}

fs2brain: fs2brain.c ${OBJS}
	${CC} -o $@ fs2brain.c ${OBJS} ${LIBS}

loganREFROI: loganREFROI.c ${OBJS}
	${CC} -o $@ loganREFROI.c ${OBJS} ${LIBS}

ALL: ${PROG}.c ${OBJS}
	${CC} -o calrsfmat calrsfmat.c ${OBJS} ${LIBS}
	${CC} -o PreprocessRSFtest2 PreprocessRSFtest2.c ${OBJS} ${LIBS}		
	${CC} -o rsfroieval2 rsfroieval2.c ${OBJS} ${LIBS}
	${CC} -o roieval2 roieval2.c ${OBJS} ${LIBS}
	${CC} -o roieval3 roieval3.c ${OBJS} ${LIBS}
	${CC} -o fs2brain fs2brain.c ${OBJS} ${LIBS}
	${CC} -o loganREFROI loganREFROI.c ${OBJS} ${LIBS}

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	rm ${OBJS} ${PROG}

