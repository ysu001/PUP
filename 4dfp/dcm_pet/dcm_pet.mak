# $Header: dcm_pet.mak,v 1.4 2018/10/16 10:13:00 ysu001 Exp $
#	Makefile for dcm_pet
# $Log: dcm_pet.mak,v $
# Revision 1.4 2018/10/16 10:13:00 ysu001
# Added release target to allow installation to $RELEASE folder
#
# Revision 1.3  2011/06/07 21:57:08  jon
# linux version
#


PROG	= dcm_pet
CSRCS	= ${PROG}.c

TRX = ${NILSRC}/TRX
LSRCS	= pixel_ops.c lst.c dcm.c condition.c utility.c \
ctnthread.c dcmsupport.c dcmdict.c dcmcond.c sequences.c \
sqcond.c
 
OBJS    = ${CSRCS:.c=.o} ${LSRCS:.c=.o}
LOBJS   = ${TRX}/rec.o ${TRX}/Getifh.o ${TRX}/endianio.o
LIBS	= -lm

ARCH	= $(shell uname -p)
ifeq (${ARCH},sparc)
	ARCHITECTURE	= BIG_ENDIAN_ARCHITECTURE
else
	ARCHITECTURE	= LITTLE_ENDIAN_ARCHITECTURE
endif
LONGSIZE	= 32
INTSIZE		= 32
SHORTSIZE	= 16

C_OPTS	= -DDEBUG -DMALLOC_DEBUG -DSHARED_MEMORY \
-DSEMAPHORE -DX11 -DATHENA -DMOTIF -DX11R4 -DUSLEEP -DMSQL \
-DTBL_REQUIRES_HAT_ESCAPE -D$(ARCHITECTURE) -DLONGSIZE=$(LONGSIZE) \
-DINTSIZE=$(INTSIZE) -DSHORTSIZE=$(SHORTSIZE)

CFLAGS	= -O -Iinclude -I/usr/include -I${TRX} ${C_OPTS} -lm

CC = gcc ${CFLAGS}

${PROG}: ${OBJS} 
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS}

release: ${PROG}
	chmod 771 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

.c.o:
	$(CC) -c $<


