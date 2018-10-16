#$Id: dcm_dump_file.mak,v 1.1 2008/07/24 06:48:16 avi Exp $
#$Log: dcm_dump_file.mak,v $
# Revision 1.1  2008/07/24  06:48:16  avi
# Initial revision
#

PROG	= dcm_dump_file
CSRCS	= ${PROG}.c
LSRCS	= lst.c dcm.c condition.c utility.c ctnthread.c dcmsupport.c dcmdict.c dcmcond.c
OBJS    = ${CSRCS:.c=.o} ${LSRCS:.c=.o}

ARCH	= $(shell uname -p)
ifeq (${ARCH},sparc)
	ARCHITECTURE	= BIG_ENDIAN_ARCHITECTURE
else
	ARCHITECTURE	= LITTLE_ENDIAN_ARCHITECTURE
endif
LONGSIZE	= 32
INTSIZE		= 32
SHORTSIZE	= 16

C_OPTS	=  -g -DDEBUG -DSHARED_MEMORY \
-DSEMAPHORE -DX11 -DATHENA -DMOTIF -DX11R4 -DUSLEEP -DMSQL \
-DTBL_REQUIRES_HAT_ESCAPE -D$(ARCHITECTURE) \
-DLONGSIZE=$(LONGSIZE) -DINTSIZE=$(INTSIZE) \
-DSHORTSIZE=$(SHORTSIZE)

CCFLAGS		= -O -I./include ${C_OPTS}
CC		= cc ${CCFLAGS}

.c.o:
	$(CC) -c $<

${PROG}: ${OBJS} 
	${CC} -o $@ ${OBJS}

clean:
	/bin/rm ${PROG} ${OBJS}

test:
	echo ${ARCH}
	echo ${ARCHITECTURE}

release: ${PROG}
	chmod 751 ${PROG}
	/bin/mv ${PROG} ${RELEASE}
