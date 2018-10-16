#$Header: /data/petsun4/data1/src_solaris/ecat/RCS/ecat_header.mak,v 1.2 2010/09/08 20:52:27 avi Exp $
#$Log: ecat_header.mak,v $
# Revision 1.2  2010/09/08  20:52:27  avi
# remove superfluous (and evidently problematic) compile directives
#
# Revision 1.1  2010/09/01  17:38:01  larsc
# Initial revision
#
PROG	= ecat_header
CSRCS	= ${PROG}.c
ECAT	= libecat7-1.5/src
LOBJS	= ${ECAT}/analyze.o ${ECAT}/convert_64.o ${ECAT}/convert_70.o ${ECAT}/crash.o \
	  ${ECAT}/interfile.o ${ECAT}/load_volume7.o ${ECAT}/machine_indep.o \
	  ${ECAT}/matrix.o ${ECAT}/matrix_64.o ${ECAT}/matrix_extra.o \
	  ${ECAT}/matrix_slice.o ${ECAT}/matrix_xdr.o ${ECAT}/num_sort.o \
	  ${ECAT}/rfa_xdr.o ${ECAT}/rts_cmd.o ${ECAT}/save_volume7.o \
	  ${ECAT}/sino_dets.o
	  
OBJS	= ${CSRCS:.c=.o}

.c.o:
	${CC} -c $<

CFLAGS	= -I. -I${ECAT} -O
ifeq (${OSTYPE}, linux)
	CC	= gcc ${CFLAGS}
	LIBS	= -lm
else
	CC	= cc ${CFLAGS}
	LIBS	= -lm -lrpcsvc -lnsl
endif

${PROG}: ${OBJS}
	${CC} -o $@ ${OBJS} ${LOBJS} ${LIBS} 

release: ${PROG}
	chmod 775 ${PROG}
	/bin/mv ${PROG} ${RELEASE}

clean:
	/bin/rm ${OBJS} ${PROG}

checkout:
	co ${CSRCS}

checkin:
	ci ${CSRCS}
