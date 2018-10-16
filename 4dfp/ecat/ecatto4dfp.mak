#$Header: /data/petsun4/data1/src_solaris/ecat/RCS/ecatto4dfp.mak,v 1.1 2010/09/01 17:38:34 larsc Exp $
#$Log: ecatto4dfp.mak,v $
# Revision 1.1  2010/09/01  17:38:34  larsc
# Initial revision
#
PROG	= ecatto4dfp
CSRCS	= ${PROG}.c
TRX	= ../TRX
ECAT	= libecat7-1.5/src
ECATUTIL	= libecat7-1.5/utils
LOBJS	= ${ECAT}/analyze.o ${ECAT}/convert_64.o ${ECAT}/convert_70.o ${ECAT}/crash.o \
	  ${ECAT}/interfile.o ${ECAT}/load_volume7.o ${ECAT}/machine_indep.o \
	  ${ECAT}/matrix.o ${ECAT}/matrix_64.o ${ECAT}/matrix_extra.o \
	  ${ECAT}/matrix_slice.o ${ECAT}/matrix_xdr.o ${ECAT}/num_sort.o \
	  ${ECAT}/rfa_xdr.o ${ECAT}/rts_cmd.o ${ECAT}/save_volume7.o \
	  ${ECAT}/sino_dets.o ${TRX}/endianio.o ${TRX}/Getifh.o ${TRX}/Inithdr.o \
	  ${TRX}/rec.o
OBJS	= ${CSRCS:.c=.o}

.c.o:
	${CC} -c $<

CFLAGS	= -O -I. -I${ECAT} -I${ECATUTIL} -I${TRX}
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
