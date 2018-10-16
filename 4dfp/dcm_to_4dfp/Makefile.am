SUBDIRS = include

bin_PROGRAMS = dcm_to_4dfp

AM_CPPFLAGS = -I@srcdir@/include -I@TRX@

dcm_to_4dfp_SOURCES = dcm_to_4dfp.c config.h \
	pixel_ops.c lst.c dcm.c condition.c utility.c \
	ctnthread.c dcmsupport.c dcmdict.c dcmcond.c \
	gitversion.c

gitversion.c: @GITPREREQ@
	echo "const char *gitrev = \""`git rev-parse HEAD`"\";" >$@
	echo "const char *gittag = \""`git describe --always --tag`"\";" >>$@

dcm_to_4dfp_LDADD = @TRX@/rec.o @TRX@/Getifh.o @TRX@/endianio.o

dist-hook:
	@test -f $(top_srcdir)/generate-changelog.sh && sh $(top_srcdir)/generate-changelog.sh $(top_srcdir) $(distdir) || true
	tar cvfz $(PACKAGE).tar.gz --transform "s!$(distdir)!$(PACKAGE)!" \
		$(distdir)

DISTCLEANFILES=$(PACKAGE).tar.gz

TESTS = test_dcm_to_4dfp test_dcm
check_PROGRAMS = test_dcm_to_4dfp test_dcm

test_dcm_to_4dfp_SOURCES = $(dcm_to_4dfp_SOURCES) test_dcm_to_4dfp.c 
test_dcm_to_4dfp_CPPFLAGS = $(AM_CPPFLAGS) -DBUILD_UNIT_TEST
test_dcm_to_4dfp_LDADD = $(dcm_to_4dfp_LDADD)

test_dcm_SOURCES = dcm.c test_dcm.c
test_dcm_CPPFLAGS = $(AM_CPPFLAGS) -DBUILD_UNIT_TEST
test_dcm_LDADD = lst.o condition.o dcmdict.o dcmcond.o dcmsupport.o \
	@TRX@/endianio.o @TRX@/Getifh.o
