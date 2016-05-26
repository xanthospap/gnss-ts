noinst_PROGRAMS = \
		testDatetime \
		readcts \
		readcat \
		vincenty

MCXXFLAGS = \
	-std=c++14 \
	-g \
	-pg \
	-Wno-unused-variable \
	-Wall \
	-Wextra \
	-Werror \
	-pedantic \
	-W \
	-Wshadow \
	-Wdisabled-optimization \
	-DDEBUG

testDatetime_SOURCES   = test_datetime.cpp
testDatetime_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
testDatetime_LDADD     = $(top_srcdir)/src/libdatetime.la

readcts_SOURCES   = read_cts.cpp
readcts_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
readcts_LDADD     = $(top_srcdir)/src/libdatetime.la

readcat_SOURCES   = test_earthqcat.cpp
readcat_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
readcat_LDADD     = $(top_srcdir)/src/libdatetime.la

vincenty_SOURCES   = test_vincenty.cpp
vincenty_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
vincenty_LDADD     = $(top_srcdir)/src/libdatetime.la
