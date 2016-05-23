noinst_PROGRAMS = \
		testDatetime \
		readcts

MCXXFLAGS = \
	-std=c++14 \
	-g \
	-pg \
	-Wall \
	-Wextra \
	-Werror \
	-pedantic \
	-W \
	-Wshadow \
	-Wdisabled-optimization \
	-DDEBUG \
	-Wno-unused-variable

testDatetime_SOURCES   = test_datetime.cpp
testDatetime_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
testDatetime_LDADD     = $(top_srcdir)/src/libdatetime.la

readcts_SOURCES   = read_cts.cpp
readcts_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
readcts_LDADD     = $(top_srcdir)/src/libdatetime.la
