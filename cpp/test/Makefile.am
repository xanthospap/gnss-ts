noinst_PROGRAMS = \
		readcts \
		readcat

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

AM_LDFLAGS = -lggeodesy -lggdatetime

readcts_SOURCES   = read_cts.cpp
readcts_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
readcts_LDADD     = $(top_srcdir)/src/libgtms.la
#readcts_LDFLAGS   = -llibggeodesy.so -llibggdatetime.so

readcat_SOURCES   = test_earthqcat.cpp
readcat_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
readcat_LDADD     = $(top_srcdir)/src/libgtms.la
#readcat_LDFLAGS   = -llibggeodesy.so -llibggdatetime.so
