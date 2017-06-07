noinst_PROGRAMS = \
		readcts \
		readcat \
		period \
		lomb_scargle

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

lomb_scargle_SOURCES   = lomb_scargle.cpp
lomb_scargle_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
lomb_scargle_LDADD     = $(top_srcdir)/src/libgtms.la
#readcts_LDFLAGS   = -llibggeodesy.so -llibggdatetime.so

period_SOURCES   = test_periods.cpp
period_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
period_LDADD     = $(top_srcdir)/src/libgtms.la
#readcts_LDFLAGS   = -llibggeodesy.so -llibggdatetime.so

readcts_SOURCES   = read_cts.cpp
readcts_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
readcts_LDADD     = $(top_srcdir)/src/libgtms.la
#readcts_LDFLAGS   = -llibggeodesy.so -llibggdatetime.so

readcat_SOURCES   = test_earthqcat.cpp
readcat_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
readcat_LDADD     = $(top_srcdir)/src/libgtms.la
#readcat_LDFLAGS   = -llibggeodesy.so -llibggdatetime.so
