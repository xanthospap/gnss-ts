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

readcts_SOURCES   = read_cts.cpp
readcts_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
readcts_LDADD     = $(top_srcdir)/src/libgtms.la

readcat_SOURCES   = test_earthqcat.cpp
readcat_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
readcat_LDADD     = $(top_srcdir)/src/libgtms.la
