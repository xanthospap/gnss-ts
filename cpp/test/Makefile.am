noinst_PROGRAMS = 

check_PROGRAMS = testTsFlag \
  testEarthquakeCatalogue \
  testIgsLogFile \
  testEventList

MCXXFLAGS = \
	-std=c++17 \
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

#test_ts_SOURCES   = test_ts.cpp
#test_ts_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
#test_ts_LDADD     = $(top_srcdir)/src/libgtms.la

testTsFlag_SOURCES   = test_ts_flag.cc
testTsFlag_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src ##-L$(top_srcdir)/src
testTsFlag_LDADD     = $(top_srcdir)/src/libgtms.la

testEarthquakeCatalogue_SOURCES   = test_earthquake_catalogue.cc
testEarthquakeCatalogue_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src ##-L$(top_srcdir)/src
testEarthquakeCatalogue_LDADD     = $(top_srcdir)/src/libgtms.la

testIgsLogFile_SOURCES   = test_igs_log_file.cc
testIgsLogFile_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src ##-L$(top_srcdir)/src
testIgsLogFile_LDADD     = $(top_srcdir)/src/libgtms.la

testEventList_SOURCES   = test_event_list.cc
testEventList_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src ##-L$(top_srcdir)/src
testEventList_LDADD     = $(top_srcdir)/src/libgtms.la

#lomb_scargle_SOURCES   = lomb_scargle.cpp
#lomb_scargle_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
#lomb_scargle_LDADD     = $(top_srcdir)/src/libgtms.la
#readcts_LDFLAGS   = -llibggeodesy.so -llibggdatetime.so
