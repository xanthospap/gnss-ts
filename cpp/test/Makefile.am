noinst_PROGRAMS = 

check_PROGRAMS = testTsFlag \
  testEarthquakeCatalogue \
  testIgsLogFile \
  testEventList \
  testTimeSeries1 \
  testTimeSeries2

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
	-DTS_DEBUG

AM_LDFLAGS = -lgeodesy -ldatetime

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

testTimeSeries1_SOURCES   = test_timeseries_1.cc
testTimeSeries1_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src ##-L$(top_srcdir)/src
testTimeSeries1_LDADD     = $(top_srcdir)/src/libgtms.la

testTimeSeries2_SOURCES   = test_timeseries_2.cc
testTimeSeries2_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src ##-L$(top_srcdir)/src
testTimeSeries2_LDADD     = $(top_srcdir)/src/libgtms.la

#lomb_scargle_SOURCES   = lomb_scargle.cpp
#lomb_scargle_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
#lomb_scargle_LDADD     = $(top_srcdir)/src/libgtms.la
#readcts_LDFLAGS   = -llibggeodesy.so -llibggdatetime.so
