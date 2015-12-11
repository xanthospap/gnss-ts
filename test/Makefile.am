bin_PROGRAMS = testGeodesy testFlags

testGeodesy_SOURCES  = test_geodesy.cpp
testGeodesy_CXXFLAGS = -std=c++14 -g -Wall -Wextra -pedantic -W -Wshadow \
	-Winline -Wdisabled-optimization -DDEBUG \
	-I$(top_srcdir)/src -L$(top_srcdir)/src
testGeodesy_LDADD    = $(top_srcdir)/src/libgnssts.a

testFlags_SOURCES  = test_flags.cpp
testFlags_CXXFLAGS = -std=c++14 -g -Wall -Wextra -pedantic -W -Wshadow \
	-Winline -Wdisabled-optimization -DDEBUG \
	-I$(top_srcdir)/src -L$(top_srcdir)/src
testFlags_LDADD    = $(top_srcdir)/src/libgnssts.a
