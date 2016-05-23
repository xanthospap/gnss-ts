lib_LTLIBRARIES = libdatetime.la

## Only use the following flags for debuging purposes
libdatetime_la_CXXFLAGS = \
	-std=c++14 \
	-g \
	-pg \
	-Wall \
	-Wextra \
	-Werror \
	-pedantic \
	-W \
	-Wshadow \
	-Winline \
	-Wdisabled-optimization \
	-DDEBUG

dist_include_HEADERS = \
	geoconst.hpp \
	dtfund.hpp \
	dtcalendar.hpp \
	datetime_read.hpp \
	ellipsoid.hpp \
	geodesy.hpp \
	car2ell.hpp \
	car2top.hpp \
	ell2car.hpp \
	genflags.hpp \
	tsflagenum.hpp \
	timeseries.hpp \
	crdts.hpp \
	cts_read.hpp

dist_libdatetime_la_SOURCES = \
	dtfund.cpp \
	top2daz.cpp
