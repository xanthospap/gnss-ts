lib_LTLIBRARIES = libgtms.la

##
## Only use the following flags for debuging purposes
## --------------------------------------------------
##
libgtms_la_CXXFLAGS = \
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
	-march=native \
	-DDEBUG

dist_include_HEADERS = \
        cfft.hpp \
	genflags.hpp \
	timeseries.hpp \
	period.hpp \
	cts_read.hpp \
	earthquake_cat.hpp \
	event_list.hpp \
	model.hpp \
	crdts.hpp \
	artificial.hpp

dist_libgtms_la_SOURCES = \
        cfft.cpp
