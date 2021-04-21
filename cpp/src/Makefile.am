lib_LTLIBRARIES = libgtms.la

##
## Only use the following flags for debuging purposes
## --------------------------------------------------
##
libgtms_la_CXXFLAGS = \
	-std=c++17 \
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
	genflags.hpp \
	ts_flag.hpp \
	earthquake.hpp \
	earthquake_catalogue.hpp \
        igs_sta_log.hpp \
	event_list.hpp \
        timeseries.hpp

dist_libgtms_la_SOURCES = \
	earthquake.cpp \
	earthquake_catalogue.cpp \
        igs_sta_log.cpp \
	event_list.cpp \
        timeseries.cpp
