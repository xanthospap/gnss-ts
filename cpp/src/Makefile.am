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
	-DTS_DEBUG

dist_include_HEADERS = \
	genflags.hpp \
	ts_flag.hpp \
	earthquake.hpp \
	earthquake_catalogue.hpp \
        igs_sta_log.hpp \
	event_list.hpp \
        data_point.hpp \
        timeseries.hpp \
        md_jump.hpp \
        md_harmonics.hpp \
        md_velocity_change.hpp \
        model.hpp \
	kalman.hpp \
	cts_read.hpp

dist_libgtms_la_SOURCES = \
	earthquake.cpp \
	earthquake_catalogue.cpp \
        igs_sta_log.cpp \
	event_list.cpp \
        timeseries.cpp \
        model.cpp \
	kalman.cpp
