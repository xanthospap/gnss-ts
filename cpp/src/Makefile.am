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
        cfft.hpp \
	genflags.hpp \
	tsflag.hpp \
	ts_flag.hpp \
	timeseries.hpp \
	period.hpp \
	cts_read.hpp \
	earthquake_cat.hpp \
	earthquake.hpp \
	earthquake_catalogue.hpp \
	event_list.hpp \
	psd.hpp \
	model.hpp \
	crdts.hpp \
	artificial.hpp \
	stalogrd.hpp \
        igs_sta_log.hpp \
	procalgo.hpp

dist_libgtms_la_SOURCES = \
	tsflag.cpp \
        cfft.cpp \
	earthquake.cpp \
	earthquake_catalogue.cpp \
        igs_sta_log.cpp
