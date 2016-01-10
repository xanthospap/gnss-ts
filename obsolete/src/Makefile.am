lib_LIBRARIES = libgnssts.a

## Only use the following flags for debuging purposes
libgnssts_a_CXXFLAGS = -std=c++14 
	-g \
	-Wall \
	-Wextra \
	-Werror \
	-pedantic \
	-W \
	-Wshadow \
	-Winline \
	-Wdisabled-optimization \
	-DDEBUG

dist_include_HEADERS = geodesy.hpp \
	geoconst.hpp \
	car2ell.hpp \
	car2top.hpp \
	ell2car.hpp \
	datetime.hpp \
	flags.hpp

dist_libgnssts_a_SOURCES = top2daz.cpp \
	datetime.cpp
