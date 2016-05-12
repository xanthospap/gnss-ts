#ifndef _CPPBASIC_HPP
#define _CPPBASIC_HPP

/* ************************************************************************** **
**                NATIONAL TECHNICAL UNIVERSITY OF ATHENS                     **
**                    DIONYSOS SATELLITE OBSERVATORY                          **
**                      HIGHER GEODESY LABORATORY                             **
** -------------------------------------------------------------------------- **
** HEADER FILE:                                                               **
** PURPOSE    :                                                               **
** DEP/CIES   :                                                               **
** COMP. FLAGS:                                                               **
** CREATED    :                                                               **
** MODIFIED   :                                                               **
** VERSION    :                                                               **
** TODO       :                                                               **
** REFERENCES :                                                               **
** ************************************************************************** */
/* FURTHER INFORMATION:                                                       */

#include <iostream>
#include <sstream>
#include <vector>

template<class T>
inline std::string numberToStdString(const T& num)
{
	std::ostringstream o;
	o << num;
	return o.str();
}
/* split a string to a vector of strings,
using char c as a */
std::vector<std::string> lineToStringVector(const std::string&, const char c=' ');
/* remove all instances of char c, both from the
begining and the end of a string */
std::string              strip             (const std::string&, const char c=' ');
/* convert a FORTRAN floating point string (eg 8934732D-5)
to a C-style floating point string (eg 8934732e-5) */
std::string              for2cppfloat      (const std::string&);

#endif
