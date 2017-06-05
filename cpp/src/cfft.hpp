#ifndef __NGPT_FFT_HPP__
#define __NGPT_FFT_HPP__

#include <stdlib.h>

namespace ngpt
{

void
four1__(double data[], std::size_t nn, int isign);

void
realft__(double data[], std::size_t n, int isign);

} // namespace ngpt

#endif
