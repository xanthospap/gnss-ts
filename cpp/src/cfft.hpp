#ifndef __NGPT_FFT_HPP__
#define __NGPT_FFT_HPP__

namespace ngpt
{

/// Discrete Fourier Transform (DFT).
void
four1__(double data[], std::size_t nn, int isign);

/// (Real) Fast Fourier Transform (FFT).
void
realft__(double data[], std::size_t n, int isign);

} // namespace ngpt

#endif
