// c++ standard headers
#include <stdexcept>
#include <cmath>
#ifdef DEBUG
#include <cassert>
#include <fenv.h>
#include <iostream>
#endif

// Fast Fourier Transform header
#include "cfft.hpp"

namespace ngpt
{

/// @brief Discrete Fourier Transform
///
/// Replaces data[0..2*nn-1] by its discrete Fourier transform, if isign is 
/// input as 1; or replaces data[0..2*nn-1] by nn times its inverse discrete 
/// Fourier transform, if isign is input as −1. data is a complex array of 
/// length nn or, equivalently, a real array of length 2*nn. nn MUST be an 
/// integer power of 2 (this is not checked for!).
///
/// @param[in,out] data  Array of size 2*nn; at input this is the data array. At
///                      output it is overwritten by its DFT (or inverse DFT
///                      depending on the value of sign).
/// @param[in]     nn    (Half) the size of data array.
/// @param[in]     sign  If +1 is given the DFT is computed; if it is -1 then
///                      the inverse DFT is computed instead.
///
/// @note This function is used within realft__.
///
/// Reference, Numerical Recipes in C, Ch. 12.2, pg. 507
void
four1__(double data[], std::size_t nn, int isign)
{
#ifdef DEBUG
    // enable catching of floating point exceptions in debug mode.
    feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
#endif

    std::size_t n,mmax,m,j,istep,i;
    double      wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

    n = nn << 1;
    j = 1;
    // This is the bit-reversal section of the function.
    for (i = 1; i < n; i += 2) {
        if (j > i ) {
            // Exchange the two complex numbers.
            std::swap(data[j-1], data[i-1]);
            std::swap(data[j], data[i]);
        }
        m = nn;
        while (m >= 2 && j > m ) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    // Here begins the Danielson-Lanczos section of the routine.
    mmax = 2;
    while (n > mmax) { // Outer loop executed log_2(nn) times.
        istep = mmax << 1;
        theta = isign*(6.28318530717959/mmax); // Initialize the trigonometric recurrence.
        wtemp = std::sin(0.5*theta);
        wpr   = -2.0*wtemp*wtemp;
        wpi   = std::sin(theta);
        wr    = 1e0;
        wi    = 0e0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j         = i+mmax;
                tempr     = wr*data[j-1]-wi*data[j];
                tempi     = wr*data[j]+wi*data[j-1];
                data[j-1] = data[i-1]-tempr;
                data[j]   = data[i]-tempi;
                data[i-1] += tempr;
                data[i]   += tempi;
            }
            wr = (wtemp=wr)*wpr-wi*wpi+wr;
            wi = wi*wpr+wtemp*wpi+wi;
        }
        mmax = istep;
    }

    // All done
    return;
}

/// @brief (Real) Fast Fourier Transform.
///
/// Calculates the Fourier transform of a set of n real-valued data points.
/// Replaces this data (which is stored in array data[0..n-1]) by the positive
/// frequency half of its complex Fourier transform. The real-valued first
/// and last components of the complex transform are returned as elements
/// data[0] and data[1], respectively. n must be a power of 2. This routine
/// also calculates the inverse transform of a complex data array if it is
/// the transform of real data. (Result in this case must be multiplied by 2/n.)
///
/// @param[in,out] data  Array of size n, containing at input the data values.
///                      At output it is overwritten by the its Fourier transform.
/// @param[in]     n     Size of data array.
/// @param[out]    sign  +1 for the forward transform, -1 for the inverse.
///
/// Reference: Numerical Recipes in C, Ch. 12.3
void
realft__(double data[], std::size_t n, int isign)
{
#ifdef DEBUG
    // enable catching of floating point exceptions in debug mode.
    feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
#endif

    std::size_t i,i1,i2,i3,i4,np3;
    double      c1=0.5e0,
                c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;

    theta = 3.141592653589793/(double)(n>>1); // Initialize the recurrence.

    if (isign == 1) {
        c2 = -.5e0;
        four1__(data, n>>1, 1);  // The forward transform is here.
    } else {
        c2    = .5e0;            // Otherwise set up for an inverse transform.
        theta = -theta;
    }

    wtemp = std::sin(0.5e0*theta);
    wpr   = -2e0*wtemp*wtemp;
    wpi   = std::sin(theta);
    wr    = 1e0+wpr;
    wi    = wpi;
    np3   = n+1;
    for (i = 1; i < (n>>2); i++) { // Case i=1 done separately below.
        i4  = 1+(i3=np3-(i2=1+(i1=i+i)));
        // The two separate transforms are separated out of data.
        h1r  =  c1*(data[i1]+data[i3]);
        h1i  =  c1*(data[i2]-data[i4]);
        h2r  = -c2*(data[i2]+data[i4]);
        h2i  =  c2*(data[i1]-data[i3]);
        // Here they are recombined to form the true transform of the original
        // true data
        data[i1] =  h1r+wr*h2r-wi*h2i;
        data[i2] =  h1i+wr*h2i+wi*h2r;
        data[i3] =  h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        // The recurrence.
        wr  = (wtemp=wr)*wpr-wi*wpi+wr;
        wi  = wi*wpr+wtemp*wpi+wi;
#ifdef DEBUG
        assert( std::min(std::min(i1, i2), std::min(i3, i4)) > 1 );
        assert( i1 < n && i2 < n && i3 < n && i4 < n );
#endif
    }

    if (isign == 1) {
        // Squeeze the first and last data together to get them all within the
        // original array.
        data[0] = (h1r=data[0])+data[1];
        data[1] = h1r-data[1];
    } else {
        // This is the inverse transform for the case isign=-1.
        data[0] = c1*((h1r=data[0])+data[1]);
        data[1] = c1*(h1r-data[1]);
        four1__(data,n>>1,-1);
    }
 
    // All done
    return;
}

} // namespace ngpt
