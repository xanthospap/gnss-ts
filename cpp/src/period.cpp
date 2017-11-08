#include "period.hpp"

/// @brief Extirpolation
///
/// Given an array yy[0..n-1], extirpolate (spread) a value y into m actual
/// array elements that best approximate the "fictional" (i.e. possibly 
/// non-integer) array element number x. The weights used are coefficients of 
/// the Lagrange interpolating polynomial.
/// This function is used for the fast Lomb-Scargle algorithm (i.e. via FFT).
/// If (the input) is an integer, then the function will set yy[x] to yy[x]+y.
/// Else, it will fill the entries yy[x-m/2+1 : x-m/2-1].
///
/// @param[in] y  The value to "spread"
/// @param[in] yy The array to extirpolate at
/// @param[in] n  The array size
/// @param[in] m  To be used for Lagrange polynomial interpolation
///
/// Reference: Numerical Recipes in C, Ch. 13.8, pg. 583
void
ngpt::spread__(double y, double yy[], std::size_t n, double x, int m)
{

/*#ifdef DEBUG
    // enable catching of floating point exceptions in debug mode.
    feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
#endif*/

    long        ihi,ilo,ix,j,nden;
    static long nfac[] = {0,1,1,2,6,24,120,720,5040,40320,362880};
    double      fac,intpart;

    if (m > 10) {
        throw std::runtime_error
            {"spread__: factorial table too small in spread"};
    }

    ix = static_cast<long>(x);            // value x as index
    if (std::modf(x, &intpart) == 0e0 ) { // if x is integer ....
        yy[ix] += y;
    } else {                              // x is not an integer ...
        // lowest index to fill
        ilo  = std::min(std::max(static_cast<long>(x-.5e0*m+1e0), 0L),
                        (long)(n-m));
        // highest index to fill
        ihi  = ilo + m - 1;
#ifdef DEBUG
        assert( x   >= 0 );
        assert( ihi <  (long)n && ihi >= 0 );
        assert( ilo >= 0 );
#endif
        nden = nfac[m];
        fac  = x-ilo;
        for (j = ilo+1; j <= ihi; j++) fac *= (x-j);
        yy[ihi] += y*fac/(nden*(x-ihi));
        for (j = ihi-1; j >= ilo; j--) {
            nden   = (nden/(j+1-ilo))*(j-ihi);
            yy[j] += y*fac/(nden*(x-j));
        }
    }

    // all done
    return;
}
