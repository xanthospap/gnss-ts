#include <numeric>
#include <algorithm>
#include <cassert>

/*
 * Find the bit (i.e. power of 2) immidiately greater than of equal to N.
 * @warning This works for N in range [0, 2**64]
 */
int
bitceil(int N)
{
    N--;
    N |= N >> 1;
    N |= N >> 2;
    N |= N >> 4;
    N |= N >> 8;
    N |= N >> 16;
    N |= N >> 32;
    N++;
    return N;
}

void
extirpolate(const double* x, const double* y, std::size_t sz, int N=0, int M=4)
{
    static long nfac[/*11*/] = {0,1,1,2,6,24,120,720,5040,40320,362880};

    if ( !N ) {
        N = static_cast<int>(x[sz-1] + .5e0 * M + 1);
    }
    
    // Now use legendre polynomial weights to populate the results array;
    // This is an efficient recursive implementation (See Press et al. 1989)
    int low {0}, high{N-M};
    for (std::size_t i = 0; i< sz; i++) {
        if (static_cast<double>(std::floor(x[i])) == x[i]) {
            result[i] = y[i];
        } else {
            long tmp { x[i]-(M/2) };
            if (tmp < low) {
                tmp = low
            } else if (tmp > high) {
                tmp = high;
            }
            result[i] = tmp;
        }
    }


}

/*
 * Compute (approximate) trigonometric sums for a number of frequencies
 * 
 * @details This routine computes weighted sine and cosine sums:
 *
 * S_j = sum_i { h_i * sin(2 pi *f_j * t_i) }
 * C_j = sum_i { h_i * cos(2 pi *f_j * t_i) }
 *
 * Where f_j = freq_factor * (f0 + j * df) for values j in 1 ... N.
 * The sums can be computed either by a brute force O[N^2] method, or
 * by an FFT-based O[Nlog(N)] method].
 *
 * @param[in] t            Array of input times (of size M)
 * @param[in] h            Array of weights for the sum (of size M)
 * @param[in] df           Frequency spacing
 * @param[in] M            Number of entries in input arrays (t and h)
 * @param[in] N            Number of frequency bins to return
 * @param[in] f0           The low frequency to use (default = 0)
 * @param[in] freq_factor  Factor which multiplies the frequency (default = 1)
 * @param[in] use_fft      If true, use the approximate FFT algorithm to compute
 *                         the result. This uses the Press & Rybicki's
 *                         Lagrangian extirpolation.
 * @param[in] oversampling Oversampling freq_factor for the approximation;
 *                         roughly the number of time samples across the 
 *                         highest-frequency sinusoid. This parameter contains
 *                         the tradeoff between accuracy and speed. Not
 *                         referenced if use_fft is false (default=5).
 * @param[in] Mfft         The number of adjacent points to use in the FFT
 *                         approximation. Not referenced if use_fft is false
 *                         (default = 4)
 * @param[out] S           Array of size N
 * @param[out] C           Array of size N
 *
 * Reference: astroML, 
 */
void
trig_sum(const double* t, const double* h, double df, std::size_t M, std::size_t N,
    double* S, double* C,
    double f0=0e0, double freq_factor=1e0, int oversampling=5,
    bool use_fft=true, int Mfft=4)
{
    df *= freq_factor;
    f0 *= freq_factor;
    assert(df > 0);

    if ( !use_fft ) {
        double f_j;
        for (std::size_t j = 0; j < N; j++) {
            f_j  = freq_factor * (f0 + j*df);
            S[i] = 0e0;
            C[i] = 0e0;
            for (std::size_t i = 0; i < M; i++) {
                S[i] += h[i] * std::sin(2e0*M_PI*f_j*t[i]);
                C[i] += h[i] * std::cos(2e0*M_PI*f_j*t[i]);
            }
        }
    } else {
        assert( Mfft > 0 );
        // required size of fft is the power of 2 above the oversampling rate
        int fft_pow2 = static_cast<int>(N*oversampling);
        assert( fft_pow2 > 0 && fft_pow2 < std::pow(2, 64) );
        int Nfft = bitceil(fft_pow2);
        auto t0  = t[0];
        if (f0 > 0):

    }

    return;
}

void
lscargle(double* t, double* y, double* dy=nullptr, std::size_t N,
    double f0=0, double df=-1, int Nf=-1,
    bool center_data=true, bool fit_offset=true,
    bool use_fft=true, int freq_oversampling=5,
    int nyquist_factor=2)
{

    // setup data weights
    double* w = new double[N];
    if (dy) {
        double wsum = 0;
        for (std::size_t i = 0; i< N; i++) {
            w[i] = 1e0 / (dy[i]*dy[i]);
            wsum += w[i];
        }
        for (std::size_t i = 0; i< N; i++) { w[i] /= wsum; }
    } else {
        for (std::size_t i = 0; i< N; i++) { w[i] =1e0; }
    }

    // setup frequency grid
    if (df == -1) {
        double peak_width = 1e0 / (t[N-1] - t[0]);
        df = peak_width / freq_oversampling;
    }
    if (Nf == -1) {
        double avg_nyquist = (.5 * N) / (t[N-1] - t[0]);
        Nf = std::max(16, static_cast<int>((nyquist_factor*avg_nyquist-f0)/df));
    }
    assert( df > 0 && Nf > 0 );

    // Center the data
    if ( center_data || fit_offset ) {
        double dot = std::inner_product(w, w+N, y, 0e0);
        std::transform(y, y+N, y, [dot](double a){return a-dot;} );
    }
