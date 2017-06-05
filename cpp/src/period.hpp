#ifndef __NGPT_PERIOD_HPP__
#define __NGPT_PERIOD_HPP__

// c++ standard headers
#include <stdexcept>
#include <cmath>

// gtms headers
#include "timeseries.hpp"
#include "cfft.hpp"

/*
 * TRIGONOMETRIC RECURRENCE FORMULA
 * ================================
 *
 * Trig functions whose arguments form a linear sequence:
 * θ = θ(0) + nδ, n = 0, 1, 2, ...
 * are efficiently calculated by the following  recurrence:
 * cos(θ+δ) = cos(θ) - [αcos(θ)+βsin(θ)]
 * sin(θ+δ) = sin(θ) - [αsin(θ)+βcos(θ)]
 * where α and β are the precomputed coefficients:
 * α = 2sin^2(δ/2), β=sinδ
 * The reason for doing things this way, rather than with the standard (and 
 * equivalent) identities for sums of angles, is that here α and β do not lose
 * significance if the δ is small.
 */

namespace ngpt
{

/// Given an array yy[0..n-1], extirpolate (spread) a value y into m actual array
/// elements that best approximate the "fictional" (i.e. possibly non-integer)
/// array element number x. The weights used are coefficients of the Lagrange
/// interpolating polynomial.
/// This function is used for the fast Lomb-Scargle algorithm.
///
///  Reference: Numerical Recipes in C, Ch. 13.8
void
spread__(double y, double yy[], std::size_t n, double x, int m)
{
    long ihi,ilo,ix,j,nden;
    static long nfac[] = {0,1,1,2,6,24,120,720,5040,40320,362880};
    double fac;

    if (m > 10) {
        throw std::runtime_error{"spread__: factorial table too small in spread"};
    }

    ix = static_cast<long>(x);
    if (x == static_cast<double>(ix)) {
        yy[ix] += y;
    } else {
        ilo  = std::min( std::max(static_cast<long>(x-0.5*m+1.0), (long)1),
                        (long)(n-m+1));
        ihi  = ilo + m - 1;
#ifdef DEBUG
        assert( x   >= 0 );
        assert( ihi <  (long)n );
        assert( ilo >= 0 );
#endif
        nden = nfac[m];
        fac  = x-ilo;
        for (j = ilo+1;j <= ihi; j++) {
            fac *= (x-j);
        }
        yy[ihi] += y*fac/(nden*(x-ihi));
        for (j = ihi-1;j >= ilo; j--) {
            nden=(nden/(j+1-ilo))*(j-ihi);
            yy[j] += y*fac/(nden*(x-j));
        }
    }
    return;
}

std::size_t
lomb_scargle_fast_workspace(std::size_t N, double ofac, double hifac)
{
    constexpr int MACC {4};
    std::size_t nfreqt = ofac*hifac*N*MACC;
    std::size_t nfreq  = 64;
    while (nfreq < nfreqt) { nfreq <<= 1; }
    return nfreq << 1;
}

template<class T, class F>
void lomb_scargle_fast(const timeseries<T,F>& ts, double ofac, double hifac, 
    double wk1[], double wk2[], std::size_t nwk, std::size_t& nout, 
    std::size_t& jmax, double& prob)
{
    // Number of interpolation points per 1/4 cycle of highest frequency.
    constexpr int MACC {4};

    // real data size (i.e. ommiting outliers & skipped data points
    std::size_t  N = ts.data_pts() - ts.skipped_pts();
    
    std::size_t  j,k,ndim,nfreq,nfreqt;
    double       ave,ck,ckk,cterm,cwt,den,df,effm,expy,fac,fndim,hc2wt,
                 hs2wt,hypo,pmax,sterm,swt,var,xdif,xmax,xmin,*ts_epochs,*ts_vals;

    nout   = 0.5*ofac*hifac*N;
    // Size the FFT as next power of 2 above nfreqt
    nfreqt = ofac*hifac*N*MACC;
    nfreq  = 64;
    while (nfreq < nfreqt) { nfreq <<= 1; }
    ndim   = nfreq << 1;
    std::cout<<"\n\t(fasper) workspace = "<<ndim;
    if (ndim > nwk) {
        auto s_ndim = std::to_string(ndim);
        auto s_nwk  = std::to_string(nwk);
        throw std::runtime_error
        {"lomb_scargle_fast: workspaces too small ("+s_ndim+">"+s_nwk};
    }

    /// allocate memory
    double* MEM;
    try {
        MEM = new double[N*2];
    } catch (std::bad_alloc&) {
        throw 1;
    }
    ts_epochs = MEM;
    ts_vals   = MEM+N;

    auto ts_start = ts.cbegin(),
         ts_stop  = ts.cend();
    std::size_t index = 0;

    /// get mean and variance of the input data; also copy the ts data to
    /// the arrays ts_epochs & ts_vals (i.e. x and y data points). Only valid
    /// data points are considered.
    ave = var = 0;
    double prev_ave {0};
    for (auto it = ts_start; it != ts_stop; ++it) {
        if ( !it.data().skip() ) {
            ts_epochs[index] = it.epoch().as_mjd();
            ts_vals[index]   = it.data().value();
            prev_ave = ave;
            ave += (ts_vals[index]-ave)/(index+1);
            var += (ts_vals[index]-prev_ave)*(ts_vals[index]-ave);
            ++index;
        }
    }
    var /= (N-1);
    assert( index == N );

    xmin = ts_epochs[0];    // min epoch (MJD)
    xmax = ts_epochs[N-1];  // max epoch (MJD)
    xdif = xmax - xmin;     // difference in days
    
    //  Zero the workspaces
    for (j = 0; j< ndim; j++) { wk1[j] = wk2[j] = 0e0; }

    fac   = ndim/(xdif*ofac);
    fndim = ndim;
    // Extirpolate the data into the workspaces.
    for (j = 0; j < N; j++) {
        ck  = (ts_epochs[j]-xmin)*fac;
        while (ck > fndim) { ck -= fndim; }
        ckk = 2e0*(ck++);
        while (ckk > fndim) { ckk -= fndim; }
        ++ckk;
        spread__(ts_vals[j]-ave, wk1, ndim, ck, MACC);
        spread__(1e0, wk2, ndim, ckk, MACC);
    }

    // Take the Fourier Transforms
    realft__(wk1, ndim, 1);
    realft__(wk2, ndim, 1);
    df   = 1e0/(xdif*ofac);
    pmax = -1e0;

    // Compute the Lomb value for each frequency.
    for (k = 2, j = 0; j < nout; j++, k += 2) {
        hypo  = std::sqrt(wk2[k]*wk2[k]+wk2[k+1]*wk2[k+1]);
        hc2wt = 0.5e0*wk2[k]/hypo;
        hs2wt = 0.5e0*wk2[k+1]/hypo;
        cwt   = std::sqrt(0.5e0+hc2wt);
        swt   = std::copysign(std::sqrt(0.5e0-hc2wt),hs2wt);
        den   = 0.5e0*N+hc2wt*wk2[k]+hs2wt*wk2[k+1];
        cterm = std::sqrt(cwt*wk1[k]+swt*wk1[k+1])/den;
        sterm = std::sqrt(cwt*wk1[k+1]-swt*wk1[k])/(N-den);
        wk1[j]= j*df;
        wk2[j]= (cterm+sterm)/(2.0e0*var);
        if (wk2[j] > pmax) { pmax = wk2[(jmax=j)]; }
    }
    // Estimate significance of largest peak value.
    expy = exp(-pmax);
    effm = 2e0*(nout)/ofac;
    prob = effm*expy;
    if (prob > .01e0) { prob = 1e0-pow(1e0-expy,effm); }

    delete[] MEM;

}

///  Given n data points with abscissas x[0..n-1] (which need not be equally 
///  spaced) and ordinates y[0..n-1], and given a desired oversampling factor 
///  ofac (a typical value being 4 or larger), this routine fills array 
///  px[0..np-1] with an increasing sequence of frequencies (not angular 
///  frequencies) up to hifac times the “average” Nyquist frequency, and fills
///  array py[0..np-1] with the values of the Lomb normalized periodogram at 
///  those frequencies. The arrays x and y are not altered. np, the dimension
///  of px and py, must be large enough to contain the output, or an error 
///  results. The routine also returns jmax such that py[jmax] is the maximum 
///  element in py, and prob, an estimate of the significance of that maximum
///  against the hypothesis of random noise. A small value of prob indicates
///  that a significant periodic signal is present.
///
///  Reference: Numerical Recipes in C, Ch. 13.8
///
template<class T, class F>
    void lomb_scargle_period(const timeseries<T,F>& ts, double ofac, double hifac,
    double px[], double py[], int np, int& nout, int& jmax, double& prob)
{
    /// real data size (i.e. ommiting outliers & skipped data points
    std::size_t N = ts.data_pts() - ts.skipped_pts();
    
    double ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
           sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy;
    double arg,wtemp,*wi,*wpi,*wpr,*wr,*ts_epochs,*ts_vals;
    
    /// size of output arrays (# of frequencies to be examined)
    nout = 0.5 * ofac * hifac * N;
    if (nout > np) {
        throw std::out_of_range
            {"lomb_scargle_period: [ERROR] output arrays too short in period"};
    }

    /// allocate memory
    double* MEM;
    try {
        MEM = new double[N*6];
    } catch (std::bad_alloc&) {
        throw 1;
    }
    wi        = MEM;
    wpi       = MEM+N;
    wpr       = MEM+2*N;
    wr        = MEM+3*N;
    ts_epochs = MEM+4*N;
    ts_vals   = MEM+5*N;

    auto ts_start = ts.cbegin(),
         ts_stop  = ts.cend();
    std::size_t index = 0;

    /// get mean and variance of the input data; also copy the ts data to
    /// the arrays ts_epochs & ts_vals (i.e. x and y data points). Only valid
    /// data points are considered.
    ave = var = 0;
    double prev_ave {0};
    for (auto it = ts_start; it != ts_stop; ++it) {
        if ( !it.data().skip() ) {
            ts_epochs[index] = it.epoch().as_mjd();
            ts_vals[index]   = it.data().value();
            prev_ave = ave;
            ave += (ts_vals[index]-ave)/(index+1);
            var += (ts_vals[index]-prev_ave)*(ts_vals[index]-ave);
            ++index;
        }
    }
    var /= (N-1);
    assert( index == N );

    xmin = ts_epochs[0];    // min epoch (MJD)
    xmax = ts_epochs[N-1];  // max epoch (MJD)
    xdif = xmax - xmin;     // difference in days
    xave = 0.5*(xmax+xmin); // mean epoch
    pymax= 0.0;             // max frequency
    pnow = 1.0/(xdif*ofac); // starting frequency & frequency step

    // Initialize values for the trigonometric recurrences at each data point.
    for (std::size_t j = 0; j < N; j++) { 
        arg    = D2PI*((ts_epochs[j]-xave)*pnow);
        wpr[j] = -2.0*std::sin(0.5*arg)*std::sin(0.5*arg);
        wpi[j] = std::sin(arg);
        wr[j]  = std::cos(arg);
        wi[j]  = wpi[j];
    }
    // Main loop over the frequencies to be evaluated
    for (int i = 0; i < nout; i++) {
        px[i] = pnow;
        sumsh = sumc = 0.0;
        // First, loop over the data to get τ and related quantities.
        for (std::size_t j = 0; j < N; j++) {
            c = wr[j];
            s = wi[j];
            sumsh += s*c;
            sumc  += (c-s)*(c+s);
        }
        wtau  = 0.5*std::atan2(2.0*sumsh, sumc);
        swtau = std::sin(wtau);
        cwtau = std::cos(wtau);
        sums  = sumc = sumsy = sumcy = 0.0;
        // Then, loop over the data again to get the periodogram value.
        for (std::size_t j = 0; j < N; j++) {
            s      = wi[j];
            c      = wr[j];
            ss     = s*cwtau-c*swtau;
            cc     = c*cwtau+s*swtau;
            sums  += ss*ss;
            sumc  += cc*cc;
            yy     = ts_vals[j]-ave;
            sumsy += yy*ss;
            sumcy += yy*cc;
            // Update the trigonometric recurrences.
            wr[j] = ((wtemp=wr[j])*wpr[j]-wi[j]*wpi[j])+wr[j];
            wi[j] = (wi[j]*wpr[j]+wtemp*wpi[j])+wi[j];
        }
        py[i] = 0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var;
        if (py[i] >= pymax)
            pymax=py[(jmax=i)];
        pnow += 1.0/(ofac*xdif); // The next frequency.
    }

    // Evaluate statistical significance of the maximum.
    expy = std::exp(-pymax);
    effm = 2.0*nout/ofac;
    prob = effm*expy;
    if (prob > 0.01) prob=1.0-std::pow(1.0-expy,effm);

    // Free memory
    delete[] MEM;

    // The end
    return;
}

///
///  Given n data points with abscissas x[1..n] (which need not be equally 
///  spaced) and ordinates y[1..n], and given a desired oversampling factor 
///  ofac (a typical value being 4 or larger), this routine fills array 
///  px[1..np] with an increasing sequence of frequencies (not angular 
///  frequencies) up to hifac times the “average” Nyquist frequency, and fills
///  array py[1..np] with the values of the Lomb normalized periodogram at 
///  those frequencies. The arrays x and y are not altered. np, the dimension
///  of px and py, must be large enough to contain the output, or an error 
///  results. The routine also returns jmax such that py[jmax] is the maximum 
///  element in py, and prob, an estimate of the significance of that maximum
///  against the hypothesis of random noise. A small value of prob indicates
///  that a significant periodic signal is present.
///
///  Reference: Numerical Recipes in C, Ch. 13.8
///
template<class T, class F>
    void lomb_scargle_period(const timeseries<T,F>& ts, double minfreq,
    double maxfreq, double dfreq, double px[], double py[], int np, int& nout,
    int& jmax, double& prob)
{
    /// real data size (i.e. ommiting outliers & skipped data points
    std::size_t N = ts.data_pts() - ts.skipped_pts();
    
    double ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
           sumsy,swtau,var,wtau,xave,xmax,xmin,yy,wtemp;
    double arg,wctmp,wstmp,*wi,*wpi,*wpr,*wr,*ts_epochs,*ts_vals;
    
    /// size of output arrays (# of frequencies to be examined)
    nout = static_cast<int>( (maxfreq-minfreq)/dfreq )+1;
    if (nout > np) {
        throw std::out_of_range
            {"lomb_scargle_period: [ERROR] output arrays too short in period"};
    }

    /// allocate memory
    double* MEM;
    try {
        MEM = new double[N*6];
    } catch (std::bad_alloc&) {
        throw 1;
    }
    wi        = MEM;
    wpi       = MEM+N;
    wpr       = MEM+2*N;
    wr        = MEM+3*N;
    ts_epochs = MEM+4*N;
    ts_vals   = MEM+5*N;

    auto ts_start = ts.cbegin(),
         ts_stop  = ts.cend();
    std::size_t index = 0;

    /// get mean and variance of the input data; also copy the ts data to
    /// the arrays ts_epochs & ts_vals (i.e. x and y data points). Only valid
    /// data points are considered.
    ave = var = 0e0;
    double prev_ave {0e0};
    for (auto it = ts_start; it != ts_stop; ++it) {
        if ( !it.data().skip() ) {
            ts_epochs[index] = it.epoch().as_mjd();
            ts_vals[index]   = it.data().value();
            prev_ave         = ave;
            ave             += (ts_vals[index]-ave)/(index+1);
            var             += (ts_vals[index]-prev_ave)*(ts_vals[index]-ave);
            ++index;
        }
    }
    var = var/(N-1);
    assert( index == N );

    xmin  = ts_epochs[0];    // min epoch (MJD)
    xmax  = ts_epochs[N-1];  // max epoch (MJD)
    xave  = 0.5*(xmax+xmin); // mean epoch
    pnow  = minfreq;
    pymax = 0e0;

    // Initialize values for the trigonometric recurrences at each data point.
    /*
     * ω = 2πf => f = ω/2π
     * tμ  <- t_mean
     * arg <- 2π(t(i)-tμ)f = (ω2π(t(i)-tμ))/2π = ω(t(i)-tμ)
     * wpr <- -2 sin^2{ω(t(i)-tμ)/2} : for trig. recurrence (1)
     * wpi <- sin{ω(t(i)-tμ)}        : for trig. recurrence (1)
     * wr  <- cos{ω(t(i)-tμ)}
     * wi  <- sin{ω(t(i)-tμ)}
     *
     * (1) freq(n) = freq(0) + 2π(t(i)-tμ)*δf*n, n=0,1,2,....N-1
     */
    double rarg;
    for (std::size_t j = 0; j < N; j++) {
        arg    = D2PI*((ts_epochs[j]-xave)*pnow);
        rarg   = D2PI*((ts_epochs[j]-xave)*dfreq);
        wpr[j] = /*-2.0*std::sin(0.5*arg)*std::sin(0.5*arg);*/
                 -2.0*std::sin(0.5*rarg)*std::sin(0.5*rarg);
        wpi[j] = /*std::sin(arg);*/
                 std::sin(rarg);
        wr[j]  = std::cos(arg);
        wi[j]  = /*wpi[j];*/std::sin(arg);
    }
    // Main loop over the frequencies to be evaluated
    /*
     *  px[i] = f
     */
    for (int i = 0; i < nout; i++) {
        px[i] = pnow;
        sumsh = sumc = 0.0;
        // First, loop over the data to get τ and related quantities.
        /*
         * c     <- cos{ω(t(i)-tμ)}
         * s     <- sin{ω(t(i)-tμ)}
         * sumsh <- Σcos{ω(t(i)-tμ)}*sin{ω(t(i)-tμ)}
         * sumc  <- Σ(c-s)*(c+s)
         */
        for (std::size_t j = 0; j < N; j++) {
            c = wr[j];
            s = wi[j];
            sumsh += s*c;
            sumc  += (c-s)*(c+s);
        }
        /* 
         *  Lomb-Scargle: tan(2ωτ) = Σsin(2ωt_i) / Σcos(2ωt_j) [1]
         *  taking into account the trig. identities:
         *  sin(2θ) = 2 sinθ cosθ [2]
         *  cos(2θ) = cos^2(θ) - sin^2(θ) = {cos(θ)-sin(θ)}*{cos(θ)+sin(θ)} [3]
         *  [1]: ωτ = 0.5 * atan{Σsin(2ωt_i) / Σcos(2ωt_j)}
         *       ωτ = 0.5 * atan{Σ[2sin(ωt_i)*cos(ωt_i)] / 
         *                       Σ[{cos(ωt_i)-sin(ωt_i)}*{cos(ωt_i)+sin(ωt_i)}]}
         *       ωτ = 0.5 * atan(N*2*sumsh/sumc)
         */
        wtau  = 0.5*std::atan2(2.0*sumsh, sumc);
        swtau = std::sin(wtau);
        cwtau = std::cos(wtau);
        sums  = sumc = sumsy = sumcy = 0.0;
        // Then, loop over the data again to get the periodogram value.
        for (std::size_t j = 0; j < N; j++) {
            /*
             * Trig. identities:
             * sin(α-β) = sin(α)*cos(β) - cos(α)*sin(β)  [1]
             * cos(α-β) = cos(α)*cos(β) + sin(α)*sin(β)  [2]
             * Remember: ω(t(i)-tμ)
             *
             * s  = sin(arg)
             * c  = cos(arg)
             * ss = sin(arg)*cos(ωτ) - cos(arg)*sin(ωτ) -[1]-> ss = sin(arg-ωτ)
             * cc = cos(arg)*cos(ωτ) + sin(arg)*sin(ωτ) -[2]-> cc = cos(arg-ωτ)
             * where arg-ωτ = ω(t(i)-tμ)-ωτ = ω(t(i)-tμ-τ)
             *
             * sums = Σ(ss*ss)
             * sumc = Σ(cc*cc)
             * yy   = y_y - y_mean
             * sumsy = Σ(yy*ss)
             * sumcy = Σ(yy*cc)
             */
            s      = wi[j];
            c      = wr[j];
            ss     = s*cwtau-c*swtau;
            cc     = c*cwtau+s*swtau;
            sums  += ss*ss;
            sumc  += cc*cc;
            yy     = ts_vals[j]-ave;
            sumsy += yy*ss;
            sumcy += yy*cc;
            // Update the trigonometric recurrences.
            wtemp=wr[j];
            wr[j] = (wtemp*wpr[j]-wi[j]*wpi[j])+wr[j];
            wi[j] = (wi[j]*wpr[j]+wtemp*wpi[j])+wi[j];
        }
        py[i] = 0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var;
        if (py[i] >= pymax)
            pymax=py[(jmax=i)];
        pnow += dfreq; // The next frequency.
    }
#ifdef DEBUG
    std::cout<<"\n\tLast computed frequency = "<<pnow-dfreq<<"; max freq = "<<maxfreq;
#endif

    // Evaluate statistical significance of the maximum.
    expy = std::exp(-pymax);
    effm = 2.0*nout;
    prob = effm*expy;
    if (prob > 0.01) prob=1.0-std::pow(1.0-expy,effm);

    // Free memory
    delete[] MEM;

    // The end
    return;
}

// see https://arxiv.org/pdf/0901.2573.pdf
template<class T, class F>
    void lomb_scargle_period2(const timeseries<T,F>& ts, double minfreq,
    double maxfreq, double dfreq, double px[], double py[], int np, int& nout, int& jmax)
{
    /// real data size (i.e. ommiting outliers & skipped data points
    std::size_t N = ts.data_pts() - ts.skipped_pts();
    
    /// size of output arrays (# of frequencies to be examined)
    nout = static_cast<int>( (maxfreq-minfreq)/dfreq )+1;
    if (nout > np) {
        throw std::out_of_range
            {"lomb_scargle_period: [ERROR] output arrays too short in period"};
    }
    
    /// allocate memory
    double* MEM;
    try {
        MEM = new double[N*3];
    } catch (std::bad_alloc&) {
        throw 1;
    }
    double* ts_epochs = MEM;
    double* ts_vals   = MEM+N;
    double* ts_wght   = MEM+2*N;

    auto ts_start = ts.cbegin(),
         ts_stop  = ts.cend();
    std::size_t index = 0;

    /// get mean and variance of the input data; also copy the ts data to
    /// the arrays ts_epochs & ts_vals (i.e. x and y data points). Only valid
    /// data points are considered.
    double ave{0e0},prev_ave{0e0},W{0e0},Y{0e0},YY_hat{0e0},sigma;
    for (auto it = ts_start; it != ts_stop; ++it) {
        if ( !it.data().skip() ) {
            ts_epochs[index] = it.epoch().as_mjd();
            ts_vals[index]   = it.data().value();
            sigma            = it.data().sigma();
            ts_wght[index]   = 1e0/(sigma*sigma);
            W               += 1e0/(sigma*sigma);
            ave             += it.data().value()/sigma;
        }
    }
    ave /= N;
    double Wfac {W/N};

    /*
    for (std::size_t i = 0; i < N; i++) {
        Y               += ts_vals[i]/(sigma*sigma);
        YY_hat          += (it.data().value()*it.data().value())/(sigma*sigma);
        prev_ave         = ave;
        ave             += (ts_vals[index]-ave)/(index+1);
        ++index;
    }
    */
    //var = var/(N-1);
    //assert( index == N );
    double xmin = ts_epochs[0];    // min epoch (MJD)
    double xmax = ts_epochs[N-1];  // max epoch (MJD)
    double xave = 0.5*(xmax+xmin); // mean epoch

    double arg,carg,sarg,wnorm,C,S,YC_hat,YS_hat,CC_hat,SS_hat,CS_hat,
           YY,YC,YS,CC,SS,CS,D,max_freq{-100e0};
    std::size_t idx{0};
    jmax = 0;
    
    for (double omega = minfreq; omega < maxfreq; omega += dfreq) {
        C      = 0e0;
        S      = 0e0;
        YC_hat = 0e0;
        YS_hat = 0e0;
        CC_hat = 0e0;
        SS_hat = 0e0;
        CS_hat = 0e0;
        YY     = 0e0;
        YC     = 0e0;
        YS     = 0e0;
        for (std::size_t i = 0; i < N; i++) {
            arg     = omega * (ts_epochs[i]-xave);
            carg    = std::cos(arg);
            sarg    = std::sin(arg);
            wnorm   = ts_wght[i]/Wfac; // normalized weight
            C      += wnorm * carg;
            S      += wnorm * sarg;
            YC_hat += wnorm * ts_vals[i] * carg;
            YS_hat += wnorm * ts_vals[i] * sarg;
            CC_hat += wnorm * carg * carg;
            SS_hat += wnorm * sarg * sarg;
            CS_hat += wnorm * sarg * carg;
            YY     += wnorm * (ts_vals[i]-ave) * (ts_vals[i]-ave);
            YC     += wnorm * (ts_vals[i]-ave) * carg;
            YS     += wnorm * (ts_vals[i]-ave) * sarg;
        }
        //YY = YY_hat - Y*Y;
        //YC = YC_hat - Y*C;
        //YS = YS_hat - Y*S;
        CC = CC_hat - C*C;
        SS = SS_hat - S*S;
        CS = CS_hat - C*S;
        D  = CC*SS  - CS*CS;
        px[idx] = omega;
        py[idx] = (SS*YC*YC + CC*YS*YS - 2e0*CS*YC*YS)/(YY*D);
        if ( py[idx] > max_freq ) jmax = idx;
        ++idx;
    }

    delete[] MEM;

    return;
}


}// namespace ngpt

#endif
