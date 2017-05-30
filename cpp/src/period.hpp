#ifndef __NGPT_PERIOD_HPP__
#define __NGPT_PERIOD_HPP__

// c++ standard headers
#include <stdexcept>
#include <cmath>

// gtms headers
#include "timeseries.hpp"

namespace ngpt
{

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
    void lomb_scargle_period(const timeseries<T,F>& ts, double ofac, double hifac,
    double px[], double py[], int np, int& nout, int& jmax, double& prob)
{
    /// real data size (i.e. ommiting outliers & skipped data points
    std::size_t N = ts.data_pts() - ts.skipped_pts();
    
    double ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
           sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy;
    double arg,wtemp,*wi,*wpi,*wpr,*wr, *ts_epochs, *ts_vals;
    
    /// size of output arrays (# of frequencies to be examined)
    nout = 0.5 * ofac * hifac * N;
#ifdef DEBUG
    std::cout<<"\nComputed nout="<<nout;
#endif
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
    wi        = MEM;      /*new double[N];*/
    wpi       = MEM+N;    /*new double[N];*/
    wpr       = MEM+2*N;  /*new double[N];*/
    wr        = MEM+3*N;  /*new double[N];*/
    ts_epochs = MEM+4*N;  /*new double[N];*/
    ts_vals   = MEM+5*N;  /*new double[N];*/

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
/*
#ifdef DEBUG
    std::cout<<"\nAverage="<<xave;
    for (std::size_t j = 0; j < N; j++) {
        std::cout<<"\n"<<wpr[j]<<" "<<wpi[j]<<" "<<wr[j]<<" "<<wi[j];
    }
#endif
*/
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

//
// THIS DOES NOT WORK
template<class T, class F>
    void lomb_scargle_period(const timeseries<T,F>& ts, double minfreq, double maxfreq, double dfreq,
    double px[], double py[], int np, int& nout, int& jmax, double& prob)
{
    /// real data size (i.e. ommiting outliers & skipped data points
    std::size_t N = ts.data_pts() - ts.skipped_pts();
    
    double ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
           sumsy,swtau,var,wtau,xave,xmax,xmin,yy;
    double arg,wtemp,*wi,*wpi,*wpr,*wr, *ts_epochs, *ts_vals;
    
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
    var = var/(N-1);
    assert( index == N );

    xmin = ts_epochs[0];    // min epoch (MJD)
    xmax = ts_epochs[N-1];  // max epoch (MJD)
    xave = 0.5*(xmax+xmin); // mean epoch
    pnow  = minfreq;
    pymax = 0e0;

    // Initialize values for the trigonometric recurrences at each data point.
    for (std::size_t j = 0; j < N; j++) { 
        arg    = D2PI*((ts_epochs[j]-xave)*pnow);
        wpr[j] = -2.0*std::sin(0.5*arg)*std::sin(0.5*arg);
        wpi[j] = std::sin(arg);
        wr[j]  = std::cos(arg);
        wi[j]  = wpi[j];
    }
/*
#ifdef DEBUG
    std::cout<<"\nAverage="<<xave;
    for (std::size_t j = 0; j < N; j++) {
        std::cout<<"\n"<<wpr[j]<<" "<<wpi[j]<<" "<<wr[j]<<" "<<wi[j];
    }
#endif
*/
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
    void lomb_scargle_period2(const timeseries<T,F>& ts, double minfreq, double maxfreq, double dfreq,
    double px[], double py[], int np, int& nout, int& jmax, double& prob)
{
    /// real data size (i.e. ommiting outliers & skipped data points
    std::size_t N = ts.data_pts() - ts.skipped_pts();
    
    double ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
           sumsy,swtau,var,wtau,xave,xmax,xmin,yy;
    double arg,wtemp,*wi,*wpi,*wpr,*wr, *ts_epochs, *ts_vals;
    
    /// size of output arrays (# of frequencies to be examined)
    nout = static_cast<int>( (maxfreq-minfreq)/dfreq )+1;
    if (nout > np) {
        throw std::out_of_range
            {"lomb_scargle_period: [ERROR] output arrays too short in period"};
    }

    /// allocate memory
    double* MEM;
    try {
        MEM = new double[N*7];
    } catch (std::bad_alloc&) {
        throw 1;
    }
    wi        = MEM;
    wpi       = MEM+N;
    wpr       = MEM+2*N;
    wr        = MEM+3*N;
    ts_epochs = MEM+4*N;
    ts_vals   = MEM+5*N;
    ts_wght   = MEM+6*N;

    auto ts_start = ts.cbegin(),
         ts_stop  = ts.cend();
    std::size_t index = 0;

    /// get mean and variance of the input data; also copy the ts data to
    /// the arrays ts_epochs & ts_vals (i.e. x and y data points). Only valid
    /// data points are considered.
    ave = var = 0;
    double prev_ave {0e0}, W{0e0}, Y{0e0}, YY_hat{0e0};
    for (auto it = ts_start; it != ts_stop; ++it) {
        if ( !it.data().skip() ) {
            ts_epochs[index] = it.epoch().as_mjd();
            ts_vals[index]   = it.data().value();
            ts_wght[index]   = 1e0/(it.data().sigma()*it.data().sigma());
            W               += 1e0/(it.data().sigma()*it.data().sigma());
            Y               += it.data().value()/(it.data().sigma()*it.data().sigma());
            YY_hat          += (it.data().value()*it.data().value())/(it.data().sigma() * it.data().sigma());
            prev_ave         = ave;
            ave             += (ts_vals[index]-ave)/(index+1);
            var             += (ts_vals[index]-prev_ave)*(ts_vals[index]-ave);
            ++index;
        }
    }
    var = var/(N-1);
    assert( index == N );

    double Wfac {W/N};
    double arg,carg,sarg,wnorm,C,S,YC_hat,YS_hat,CC_hat,SS_hat,CS_hat,
           YY,YC,YS,CC,SS,CS,D;
    std::size_t idx{0};
    for (double omega = low; omega < high; omega += step) {
        for (std::size_t i = 0; i < N; i++) {
            arg     = omega * ts_epochs[i];
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
        }
        YY = YY_hat - Y*Y;
        YC = YC_hat - Y*C;
        YS = YS_hat - Y*S;
        CC = CC_hat - C*C;
        SS = SS_hat - S*S;
        CS = CS_hat - C*S;
        D  = CC*SS  - CS*CS;
        py[idx++] = (SS*YC*YC + CC*YS*YS - 2e0*CS*YC*YS)/(YY*D);
    }
}


}// namespace ngpt

#endif
