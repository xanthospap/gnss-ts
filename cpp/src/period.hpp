#ifndef __NGPT_PERIOD_HPP__
#define __NGPT_PERIOD_HPP__

// c++ standard headers
#include <stdexcept>
#include <cmath>

// gtms headers
#include "timeseries.hpp"

namespace ngpt
{

template<class T, class F>
    void lomb_scargle_period(const timeseries<T,F>& ts, double ofac, double hifac,
    double px[], double py[], int np, int& nout, int& jmax, double& prob)
{
    std::size_t N = ts.data_pts() - ts.skipped_pts();
    
    double ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
           sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy;
    double arg,wtemp,*wi,*wpi,*wpr,*wr, *ts_epochs, *ts_vals;
    
    nout = 0.5*ofac*hifac*N;
    if (nout > np) {
        throw std::out_of_range{"lomb_scargle_period: [ERROR] output arrays too short in period"};
    }

    // Allocate memory
    try {
        wi  = new double[N];
        wpi = new double[N];
        wpr = new double[N];
        wr  = new double[N];
        ts_epochs = new double[N];
        ts_vals   = new double[N];
    } catch (std::bad_alloc&) {
        throw 1;
    }

    auto ts_start = ts.cbegin(),
         ts_stop  = ts.cend();
    std::size_t index = 0;
    ave = var = 0;
    double prev_ave = 0;
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

    xmin = ts_epochs[0];
    xmax = ts_epochs[N-1];
    xdif = xmax - xmin;
    xave = 0.5*(xmax+xmin);
    pymax= 0.0;
    pnow = 1.0/(xdif*ofac); // starting frequency
    // Initialize values for the trigonometric recurrences at each data point.
    for (std::size_t j = 0; j < N; j++) { 
        arg    = D2PI*((ts_epochs[j]-xave)*pnow);
        wpr[j] = -2.0*sin(0.5*arg)*sin(0.5*arg);
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
    delete[] wr;
    delete[] wpr;
    delete[] wpi;
    delete[] wi;
    delete[] ts_epochs;
    delete[] ts_vals;

    // The end
    return;
}

}// namespace ngpt

#endif
