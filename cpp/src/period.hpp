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
    double arg,wtemp,*wi,*wpi,*wpr,*wr;
    
    nout = 0.5*ofac*hifac*N;
    if (nout > np) {
        throw std::out_of_range{"lomb_scargle_period: [ERROR] output arrays too short in period"};
    }

    // Get mean and variance of the input data.
    // avevar(y, n, ave, var);
    //if (var == 0.0) {
    //    throw std::runtime_error{"zero variance in period"};
    //}

    // Allocate memory
    try {
        wi  = new double[N];
        wpi = new double[N];
        wpr = new double[N];
        wr  = new double[N];
    } catch (std::bad_alloc&) {
        throw 1;
    }

    xmin = ts.first_valid_epoch().as_mjd();
    xmax = ts.last_valid_epoch().as_mjd();

    xdif = xmax - xmin;
    xave = 0.5*(xmax+xmin);
    pymax= 0.0;
    pnow = 1.0/(xdif*ofac); // starting frequency
    // Initialize values for the trigonometric recurrences at each data point.
    for (int j = 0; j < n; j++) { 
        arg    = D2PI*((x[j]-xave)*pnow);
        wpr[j] = -2.0*sin(0.5*arg)*sin(0.5*arg);
        wpi[j] = std::sin(arg);
        wr[j]  = std::cos(arg);
        wi[j]  = wpi[j];
    }
}

}// namespace ngpt

#endif
