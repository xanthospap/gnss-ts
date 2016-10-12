#include <iostream>
#include <fstream>
#include "artificial.hpp"

using namespace ngpt;

int
main(/*int argc, char* argv[]*/)
{
    datetime<milliseconds> start { modified_julian_day{53005}, milliseconds{0} };
    datetime<milliseconds> stop  { modified_julian_day{55562}, milliseconds{0} };
    datetime_interval<milliseconds> step { modified_julian_day{1}, milliseconds{0} };
    long mean_mjd = (stop.as_mjd()-start.as_mjd()/2)+start.as_mjd();
    datetime<milliseconds> mean { modified_julian_day{mean_mjd}, milliseconds{0} };

    // create a vector of epochs
    std::vector<datetime<milliseconds>> epochs;
    for (auto t = start; t < stop; t+=step) epochs.emplace_back(t);

    // create a model
    ngpt::ts_model<milliseconds> model;
    model.mean_epoch() = mean;
    model.x0() = 0e0;  // constant
    model.vx() = 0e0; // velocity
    // add a jump at 1/1/2009
    // datetime<milliseconds> event {modified_julian_day{54832}, milliseconds{0}};
    // model.add_jump(event, 1.235); // add a jump
    // add harmonics
    model.add_period(365.25/52, -0.40, -0.11);
    model.add_period(365.25/12, -0.80, -0.01);
    std::cout<<"\nThe model (and ts) include a frequency of 1/365.25/4 = "<<1/(365.25/2);

    /*timeseries<milliseconds,pt_marker>*/
    auto ts = synthetic_ts<milliseconds,pt_marker>(epochs, model, 0, 0.1);

    // write time-series
    std::ofstream ts_fout {"foo.ts"};
    ts.dump(ts_fout);

    // compute lom-scargle periodogram; write to output
    std::size_t N = ts.data_pts() - ts.skipped_pts();
    double ofac{4}, hifac{.9};
    int nout = 0.5*ofac*hifac*N + 1;
    double *px, *py, prob;
    int jmax;
    double days_in_year = 365.25e0;

    px = new double[nout];
    py = new double[nout];
    lomb_scargle_period(ts, ofac, hifac, px, py, nout, nout, jmax, prob);
    std::cout<<"\nDominant frequency in time-series: "<<py[jmax]<<" (at: "<<jmax<<")";
    std::cout<<"\nAs yearly fraction: 1/"<<days_in_year*px[jmax]<<"\n";
    std::ofstream ls_out {"lomb.out"};
    for (int i=0;i<nout;i++) ls_out << "\n" << px[i] << " " << py[i];
    
    delete[] px;
    delete[] py;

    return 0;
}
