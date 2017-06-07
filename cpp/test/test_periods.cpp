#include <iostream>
#include <fstream>
#include "artificial.hpp"
#ifdef DEBUG
#include <fenv.h>
#endif

using namespace ngpt;

int
main(int argc, char* argv[])
{
    int algorithm = 1;
    if (argc>1) algorithm = atoi(argv[1]);


    // Do not touch from here ----
    datetime<milliseconds> start 
        { modified_julian_day{53005}, milliseconds{0} };    // 01-01-2004
    datetime<milliseconds> stop
        { modified_julian_day{55562}, milliseconds{0} };    // 01-01-2011
    datetime_interval<milliseconds> step
        { modified_julian_day{1}, milliseconds{0} };        // 1-day step size
    long mean_mjd = (stop.as_mjd()-start.as_mjd()/2) + start.as_mjd();
    datetime<milliseconds> mean
        { modified_julian_day{mean_mjd}, milliseconds{0} }; // mean date

    // create a vector of epochs (from start to step every step)
    std::vector<datetime<milliseconds>> epochs;
    for (auto t = start; t < stop; t+=step) epochs.emplace_back(t);
    // -- to here

    // create a model
    ngpt::ts_model<milliseconds> model;
    model.mean_epoch() = mean;
    model.x0()         = 0e0;  // constant coef.
    model.vx()         = 0e0;  // velocity
    // add a jump at 1/1/2009
    // datetime<milliseconds> event {modified_julian_day{54832}, milliseconds{0}};
    // model.add_jump(event, 1.235); // add a jump
    // add harmonics
    model.add_period(365.25/52, -0.20, -0.51); // harmonic every week
    model.add_period(365.25/12, -0.60, -0.11); // harmonic every month
    std::cout<<"\n> Harmonics in signal:";
    for ( auto i : model.harmonics() ) {
        std::cout<<"\n\tPeriod: "<<i.period()<<" days"
            <<" angular frequency: "<<i.angular_frequency()
            <<" frequency: "<<i.angular_frequency()/D2PI;
    }

    //  make a synthetic time-series based on the epochs and model we already
    //+ constructed
    /* timeseries<milliseconds,pt_marker> */
    auto ts = synthetic_ts<milliseconds,pt_marker>(epochs, model, 0, 0.05);

    // write time-series to "foo.ts"
    std::cout<<"\n> TimeSeries written to \"foo.ts\"";
    std::ofstream ts_fout {"foo.ts"};
    ts.dump(ts_fout);

    if ( algorithm == 1 ) {
        // compute lomb-scargle periodogram; write to output ("lomb.out")
        std::cout<<"\n> Computing Lomb-Scargle Periodogram";
        std::size_t N = ts.data_pts() - ts.skipped_pts();
        double ofac{5}, hifac{1.5};
        int    nout = 0.5*ofac*hifac*N + 1;
        double *px, *py, prob;
        int    jmax;
        double days_in_year = 365.25e0;

        px = new double[nout];
        py = new double[nout];
        lomb_scargle_period(ts, ofac, hifac, px, py, nout, nout, jmax, prob);

        std::cout<<"\n\tDominant frequency in time-series: "<<px[jmax]<<" (at index: "<<jmax<<")"<<"; this is a period of "<<1e0/px[jmax]<<" days";
        std::cout<<"\n\tMinimum frequency examined is: "<<px[0]<<", i.e. a period of "<<1e0/px[0]<<" days.";
        std::cout<<"\n\tMaximum frequency examined is: "<<px[nout-1]<<", i.e. a period of "<<1e0/px[nout-1]<<" days.";
        std::cout<<"\n\tTesting frequency delta is "<<(px[2]-px[1]-px[0])<<" or every "<<1e0/(px[2]-px[1]-px[0])<<" days.";

        std::cout<<"\n> Writing Lomb-Scargle data to \"lomb.out\"";
        std::ofstream ls_out {"lomb.out"};
        for (int i=0;i<nout;i++) ls_out << "\n" << px[i] << " " << py[i];

        /*
        std::cout<<"\n> Iteratively searching and modeling dominant freqs.";
        // create a model
        ngpt::ts_model<milliseconds> mdl2;
        mdl2.mean_epoch() = mean;
        mdl2.x0()         = 0e0;  // constant coef.
        mdl2.vx()         = 0e0;  // velocity
        auto new_ts { ts };
        for (std::size_t i = 0; i <= model.harmonics().size(); i++) {
            new_ts.epoch_ptr() = ts.epoch_ptr();
            lomb_scargle_period(new_ts, ofac, hifac, px, py, nout, nout, jmax, prob);
            std::cout<<"\n\tIteration "<<i+1<<": Dominant frequency in time-series: "<<px[jmax]<<" (at index: "<<jmax<<")"
                <<"; this is a period of "<<1e0/px[jmax]<<" days power="<<py[jmax]<<" probability="<<prob;
            mdl2.add_period(1e0/px[jmax]);
            double std;
            new_ts = ts.qr_ls_solve(mdl2, std);
            new_ts.epoch_ptr() = ts.epoch_ptr();
            std::ofstream lot {"lomb" + std::to_string(i) + ".out"};
            for (int k=0;k<nout;k++) lot << "\n" << px[k] << " " << py[k];
            std::ofstream tot {"foo" + std::to_string(i) + ".ts"};
            new_ts.dump(tot);
        }
        */

        delete[] px;
        delete[] py;
    } else if ( algorithm == 2 ) {
        // compute lomb-scargle periodogram; write to output ("lomb.out")
        std::cout<<"\n> Computing Lomb-Scargle Periodogram";
        std::size_t idx;
        double T = ts.last_valid_epoch(idx).as_mjd() - ts.first_valid_epoch(idx).as_mjd();
        std::size_t N = ts.data_pts() - ts.skipped_pts();
        double minfreq /*{0}*/ ( 1e0/(1.5*T) );
        double maxfreq {.5};
        int    n0 {5};
        double dfreq {1e0/(T*n0)};
        int    nout = n0*T*maxfreq + 1;
        std::cout<<"\nComputed nout="<<nout<<", estimate="<<static_cast<int>((maxfreq-minfreq)/dfreq) + 1;
        // assert( nout == static_cast<int>((maxfreq-minfreq)/dfreq) + 1 );
        double *px, *py, prob;
        int    jmax;
        double days_in_year = 365.25e0;

        px = new double[nout];
        py = new double[nout];
        lomb_scargle_period(ts, minfreq, maxfreq, dfreq, px, py, nout, nout, jmax, prob);

        std::cout<<"\n\tDominant frequency in time-series: "<<px[jmax]<<" (at index: "<<jmax<<")"<<"; this is a period of "<<1e0/px[jmax]<<" days";
        std::cout<<"\n\tMinimum frequency examined is: "<<px[0]<<", i.e. a period of "<<1e0/px[0]<<" days.";
        std::cout<<"\n\tMaximum frequency examined is: "<<px[nout-1]<<", i.e. a period of "<<1e0/px[nout-1]<<" days.";
        std::cout<<"\n\tTesting frequency delta is "<<(px[2]-px[1]-px[0])<<" or every "<<1e0/(px[2]-px[1]-px[0])<<" days.";

        std::cout<<"\n> Writing Lomb-Scargle data to \"lomb.out\"";
        std::ofstream ls_out {"lomb.out"};
        for (int i=0;i<nout;i++) ls_out << "\n" << px[i] << " " << py[i];

        /*
        std::cout<<"\n> Iteratively searching and modeling dominant freqs.";
        // create a model
        ngpt::ts_model<milliseconds> mdl2;
        mdl2.mean_epoch() = mean;
        mdl2.x0()         = 0e0;  // constant coef.
        mdl2.vx()         = 0e0;  // velocity
        auto new_ts { ts };
        for (std::size_t i = 0; i <= model.harmonics().size(); i++) {
            new_ts.epoch_ptr() = ts.epoch_ptr();
            lomb_scargle_period(new_ts, minfreq, maxfreq, dfreq, px, py, nout, nout, jmax, prob);
            std::cout<<"\n\tIteration "<<i+1<<": Dominant frequency in time-series: "<<px[jmax]<<" (at index: "<<jmax<<")"
                <<"; this is a period of "<<1e0/px[jmax]<<" days power="<<py[jmax]<<" probability="<<prob;
            mdl2.add_period(1e0/px[jmax]);
            double std;
            new_ts = ts.qr_ls_solve(mdl2, std);
            new_ts.epoch_ptr() = ts.epoch_ptr();
            std::ofstream lot {"lomb" + std::to_string(i) + ".out"};
            for (int k=0;k<nout;k++) lot << "\n" << px[k] << " " << py[k];
            std::ofstream tot {"foo" + std::to_string(i) + ".ts"};
            new_ts.dump(tot);
        }
        */
        delete[] px;
        delete[] py;
    } else if (algorithm == 3) {
        // compute lomb-scargle periodogram; write to output ("lomb.out")
        std::cout<<"\n> Computing Lomb-Scargle Periodogram";
        std::size_t idx;
        double T = ts.last_valid_epoch(idx).as_mjd() - ts.first_valid_epoch(idx).as_mjd();
        std::size_t N = ts.data_pts() - ts.skipped_pts();
        double minfreq /*{0}*/ ( 1e0/(1.5*T) );
        double maxfreq {.5};
        double dfreq {1e0/(4*T)};
        int    nout = static_cast<int>((maxfreq-minfreq)/dfreq) + 1;
        std::cout<<"\nComputed nout="<<nout<<", estimate="<<static_cast<int>((maxfreq-minfreq)/dfreq) + 1;
        // assert( nout == static_cast<int>((maxfreq-minfreq)/dfreq) + 1 );
        double *px, *py, prob;
        int    jmax;

        px = new double[nout];
        py = new double[nout];
        lomb_scargle_period2(ts, minfreq, maxfreq, dfreq, px, py, nout, nout, jmax);

        std::cout<<"\n\tDominant frequency in time-series: "<<px[jmax]<<" (at index: "<<jmax<<")"<<"; this is a period of "<<1e0/px[jmax]<<" days";
        std::cout<<"\n\tMinimum frequency examined is: "<<px[0]<<", i.e. a period of "<<1e0/px[0]<<" days.";
        std::cout<<"\n\tMaximum frequency examined is: "<<px[nout-1]<<", i.e. a period of "<<1e0/px[nout-1]<<" days.";
        std::cout<<"\n\tTesting frequency delta is "<<(px[2]-px[1]-px[0])<<" or every "<<1e0/(px[2]-px[1]-px[0])<<" days.";

        std::cout<<"\n> Writing Lomb-Scargle data to \"lomb.out\"";
        std::ofstream ls_out {"lomb.out"};
        for (int i=0;i<nout;i++) ls_out << "\n" << px[i] << " " << py[i];
        delete[] px;
        delete[] py;
    } else {
        // compute lomb-scargle periodogram; write to output ("lomb.out")
        std::cout<<"\n> Computing Lomb-Scargle Periodogram (fast)";
        std::size_t N = ts.data_pts() - ts.skipped_pts();
        double ofac{5}, hifac{1.5};
        std::size_t nout = 0.5*ofac*hifac*N + 1;
        std::size_t nwk = lomb_scargle_fast_workspace(N, ofac, hifac);
        double *px, *py, prob;
        std::size_t jmax;
        double days_in_year = 365.25e0;

        px = new double[nwk];
        py = new double[nwk];
        lomb_scargle_fast(ts, ofac, hifac, px, py, nwk, nout, jmax, prob);

        std::cout<<"\n\tDominant frequency in time-series: "<<px[jmax]<<" (at index: "<<jmax<<")"<<"; this is a period of "<<1e0/px[jmax]<<" days";
        // std::cout<<"\n\tMinimum frequency examined is: "<<px[0]<<", i.e. a period of "<<1e0/px[0]<<" days.";
        // std::cout<<"\n\tMaximum frequency examined is: "<<px[nout-1]<<", i.e. a period of "<<1e0/px[nout-1]<<" days.";
        // std::cout<<"\n\tTesting frequency delta is "<<(px[2]-px[1]-px[0])<<" or every "<<1e0/(px[2]-px[1]-px[0])<<" days.";

        std::cout<<"\n> Writing Lomb-Scargle data to \"lomb.out\"";
        std::ofstream ls_out {"lomb.out"};
        for (std::size_t i=0;i<nout;i++) ls_out << "\n" << px[i] << " " << py[i];

        /*
        std::cout<<"\n> Iteratively searching and modeling dominant freqs.";
        // create a model
        ngpt::ts_model<milliseconds> mdl2;
        mdl2.mean_epoch() = mean;
        mdl2.x0()         = 0e0;  // constant coef.
        mdl2.vx()         = 0e0;  // velocity
        auto new_ts { ts };
        for (std::size_t i = 0; i <= model.harmonics().size(); i++) {
            new_ts.epoch_ptr() = ts.epoch_ptr();
            lomb_scargle_period(new_ts, ofac, hifac, px, py, nout, nout, jmax, prob);
            std::cout<<"\n\tIteration "<<i+1<<": Dominant frequency in time-series: "<<px[jmax]<<" (at index: "<<jmax<<")"
                <<"; this is a period of "<<1e0/px[jmax]<<" days power="<<py[jmax]<<" probability="<<prob;
            mdl2.add_period(1e0/px[jmax]);
            double std;
            new_ts = ts.qr_ls_solve(mdl2, std);
            new_ts.epoch_ptr() = ts.epoch_ptr();
            std::ofstream lot {"lomb" + std::to_string(i) + ".out"};
            for (int k=0;k<nout;k++) lot << "\n" << px[k] << " " << py[k];
            std::ofstream tot {"foo" + std::to_string(i) + ".ts"};
            new_ts.dump(tot);
        }
        */

        delete[] px;
        delete[] py;
    }

    std::cout<<"\n";
    return 0;
}
