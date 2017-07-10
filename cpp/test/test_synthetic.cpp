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
    if ( argc < 2 ) {
        std::cerr<<"\nError. Need to provide a file with earthquake coeffs.\n";
        return 1;
    }
    // open the file and read PSD coeffs
    std::ifstream fin (argv[1]);
    if (!fin.is_open()) {
        std::cerr<<"\nCould not open file \""+std::string(argv[1])+"\"";
        return 1;
    }
    double a1_ref, t1_ref, a2_ref, t2_ref;
    double a1_app, t1_app, a2_app, t2_app;
    fin >> a1_ref>> t1_ref>> a2_ref>> t2_ref;
    fin >> a1_app>> t1_app>> a2_app>> t2_app;
    fin.close();

    int ITERS = 10;
    psd_model psd_type = psd_model::log;
    if ( psd_type == psd_model::pwl ) ITERS = 1;

    // Do not touch from here ----
    datetime<milliseconds> start 
        { modified_julian_day{53005}, milliseconds{0} };    // 01-01-2004
    datetime<milliseconds> stop
        { modified_julian_day{55562}, milliseconds{0} };    // 01-01-2011
    datetime_interval<milliseconds> step
        { modified_julian_day{1}, milliseconds{0} };        // 1-day step size
    long mean_mjd = ((stop.as_mjd()-start.as_mjd())/2) + start.as_mjd();
    datetime<milliseconds> mean
        { modified_julian_day{mean_mjd}, milliseconds{0} }; // mean date

    // create a vector of epochs (from start to step every step)
    std::vector<datetime<milliseconds>> epochs;
    for (auto t = start; t < stop; t+=step) epochs.emplace_back(t);
    // -- to here

    // create a model
    ngpt::ts_model<milliseconds> ref_model;
    ref_model.mean_epoch() = mean;
    ref_model.x0()         = 0e0;  // constant coef.
    ref_model.vx()         = 5e-3;  // velocity
    // add a jump at 1/1/2009
    // datetime<milliseconds> event {modified_julian_day{54832}, milliseconds{0}};
    // ref_model.add_jump(event, 1.235); // add a jump
    // add harmonics
    // ref_model.add_period(365.25/52, -0.020, -0.0051); // harmonic every week
    // ref_model.add_period(365.25/12, -0.0060, -0.011); // harmonic every month
    // add an earthquake
    datetime<milliseconds> event
        {modified_julian_day{54693}, milliseconds{0}};  // 15-08-2008
    // ref_model.add_earthquake(event, -49.21e0/1e3, 1.9498e0, -39.81e0/1e3, 0.2373e0);
    ref_model.add_earthquake(event, psd_type, a1_ref, t1_ref, a2_ref, t2_ref);
    std::cout<<"\nReference Model";
    std::cout<<"\n------------------------------------------------------------\n";
    ref_model.dump(std::cout);

    //  make a synthetic time-series based on the epochs and ref_model we already
    //+ constructed
    /* timeseries<milliseconds,pt_marker> */
    auto ts = synthetic_ts<milliseconds,pt_marker>(epochs, ref_model, 0, 1e-3);
    // --just to see that this is working --
    /*auto*/timeseries<milliseconds,pt_marker> ts2 {ts,100,150};
    assert( ts2[0]  == ts[100] );
    assert( ts2[49] == ts[149] );

    // let's dare an estimate
    ngpt::ts_model<milliseconds> estim_mdl;
    estim_mdl.mean_epoch() = ref_model.mean_epoch();
    // estim_mdl.add_period(365.25/52);
    // estim_mdl.add_period(365.25/12);
    //  instead of adding an earthquake, lets add a jump (at the time of the
    //+ earthquake)
    estim_mdl.add_earthquake(event, psd_type, a1_app, t1_app, a2_app, t2_app);
    // estim_mdl.add_jump(event);
    // estim_mdl.add_velocity_change(event);
    // estim_mdl.add_jump(event);
    std::cout<<"\n\nApproximate Model";
    std::cout<<"\n------------------------------------------------------------\n";
    estim_mdl.dump(std::cout);
    double post_std_dev;
    ITERS  = 10;
    for (int i = 0; i < ITERS; i++) {
        ts.qr_ls_solve(estim_mdl, post_std_dev);
        std::cout<<"\n\nEstimated Model, iteration: "<<i;
        std::cout<<"\n------------------------------------------------------------\n";
        estim_mdl.dump(std::cout);
    }

    // This needs fixing; it doesn't work!
    // std::size_t idx;
    // ts.upper_bound(event, idx);
    // std::cout<<"\nIdx="<<idx<<"/"<<ts.epochs();
    // timeseries<milliseconds,pt_marker> ts_earth {ts, idx};
    // ts_earth.qr_ls_solve(estim_mdl, post_std_dev);

    // write time-series to "foo.ts"
    std::cout<<"\n> TimeSeries written to \"foo.ts\"";
    std::ofstream ts_fout {"foo.ts"};
    ts.dump(ts_fout);

    // write the 'model-line'
    std::vector<ngpt::datetime<milliseconds>> mjds;
    auto yy = estim_mdl.make_model(start, stop, step, &mjds);
    std::cout<<"\n> Model written to \"bar.ts\"";
    std::ofstream md_out {"bar.ts"};
    for (std::size_t i = 0; i < yy.size(); i++)
        md_out << "\n" << mjds[i].as_mjd() << " " << yy[i];
    md_out.close();

    std::cout<<"\n";
    return 0;
}
