#include <iostream>
#include <fstream>
#include "procalgo.hpp"
#include "period.hpp"
#include "cts_read.hpp"
#include "psd.hpp"

// help message
void help()
{
    std::cout<<"\nProgram lomb-scargle";
    std::cout<<"\nPurpose Read in a (coordinate) time-series file, and estimate";
    std::cout<<"\n        a model.";
    std::cout<<"\nUsage   lomb_scargle -i <ts-file> [-l <igs-log-file> -q <noa-earthquake-catalogue> -e <event-file>]";
    return;
}

/// Given a file (with path), get the filename (after the last '/' character)
/// and return the first 4 chars as string. This function is used to get the
/// station name from a station coordinate file.
/// E.g. provided /home/user/some/path/ssss.xyz the function will return 'ssss'.
std::string
split_path(std::string s)
{
    auto pos = s.find_last_of('/');
    if (pos == std::string::npos ) {
        if ( s.size() < 4) return std::string("xxxx");
        return s.substr(0, 4);
    }
    return s.substr(pos+1, 4);
}

/// Minimum earthquake magnitude to be considered.
double MIN_ERTHQ_MAG = 4.0e0;
double C1 = -3.5e0;
double C2 = 2.0e0;

/// Parameters for harmonic analysis
double minfreq = 0e0;
double maxfreq = 0e0;
double dfreq   = 0e0;

int
main(int argc, char* argv[])
{
    if (argc < 2) {
        help();
        std::cout<<"\n";
        return 1;
    }
    
    // Various filenames
    std::string ctsf;
    const char *log_file   = nullptr,
               *event_file = nullptr,
               *erthq_file = nullptr;
    bool ctsf_found       = false;
    bool automatic_harmonic_analysis = false;
    // Parse command line options 
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-i")) {
            // input cts file
            assert( argc >= i+1 );
            ++i;
            ctsf = argv[i];
            ctsf_found = true;
        } else if (!strcmp(argv[i], "-l")) {
            // station log file
            assert( argc >= i+1 );
            ++i;
            log_file = argv[i];
        } else if (!strcmp(argv[i], "-e")) {
            // event file
            assert( argc >= i+1 );
            ++i;
            event_file = argv[i];
        } else if (!strcmp(argv[i], "-q")) {
            assert( argc >= i+1 );
            ++i;
            erthq_file = argv[i];
        } else if (!strcmp(argv[i], "-m")) {
            assert( argc >= i+1 );
            ++i;
            MIN_ERTHQ_MAG = std::atof(argv[i]);
        } else if (!strcmp(argv[i], "-c1")) {
            assert( argc >= i+1 );
            ++i;
            C1 = std::atof(argv[i]);
        } else if (!strcmp(argv[i], "-c2")) {
            assert( argc >= i+1 );
            ++i;
            C2 = std::atof(argv[i]);
        } else if (!strcmp(argv[i], "-a")) {
            automatic_harmonic_analysis = true;
        } else {
            std::cerr<<"\n[DEBUG] Fuck is that? The switch "<<argv[i]<<" is unrelevant.";
        }
    }
    if (!ctsf_found) {
        std::cerr<<"\n[ERROR] Cannot do anything without a time-series file.";
        return 1;
    }

    std::string sname = split_path(ctsf);
    std::cout<<"\nAnalysis report for station: "<<sname;
    //std::cout<<"\nAnalysis output written to file: "<< (sname + ".prd");
     
    // Read in the time-series from the cts file; time-resolution is 
    // milliseconds.
    ngpt::crdts<ngpt::milliseconds> ts = 
        ngpt::cts_read<ngpt::milliseconds>(ctsf, sname);
    auto mean_crd_xyz = ts.mean_coordinates();
    
    // Transform to topocentric rf (asumming input ts was geocentric cartesian).
    ts.cartesian2topocentric();
    
    // Print a short report
    std::cout<<"\nShort report on time-series:";
    std::cout<<"\n\tTime interval (span) from "
        << ngpt::strftime_ymd_hms(ts.first_epoch()) << " to "
        << ngpt::strftime_ymd_hms(ts.last_epoch());
    std::cout<<"\n\tNumber of epochs in time-series: "<<ts.size();

    // Apply any external jump information (event/log file).
    if (event_file) {
        std::cout<<"\nApplying event-list file: \'"<<event_file<<"\'.";
        ts.apply_event_list_file(event_file);
    }
    if (log_file) {
        std::cout<<"\nApplying (igs) log file: \'"<<log_file<<"\'.";
        ts.apply_stalog_file(log_file);
    }
    
    // If earthquake catalogue file provided, read and apply the interesting
    // earthquakes.
    if (erthq_file) {
        std::cout<<"\nApplying earthquake catalogue file: \'"<<erthq_file<<"\'.";
        ngpt::earthquake_catalogue<ngpt::milliseconds> eq_cat {erthq_file};
        ts.apply_earthquake_catalogue(eq_cat, MIN_ERTHQ_MAG, C1, C2);
    }

    // Filter the event list; two events must be at least a week apart
    ngpt::datetime_interval<ngpt::milliseconds> aweek
        {ngpt::modified_julian_day{7}, ngpt::milliseconds{0}};
    auto initial_events = ts.events().size();
    ts.events().filter_earthquake_sequences(aweek);
    std::cout<<"\nEarthquake sequences removed; Number of events: "
        <<ts.events().size()
        <<" (removed "<<initial_events-ts.events().size()<<" earthquakes)";
    /*
    int min_interval_for_events = 10;
    std::cout<<"\nFiltering event list; min interval between two events is: "
        << min_interval_for_events<<" days.";
    ts.clear_event_list(ngpt::datetime_interval<ngpt::milliseconds>
        {ngpt::modified_julian_day{1}, ngpt::milliseconds{0}}, min_interval_for_events);
    std::cout<<"\nInitial List of Events:\n";
    ts.events().dump_event_list(std::cout);
    */

    // Must remove linear trend and offsets before searching for harmonic
    // signals. Hence try an initial estimate. Make a seperate model for each
    // component. Also print the time-series to a file named 'original.ts'
    // (outliers are marked!).
    // Note that the residuals of this estimation are stored at res_ts, which
    // is a new time-series.
    std::cout<<"\nPerforming initial modeling (no harmonics modeled but all events considered).";
    ngpt::ts_model<ngpt::milliseconds> xmodel { ts.events() };
    xmodel.mean_epoch() = ts.mean_epoch();
    auto ymodel {xmodel},
         zmodel {xmodel};
    auto res_ts = ts.qr_fit(xmodel, ymodel, zmodel);
    std::cout<<"\nPrinting original time-series (with marked outliers) to file: "
        << "\'original.ts\'";
    std::ofstream f1 {"original.ts"};
    ts.dump(f1, true, false);
    f1.close();

    // Treat each component individualy (for harmonic analysis). Make a vector
    // of components so we can easily iterate through.
    std::vector<ngpt::timeseries<ngpt::milliseconds, ngpt::pt_marker>*>
        components;
    components.push_back(&res_ts.x_component());
    components.push_back(&res_ts.y_component());
    components.push_back(&res_ts.z_component());
    std::vector<std::string> cmp_names = 
        {std::string("North"), std::string("East"), std::string("Up")};

    if (automatic_harmonic_analysis) {
        double Ut = 1e-2;
        double min_ampl = 1e-3;
        xmodel = ngpt::identify_harmonics(res_ts.x_component(), xmodel, Ut, min_ampl);
        ymodel = ngpt::identify_harmonics(res_ts.y_component(), ymodel, Ut, min_ampl);
        zmodel = ngpt::identify_harmonics(res_ts.z_component(), zmodel, Ut, min_ampl);
    } else {
        // Time-span of the time-series in days and years.
        auto   tdif = res_ts.last_epoch().delta_date(res_ts.first_epoch());
        double ddif = tdif.days().as_underlying_type() / 365.25;
        double div  = 1e0/ddif;

        auto it  = components.begin();
        auto cit = cmp_names.cbegin();
        char answr;
        // Iterate through the components and perform harmonic analysis (via the
        // Lomb-Scargle periodogram).
        for (; it != components.end(); ++it) {
            ngpt::timeseries<ngpt::milliseconds, ngpt::pt_marker> tts {**it};
            answr = 'y';
            std::cout<<"\nHarmonic Analysis of Component: "<<(*cit);
            std::cout<<"\n-------------------------------------------------------";
            std::cout<<"\nComponent written to file: "<<std::string(sname + *cit + std::string(".cmp"));
            std::ofstream fouc { sname + *cit + std::string(".cmp") };
            tts.dump(fouc);
            fouc.close();
            while (answr != 'n' && answr != 'N') {
                std::size_t N = tts.data_pts() - tts.skipped_pts();
                double          ofac{4},
                                hifac{div/div},
                                *px,
                                *py,
                                prob,
                                *mempool;
                int             nout (0.5*ofac*hifac*N+1),
                                jmax;
                double          days_in_year = 365.25e0;
                if (dfreq) {
                    std::cout<<"\n[DEBUG] Total number of days: "<<ddif * 365.25e0;
                    minfreq = 1e0 / (ddif * 365.25e0); // 
                    maxfreq = 1e0 / 0.5e0;             // 0.5 days frequency
                    dfreq   = 1e-3;
                    nout    = static_cast<int>((maxfreq-minfreq)/dfreq )+1;
                }

                mempool = new double[2*nout];
                px      = mempool;
                py      = mempool + nout;
                if (!dfreq)
                    ngpt::lomb_scargle_period(tts, ofac, hifac, px, py, nout, nout, jmax, prob);
                else
                    ngpt::lomb_scargle_period(tts, minfreq, maxfreq, dfreq, px, py, nout, nout, jmax, prob);
                std::cout<<"\n\tDominant frequency in time-series: "<<px[jmax]<<" (at: "<<jmax<<")"
                    <<"; this is a period of "<<1e0/px[jmax]<<" days";
                std::cout<<"\n\tMinimum frequency examined is: "<<px[0]
                    <<", i.e. a period of "<<1e0/px[0]<<" days";
                std::cout<<"\n\tMaximum frequency examined is: "<<px[nout-1]
                    <<", i.e. a period of "<<1e0/px[nout-1]<<" days\n";
                std::cout<<"\nDo you want to apply the frequency to the model (y/n)?";
                std::cin>>answr;
                if (answr == 'y' || answr == 'Y') {
                    ngpt::ts_model<ngpt::milliseconds> *tmp_model;
                    if (*cit == "North")
                        tmp_model = &xmodel;
                    else if (*cit == "East")
                        tmp_model = &ymodel;
                    else
                        tmp_model = &zmodel;
                    tmp_model->add_period(1e0/px[jmax]);
                    std::cout<<"\nAdded period "<<1e0/px[jmax]<<" to the model.";
                    double post_stddev;
                    tts = tts.qr_ls_solve(*tmp_model, post_stddev, 1e-3, false, true); 
                }
                delete[] mempool;
            }
            ++cit;
        }
    }

    // Re-apply the models, now containing (maybe) harmonic signals.
    ts.qr_fit(xmodel, ymodel, zmodel);

    auto mdl_n = ngpt::filter_earthquakes(ts.x_component(), xmodel, 1e-3);
    auto mdl_e = ngpt::filter_earthquakes(ts.y_component(), ymodel, 1e-3);
    auto mdl_u = ngpt::filter_earthquakes(ts.z_component(), zmodel, 1e-3);

    // dump models to file
    std::cout<<"\nDumping component models to file: \'" 
        << std::string(sname + std::string(".mod")) <<"\'.";
    std::ofstream fout { sname + std::string(".mod") };
    fout <<"Approximate Station Coordinates: ";
    double lat, lon, hgt;
    ngpt::car2ell<ngpt::ellipsoid::grs80>(std::get<0>(mean_crd_xyz), 
        std::get<1>(mean_crd_xyz), std::get<2>(mean_crd_xyz), lat, lon, hgt);
    fout << ngpt::rad2deg(lat) << " " << ngpt::rad2deg(lon) << " " << hgt << "\n";
    mdl_n.dump(fout);
    fout<<"\n";
    mdl_e.dump(fout);
    fout<<"\n";
    mdl_u.dump(fout);
    fout.close();

    std::cout<<"\nModeling done; exiting!\n";
    return 0;
}
