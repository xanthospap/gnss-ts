#include <iostream>
#include <fstream>
#include "period.hpp"
#include "cts_read.hpp"
#include "psd.hpp"

ngpt::ts_model<ngpt::milliseconds>
filter_earthquakes(ngpt::timeseries<ngpt::milliseconds, ngpt::pt_marker>& ts,
    ngpt::ts_model<ngpt::milliseconds>& model, double Ut);

void help()
{
    std::cout<<"\nProgram lomb-scargle";
    std::cout<<"\nPurpose Read in a (coordinate) time-series file, and estimate";
    std::cout<<"\n        the Lomb-Scargle periodogram.";
    std::cout<<"\nUsage   lomb_scargle ts-file";
    return;
}

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

double MIN_ERTHQ_MAG = 4.0e0;
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
    
    std::string ctsf;
    const char *log_file=nullptr,
               *event_file=nullptr,
               *erthq_file=nullptr;
    bool ctsf_found       = false,
         log_file_found   = false,
         event_file_found = false,
         erthq_file_found = false;
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
            log_file_found = true;
        } else if (!strcmp(argv[i], "-e")) {
            // event file
            assert( argc >= i+1 );
            ++i;
            event_file = argv[i];
            event_file_found = true;
        } else if (!strcmp(argv[i], "-q")) {
            // earthquake catalogue file
            assert( argc >= i+1 );
            ++i;
            erthq_file = argv[i];
            erthq_file_found = true;
        } else {
            std::cerr<<"\nFuck is that? The switch "<<argv[i]<<" is unrelevant.";
        }
    }
    if (!ctsf_found) {
        std::cerr<<"\nCannot do anything without a time-series file.";
        return 1;
    }

    std::string sname = split_path(ctsf);
    std::cout<<"\nAnalysis output written to file: "<< (sname + ".prd");
     
    // Read in the time-series from the cts file; time-resolution is 
    // milliseconds.
    ngpt::crdts<ngpt::milliseconds> ts = 
        ngpt::cts_read<ngpt::milliseconds>(ctsf, sname);
    
    // Transform to topocentric rf (asumming input ts was geocentric cartesian).
    ts.cartesian2topocentric();
    
    // Print a short report
    std::cout<<"\nShort report on time-series:";
    std::cout<<"\n\tTime interval (span) from "<<ngpt::strftime_ymd_hms(ts.first_epoch())<<" to "<<ngpt::strftime_ymd_hms(ts.last_epoch());
    std::cout<<"\n\tNumber of epochs in time-series: "<<ts.size();

    // Apply any external jump information (event/log file).
    if (event_file_found) {
        std::cout<<"\nApplying event-list file: "<<event_file<<".";
        ts.apply_event_list_file(event_file);
    }
    if (log_file_found) {
        std::cout<<"\nApplying (igs) log file: "<<log_file<<".";
        ts.apply_stalog_file(log_file);
    }
    
    // If earthquake catalogue file provided, read and apply the interesting
    // earthquakes
    if (erthq_file_found) {
        std::cout<<"\nApplying earthquake catalogue file: "<<erthq_file<<".";
        ngpt::earthquake_catalogue<ngpt::milliseconds> eq_cat {erthq_file};
        ts.apply_earthquake_catalogue(eq_cat, MIN_ERTHQ_MAG);
    }

    ts.clear_event_list(ngpt::datetime_interval<ngpt::milliseconds>{ngpt::modified_julian_day{1}, ngpt::milliseconds{0}}, 10);
    std::cout<<"\n\tEvents: ";
    ts.events().dump_event_list(std::cout);

    // Must remove linear trend before searching for harmonic signals. Include
    // any (external info on) jumps.
    ngpt::ts_model<ngpt::milliseconds> xmodel { ts.events() };
    xmodel.mean_epoch() = ts.mean_epoch();
    auto ymodel {xmodel},
         zmodel {xmodel};
    auto res_ts = ts.qr_fit(xmodel, ymodel, zmodel);
    std::ofstream f1 {"original.ts"};
    ts.dump(f1, true, false);
    f1.close();
    // ts = std::move(res_ts);
    std::ofstream f2 {"clear.ts"};
    res_ts.dump(f2, true, false);
    f2.close();

    // treat each component individualy (for harmonic analysis)
    std::vector<ngpt::timeseries<ngpt::milliseconds, ngpt::pt_marker>*> components;
    components.push_back(&res_ts.x_component());
    components.push_back(&res_ts.y_component());
    components.push_back(&res_ts.z_component());
    std::vector<std::string> cmp_names = 
        {std::string("North"), std::string("East"), std::string("Up")};

    auto   tdif = res_ts.last_epoch().delta_date(res_ts.first_epoch());
    double ddif = tdif.days().as_underlying_type() / 365.25;
    double div  = 1e0/ddif;

    auto it  = components.begin();
    auto cit = cmp_names.cbegin();
    char answr;
    for (; it != components.end(); ++it) {
        std::ofstream fouc { sname + *cit + std::string(".cmp") };
        (*it)->dump(fouc);
        fouc.close();
        std::size_t N = (*it)->data_pts() - (*it)->skipped_pts();
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
            ngpt::lomb_scargle_period(**it, ofac, hifac, px, py, nout, nout, jmax, prob);
        else
            ngpt::lomb_scargle_period(**it, minfreq, maxfreq, dfreq, px, py, nout, nout, jmax, prob);
        // fout << "\n#New Component : " << *cit;
        // for (int i=0;i<nout;i++) {
        //     fout << "\n" << px[i] << " " << py[i] << " " << 1e0/px[i];
        // }
        std::cout<<"\nComponent: "<<*cit;
        std::cout<<"\n---------------------------------------------------------";
        std::cout<<"\n\tDominant frequency in time-series: "<<px[jmax]<<" (at: "<<jmax<<")";
        std::cout<<"; this is a period of "<<1e0/px[jmax]<<" days";
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
            (*it)->qr_ls_solve(*tmp_model, post_stddev, 1e-3, false, true); 
        }

        ++cit;
        delete[] mempool;
        // fout.close();
    }

    // fit.
    ts.qr_fit(xmodel, ymodel, zmodel);

    /*
    auto vn = xmodel.make_model( *(ts.epoch_vector()) );
    auto ve = ymodel.make_model( *(ts.epoch_vector()) );
    auto vu = zmodel.make_model( *(ts.epoch_vector()) );
    std::ofstream fout { sname + std::string(".mod") };
    const std::vector<ngpt::datetime<ngpt::milliseconds>>* epoch_vec = ts.epoch_vector();
    for (std::size_t i=0; i<ts.size(); i++) {
        fout << ((*epoch_vec)[i]).as_mjd() << " " << vn[i] << " " << ve[i] << " " << vu[i] <<"\n";
    }
    fout.close();
    */

    auto mdl_n = filter_earthquakes(ts.x_component(), xmodel, 1e-3);
    auto mdl_e = filter_earthquakes(ts.y_component(), ymodel, 1e-3);
    auto mdl_u = filter_earthquakes(ts.z_component(), zmodel, 1e-3);

    std::ofstream fout { sname + std::string(".mod") };
    mdl_n.dump(fout);
    fout<<"\n";
    mdl_e.dump(fout);
    fout<<"\n";
    mdl_u.dump(fout);
    fout.close();

    return 0;
}

ngpt::ts_model<ngpt::milliseconds>
filter_earthquakes(ngpt::timeseries<ngpt::milliseconds, ngpt::pt_marker>& ts,
    ngpt::ts_model<ngpt::milliseconds>& model, double Ut=1e-3)
{
    
    std::vector<ngpt::md_earthquake<ngpt::milliseconds>> erthqk_vec {model.earthquakes()};
    if (!erthqk_vec.size()) return model;
    
    double stddev_a,
           stddev_n;
    
    /// Initial estimate to get the approximate earthquake offsets
    ts.qr_ls_solve(model, stddev_n, 1e-3, false, false);

    /// Get the a-priori earthquake vector and sort it according to each earthquake's
    /// resulting offset.
    //  -- already donde that!
    // std::vector<ngpt::md_earthquake<ngpt::milliseconds>> erthqk_vec {model.earthquakes()};
    std::sort(erthqk_vec.begin(), erthqk_vec.end(),
        [](const ngpt::md_earthquake<ngpt::milliseconds>& a,
           const ngpt::md_earthquake<ngpt::milliseconds>& b)
        {return std::sqrt(a.a1()*a.a1()+a.a2()*a.a2()) > std::sqrt(b.a1()*b.a1()+b.a2()*b.a2());});
    std::cout<<"\nInitial Earthquake offsets";
    for (auto it = erthqk_vec.begin(); it != erthqk_vec.end(); ++it) {
        std::cout<<"\n\tEarthquake at "<<ngpt::strftime_ymd_hms(it->start())<<" offset="<<it->a1();
    }

    /// Initialize a new model with no earthquakes
    ngpt::ts_model<ngpt::milliseconds> amodel{model},
                                       nmodel{model};
    nmodel.clear_earthquakes();
    amodel.clear_earthquakes();
    
    ts.qr_ls_solve(amodel, stddev_a, 1e-3, false, false);
    std::vector<ngpt::md_earthquake<ngpt::milliseconds>>::iterator it = erthqk_vec.begin(),
                                                             it_end = erthqk_vec.end();
    double factor;
    std::cout<<"\nStarting testing for significant offsets";
    for (; it != it_end; ++it) {
        nmodel.add_earthquake(*it);
        ts.qr_ls_solve(nmodel, stddev_n, 1e-3, false, false);
        factor = stddev_a / stddev_n;
        std::cout<<"\nPrevious std="<<stddev_a<<", new std="<<stddev_n;
        if ((factor-1e0) > Ut) {
            amodel = nmodel;
            stddev_a = stddev_n;
            std::cout<<"\n\tEarthquake at "<<ngpt::strftime_ymd_hms(it->start())<<" added to the model; factor is "<<factor<<", value is: "<<it->a1();
        } else {
            nmodel.erase_earthquake_at( it->start() );
            std::cout<<"\n\tEarthquake at "<<ngpt::strftime_ymd_hms(it->start())<<" removed from the model; factor is "<<factor<<", value is: "<<it->a1();
        }
    }

    return amodel;
}
