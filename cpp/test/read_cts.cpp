#include <iostream>
#include <fstream>

#include "cts_read.hpp"
#include "crdts.hpp"
#include "genflags.hpp"

//
// This is a test program to check the ngpt::crdts<> class.
//

std::string
split_path(std::string s)
{
    auto pos = s.find_last_of('/');
    if (pos == std::string::npos ) {
        if ( s.size() < 4) return std::string("xxxx");
        pos = -1;
    }
    return s.substr(pos+1, 4);
}

int
main(int argc, char* argv[])
{
    if (argc < 2 || argc > 4) {
        std::cerr<<"Usage: read_cts <cts file> [events_file [catalogue]]\n";
        return 1;
    }

    std::string cts_file = std::string(argv[1]);
    std::string cts_name = split_path(std::string(argv[1]));

    // Read in the coordinate time-series from the input file
    std::cout << "Reading time-series file \"" <<cts_file<<"\", as station \""<<cts_name<<"\"\n";
    ngpt::crdts<ngpt::milliseconds> ts
        = ngpt::cts_read<ngpt::milliseconds>(cts_file, cts_name);

    // Transform cartesian to topocentric (including sigmas)
    ts.cartesian2topocentric();

    // If events file provided, read and apply the events
    if (argc > 2) {
        std::cout<<"\nApplying events list file: \""<<argv[2]<<"\".";
        ts.apply_event_list_file(argv[2]);
    }
    
    // If catalogue file provided, read and apply the interesting earthquakes
    if (argc > 3) {
        ngpt::earthquake_catalogue<ngpt::milliseconds> eq_cat
            {std::string(argv[3])};
        ts.apply_earthquake_catalogue(eq_cat);
    }

    // make a model and fit via ls (QR)
    std::vector<double> periods = { /*365.25*/ };
    ngpt::ts_model<ngpt::milliseconds> xmodel { ts.events() };
    xmodel.add_periods( periods );
    std::cout<<"\n--Manipulating dates";
    xmodel.mean_epoch() = ts.mean_epoch();
    std::cout<<"\nFrom: "<<ts.first_epoch().as_mjd()<<" to: "<<ts.last_epoch().as_mjd();
    std::cout<<"\nCentral Epoch: "<< ts.mean_epoch().as_mjd();
    std::cout<<"\nManipulating dates--";
    auto ymodel{xmodel}, zmodel{xmodel};
    std::cout<<"\n--Fitting";
    auto residual_ts = ts.qr_fit( xmodel, ymodel, zmodel );
    std::cout<<"\nFitting--";

    // test the iterator
    // ts.test_iter();
    // residual_ts.test_period();
    
    // test the running window
    ngpt::datetime_interval<ngpt::milliseconds> window {ngpt::modified_julian_day{30}, ngpt::milliseconds{0}};
    residual_ts.test_running_window(window);
    
    std::string filename = cts_name + std::string(".neu");
    // print the time-series
    std::ofstream fout_neu (filename);
    ts.dump( fout_neu );
    fout_neu.close();
    // print the residuals (time-series)
    filename = cts_name + std::string(".res");
    std::ofstream fout_res (filename);
    residual_ts.dump(fout_res);
    fout_res.close();
    // print the time-series event list
    filename = cts_name + std::string(".evn");
    std::ofstream fout_evn (filename);
    ts.dump_event_list( fout_evn );
    fout_evn.close();
    // print the time-series model line
    filename = cts_name + std::string(".mod");
    std::ofstream fout_mod (filename);
    ts.dump_model_line( fout_mod, xmodel, ymodel, zmodel );
    fout_mod.close();
    
    // --JSON--
    // print time-seres and residuals as json
    // filename = cts_name + std::string(".json");
    // std::ofstream fout_ts_json (filename);
    // ts.dump_json(fout_ts_json, residual_ts);
    // print the modelin json format
    std::ofstream fout_mod_json ("model.json");
    models_to_json(fout_mod_json, xmodel, ymodel, zmodel);
    std::ofstream fout_evn_json ("events.json");
    ts.dump_event_list_as_json(fout_evn_json);
    
    std::cout<<"\n";

    return 0;
}
