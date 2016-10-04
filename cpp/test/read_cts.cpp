#include <iostream>
#include <fstream>

#include "cts_read.hpp"
#include "crdts.hpp"
#include "genflags.hpp"

//
// This is a test program to check the ngpt::crdts<> class.
//

int
main(int argc, char* argv[])
{
    if (argc < 2 || argc > 4) {
        std::cerr<<"Usage: read_cts <cts file> [catalogue [events_file]]\n";
        return 1;
    }

    std::string cts_file = std::string(argv[1]);
    std::string cts_name = std::string("test");

    // Read in the coordinate time-series from the input file
    std::cout << "Reading time-series file \"" <<cts_file<<"\", as station \""<<cts_name<<"\"\n";
    ngpt::crdts<ngpt::milliseconds> ts
        = ngpt::cts_read<ngpt::milliseconds>(cts_file, cts_name);

    // Transform cartesian to topocentric (including sigmas)
    ts.cartesian2topocentric();
    
    // If catalogue file provided, read and apply the interesting earthquakes
    if (argc > 2) {
        ngpt::earthquake_catalogue<ngpt::milliseconds> eq_cat
            {std::string(argv[2])};
        ts.apply_earthquake_catalogue(eq_cat);
    }
    
    // If events file provided, read and apply the events
    if (argc > 3) {
        std::cout<<"\nApplying events list file: \""<<argv[3]<<"\".";
        ts.apply_event_list_file(argv[3]);
    }

    // fit model via ls (QR)
    std::vector<double> periods = { 365.25, 365.25/2 };
    ngpt::ts_model<ngpt::milliseconds> xmodel { ts.events() };
    xmodel.add_periods( periods );
    auto ymodel{xmodel}, zmodel{xmodel};
    auto residual_ts = ts.qr_fit( xmodel, ymodel, zmodel );

    // test the iterator
    // ts.test_iter();
    
    // test the running window
    // ngpt::datetime_interval<ngpt::milliseconds> window {ngpt::modified_julian_day{30}, ngpt::milliseconds{0}};
    // ts.test_running_window(window);
    
    // print the time-series
    std::ofstream fout_neu ("test.neu");
    ts.dump( fout_neu );
    fout_neu.close();
    // print the time-series event list
    std::ofstream fout_evn ("test.evn");
    ts.dump_event_list( fout_evn );
    fout_evn.close();
    // print the time-series model line
    std::ofstream fout_mod ("test.mod");
    ts.dump_model_line( fout_mod, xmodel, ymodel, zmodel );
    fout_mod.close();
    // print as json
    ts.dump_json( std::cout, residual_ts );

    /*
    ngpt::data_point<ngpt::ts_event> dp;
    std::cout<<"\nEmpty data point flag is: " << dp.flag();
    std::cout<<"\nData point is earthquake? "<<dp.flag().check(ngpt::ts_event::earthquake);
    std::cout<<"\nData point is vel.change? "<<dp.flag().check(ngpt::ts_event::velocity_change);
    std::cout<<"\nData point is outlier   ? "<<dp.flag().check(ngpt::ts_event::outlier);
    std::cout<<"\nData point is jump      ? "<<dp.flag().check(ngpt::ts_event::jump);
    std::cout<<"\nData point is skip      ? "<<dp.flag().check(ngpt::ts_event::skip);
    dp.flag().set(ngpt::ts_event::outlier);
    std::cout<<"\nAfter setting outlier flag is: "<<dp.flag();
    std::cout<<"\nData point is earthquake? "<<dp.flag().check(ngpt::ts_event::earthquake);
    std::cout<<"\nData point is vel.change? "<<dp.flag().check(ngpt::ts_event::velocity_change);
    std::cout<<"\nData point is outlier   ? "<<dp.flag().check(ngpt::ts_event::outlier);
    std::cout<<"\nData point is jump      ? "<<dp.flag().check(ngpt::ts_event::jump);
    std::cout<<"\nData point is skip      ? "<<dp.flag().check(ngpt::ts_event::skip);
    ngpt::flag<ngpt::ts_event> flg;
    std::cout<<"\nEmpty data point flag is: " << flg;
    std::cout<<"\nData point is earthquake? "<<flg.check(ngpt::ts_event::earthquake);
    std::cout<<"\nData point is vel.change? "<<flg.check(ngpt::ts_event::velocity_change);
    std::cout<<"\nData point is outlier   ? "<<flg.check(ngpt::ts_event::outlier);
    std::cout<<"\nData point is jump      ? "<<flg.check(ngpt::ts_event::jump);
    std::cout<<"\nData point is skip      ? "<<flg.check(ngpt::ts_event::skip);
    flg.set(ngpt::ts_event::outlier);
    std::cout<<"\nAfter setting outlier flag is: "<<flg;
    std::cout<<"\nData point is earthquake? "<<flg.check(ngpt::ts_event::earthquake);
    std::cout<<"\nData point is vel.change? "<<flg.check(ngpt::ts_event::velocity_change);
    std::cout<<"\nData point is outlier   ? "<<flg.check(ngpt::ts_event::outlier);
    std::cout<<"\nData point is jump      ? "<<flg.check(ngpt::ts_event::jump);
    std::cout<<"\nData point is skip      ? "<<flg.check(ngpt::ts_event::skip);
    flg.clear(ngpt::ts_event::outlier);
    std::cout<<"\nAfter clearing outlier flag is: "<<flg;
    std::cout<<"\nData point is earthquake? "<<flg.check(ngpt::ts_event::earthquake);
    std::cout<<"\nData point is vel.change? "<<flg.check(ngpt::ts_event::velocity_change);
    std::cout<<"\nData point is outlier   ? "<<flg.check(ngpt::ts_event::outlier);
    std::cout<<"\nData point is jump      ? "<<flg.check(ngpt::ts_event::jump);
    std::cout<<"\nData point is skip      ? "<<flg.check(ngpt::ts_event::skip);
    flg.set(ngpt::ts_event::skip);
    std::cout<<"\nAfter setting skip flag is: "<<flg;
    std::cout<<"\nData point is earthquake? "<<flg.check(ngpt::ts_event::earthquake);
    std::cout<<"\nData point is vel.change? "<<flg.check(ngpt::ts_event::velocity_change);
    std::cout<<"\nData point is outlier   ? "<<flg.check(ngpt::ts_event::outlier);
    std::cout<<"\nData point is jump      ? "<<flg.check(ngpt::ts_event::jump);
    std::cout<<"\nData point is skip      ? "<<flg.check(ngpt::ts_event::skip);
    flg.set(ngpt::ts_event::earthquake);
    std::cout<<"\nAfter setting earthquake flag is: "<<flg;
    std::cout<<"\nData point is earthquake? "<<flg.check(ngpt::ts_event::earthquake);
    std::cout<<"\nData point is vel.change? "<<flg.check(ngpt::ts_event::velocity_change);
    std::cout<<"\nData point is outlier   ? "<<flg.check(ngpt::ts_event::outlier);
    std::cout<<"\nData point is jump      ? "<<flg.check(ngpt::ts_event::jump);
    std::cout<<"\nData point is skip      ? "<<flg.check(ngpt::ts_event::skip);
    */

    std::cout<<"\n";

    return 0;
}
