#include <iostream>
#include <chrono>
#include "timeseries.hpp"
#include "genflags.hpp"

using ngpt::datetime;
using ngpt::year;
using ngpt::month;
using ngpt::day_of_month;
using ngpt::modified_julian_day;
using ngpt::milliseconds;

using ngpt::timeseries;
using ngpt::pt_marker;
using ngpt::flag;

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using cnsec = std::chrono::nanoseconds;
typedef std::chrono::duration<double> dsec;

void
rand_populate(timeseries<milliseconds, pt_marker>& ts, std::size_t N)
{
    for (std::size_t i = 0; i < N; i++) ts.add_point( (double)(i) );
}

int
main(/*int argc, char* argv[]*/)
{
    std::size_t idx, itmp;
    time_point<Clock> start, end;

    // Let's construct a vector of epochs, to hold 6000 epochs
    std::vector<datetime<milliseconds>> tvec;
    int NUM_EPOCHS = 6000;
    
    // start with date: 2010/01/01
    datetime<milliseconds> t0 {year{2010},
                               month{1},
                               day_of_month{1},
                               milliseconds{0}};

    // sequentialy add new time-stamps with a dt of half a day
    tvec.push_back(t0);
    milliseconds mls {12L * 60L * 60L * 1000L};
    modified_julian_day zero_days {0};
    for (int i = 1; i < NUM_EPOCHS; i++) {
        tvec.emplace_back( tvec[i-1].add(zero_days, mls) );
    }

    // now that we have a vector of epochs, let's make a time-series
    timeseries<milliseconds, pt_marker> ts {&tvec};

    // let's add the data-points (we only care about the values, not sigma/flags
    // at this point).
    for (int i = 0; i < NUM_EPOCHS; i++) {
        ts.add_point( (double)i );
    }

    // pint some basic information
    std::cout<<"\nTime-series info:";
    // print the first epoch and first valid epoch (as MJDs)
    std::cout<<"\n\tFirst epoch:         "<<ts.first_epoch().as_mjd();
    std::cout<<"\n\tFirst (valid) epoch: "<<ts.first_valid_epoch(idx).as_mjd();
    assert( idx == 0 );
    // print the last epoch and last valid epoch (as MJDs)
    std::cout<<"\n\tLast epoch:         "<<ts.last_epoch().as_mjd();
    std::cout<<"\n\tLast (valid) epoch: "<<ts.last_valid_epoch(idx).as_mjd();
    assert( (int)idx == NUM_EPOCHS-1 );
    // print the average
    std::cout<<"\n\tAverage = "<<ts.mean();
    // delta time in MJD between data points (two ways of getting an epoch at index i)
    auto dt = (ts.epoch_ptr()->operator[](1)).as_mjd() - /* get epoch at 1 */
              ts.epoch_at(0).as_mjd();                   /* get epoch at 0 */
    std::cout<<"\n\tDelta-time between epochs: "<<dt<<" days.";
    // the total time-span of the time-series in days
    auto tspan = ts.epoch_at(NUM_EPOCHS-1).as_mjd() - ts.epoch_at(0).as_mjd();
    assert( ngpt::delta_date(ts.epoch_at(NUM_EPOCHS), ts.epoch_at(0)) 
            == ts.cend().delta_time(ts.cbegin()) );
    std::cout<<"\n\tTotal time-span of time-series: "<<tspan<<" days.";
    
    // let's mark the first three data points and the last five, and re-print
    // the first and last epochs
    for (auto it = ts.begin(); it != ts.begin()+3; ++it) { 
        it.mark(pt_marker::outlier);
    }
    for (auto it = ts.end()-1; it != ts.end()-5; --it) {
        /* WARNING
         * do not use the following to mark a data-point, as this function does
         * not automaticaly update the number of skipped data-points
         * it.data().flag().set(pt_marker::outlier); 
         * The right way to go, is:
         */
        it.mark(pt_marker::outlier);
    }
    std::cout<<"\nAfter marking of observations:";
    std::cout<<"\n\tFirst epoch:         "<<ts.first_epoch().as_mjd();
    std::cout<<"\n\tFirst (valid) epoch: "<<ts.first_valid_epoch(idx).as_mjd();
    assert( idx == 3 );

    // print the last epoch and last valid epoch (as MJDs)
    std::cout<<"\n\tLast epoch:         "<<ts.last_epoch().as_mjd();
    std::cout<<"\n\tLast (valid) epoch: "<<ts.last_valid_epoch(idx).as_mjd();
    assert( (int)idx == NUM_EPOCHS-5 );

    // let's walk through all the time-series using const iterators
    idx = itmp = 0;
    for (auto cit = ts.cbegin(); cit != ts.cend(); ++cit) {
        if ( cit.data().skip() ) ++itmp;
        assert( cit.distance_from(ts.cbegin()) == (int)idx );
        idx++;
    }
    assert( (int)idx == NUM_EPOCHS );
    assert( itmp == ts.skipped_pts() );

    // the cost of using iterators Vs vector operators ...
    std::size_t iteration = 0;
    std::cout<<"\nMethod Time(ns) Size";
    for (int j = 0; j < 10; j++) {
        for (std::size_t num = 50+j; num < 1000000; num += num/2 ) {
            std::vector<datetime<milliseconds>> tmpvec (num, t0);
            timeseries<milliseconds, pt_marker> tts {&tmpvec};
            rand_populate(tts, num);
            if ((iteration++) % 2) {
                start = Clock::now();
                for (auto it = tts.begin(); it != tts.end(); ++it) {
                    if ( it.data().skip() ) ++itmp;
                    it.data().value() = 1e0;
                }
                end = Clock::now();
                auto tit = tts.end()-1;
                assert( tit.data().value() == 1e0 );
                cnsec diff1 = duration_cast<cnsec>(end - start);
                std::cout<<"\nIterators " << diff1.count() << " "<<num;
            } else {
                start = Clock::now();
                for (idx = 0; idx < (std::size_t)num; idx++) {
                    if ( tts[idx].skip() ) ++itmp;
                    tts[idx].value() = 2e0;
                }
                end = Clock::now();
                assert( tts[num-1].value() == 2e0 );
                cnsec diff2 = duration_cast<cnsec>(end - start);
                std::cout<<"\nVector acc. " << diff2.count() << " "<<num;
            }
            // auto ar  = std::chrono::duration_cast<dsec>(diff1-diff2) * 100e0;
            // auto par = std::chrono::duration_cast<dsec>(diff2);
            // std::cout<<"\nThat is a percentage difference of "<< ar/par <<"%";
            // std::cout<<"\n"<<num<<" "<<diff1.count()<<" "<<diff2.count()<<" "<<ar/par;
        }
    }

    std::cout<<"\n";
    return 0;
}
