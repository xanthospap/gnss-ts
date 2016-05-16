#include <iostream>
#include <chrono>
#include "dtcalendar.hpp"

using ngpt::datetime;
using ngpt::year;
using ngpt::month;
using ngpt::hours;
using ngpt::minutes;
using ngpt::day_of_month;
using ngpt::seconds;
using ngpt::milliseconds;
using ngpt::microseconds;

long MilliSec = 1000L;
long MicroSec = 1000000L;

int main()
{
    
    //
    // Seconds, MicroSec and Millisec
    // -----------------------------------------------------------------------
    //
    seconds sec1 {10};
    milliseconds mlsec1 {10};
    microseconds mcsec1 {10};
    seconds sec2 {mcsec1};        // casting microsec to sec is allowed
    seconds sec3 {mlsec1};        // casting millisec to sec is allowed
    milliseconds mlsec2 {mcsec1}; // casting microsec to millisec is allowed
    // However, itis not allowed to cast from lower to higher precsision, e.g.
    // milliseconds ml1 {sec1};   // ERROR!
    // microseconds mc1 {mlsec1}; // ERROR!

    //
    // Construction of datetime objects
    // -----------------------------------------------------------------------
    //
    datetime<seconds> d2(year(2015), month(12), day_of_month(30));
    std::cout << "d2  = " << d2.stringify() << " (" << d2.secs() << ")\n";
    
    // the following won't compile; the template parameter can be
    // seconds, milliseconds or microseconds.
    // datetime<year> d3(year(2015), month(12), day_of_month(30)); // ERROR
 
    // this is fine; microseconds to seconds is allowed (BUT fractional sec 
    // are ignored!)
    datetime<seconds> d21 (year(2015), month(12), day_of_month(30),
        milliseconds(MilliSec));
    std::cout << "d21 = " << d21.stringify() << " (" << d21.secs() << ")\n";
    // the opposite however id not allowed!
    // datetime<microseconds> d5(year(2015), month(12), day_of_month(30),
    //     seconds(100)); ERROR!
    // we can also use time (i.e. hours, minutes, etc..)
    datetime<seconds> d22 (year(2015), month(12), day_of_month(30), hours(12),
        minutes(50), seconds(30));
    std::cout << "d22 = " << d22.stringify() << " (" << d22.secs() << ")\n";
    // or
    datetime<seconds> d23 (year(2015), month(12), day_of_month(30), hours(12),
        minutes(50), microseconds(30000001));
    std::cout << "d23 = " << d23.stringify() << " (" << d23.secs() << ")\n";
    // or
    datetime<microseconds> d24 (year(2015), month(12), day_of_month(30), hours(12),
        minutes(50), microseconds(30000001));
    std::cout << "d24 = " << d24.stringify() << " (" << d24.secs() << ")\n";
    // but not (seconds cannot be cast to milliseconds) ERROR!
    // datetime<milliseconds> d25 (year(2015), month(12), day_of_month(30), hours(12),
    //     minutes(50), seconds(30000001));
    // std::cout << "d25 = " << d25.stringify() << " (" << d25.secs() << ")\n";

    // this is fine; use fractional seconds (which are skipped!)
    // datetime<seconds> d22 (year(2015), month(12), day_of_month(30), 30.001234);
    // std::cout << "d22 = " << d22.stringify() << " (" << d22.secs() << ")\n";
    // or, for bigger accuracy ..
    //datetime<microseconds> d23 (year(2015), month(12), day_of_month(30),
    //    30.000001);
    //std::cout << "d23 = " << d23.stringify() << " (" << d23.secs() << ")\n";
    // same as
    //datetime<microseconds> d24 (year(2015), month(12), day_of_month(30),
    //    microseconds(1));
    //std::cout << "d24 = " << d24.stringify() << " (" << d24.secs() << ")\n";
        
    // std::cout<<"\nSize of v1 class: "<< sizeof(d1);
    std::cout<<"\nSize of v2 class: "<< sizeof(d2);
    
    //
    // Manipulation of datetime objects
    // -----------------------------------------------------------------------
    //
    d2.add_seconds(seconds(10));
    std::cout << d2.stringify() << "\n";

    std::chrono::steady_clock::time_point begin, end;
    double mjd1=0, mjd2;

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 86400 * 2.5; ++i) {
        d2.add_seconds( seconds(1L) );
        mjd2 = d2.as_mjd();
    }
    end = std::chrono::steady_clock::now();
    std::cout << "\nAdding two days in v2: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
    std::cout << "\nmjd= " << mjd2;
    std::cout<<"\nDifference (in days) = " << mjd1 - mjd2 << " = " << (mjd1 - mjd2)*86400.0 << " seconds";
    std::cout<<" = " << (mjd1 - mjd2)*86400000.0 << " milliseconds";
  
    datetime<milliseconds> d3(year(2015), month(12), day_of_month(30));
    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 86400 * 2.5; ++i) {
        d3.add_seconds( milliseconds(1000L) );
        mjd2 = d3.as_mjd();
    }
    end = std::chrono::steady_clock::now();
    std::cout << "\nAdding two days in v2: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
    std::cout << "\nmjd= " << mjd2;
    std::cout<<"\nDifference (in days) = " << mjd1 - mjd2 << " = " << (mjd1 - mjd2)*86400.0 << " seconds";
    std::cout<<" = " << (mjd1 - mjd2)*86400000.0 << " milliseconds";
   
    datetime<microseconds> d4(year(2015), month(12), day_of_month(30));
    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 86400 * 2.5; ++i) {
        d4.add_seconds( microseconds(1000000L) );
        mjd2 = d4.as_mjd();
    }
    end = std::chrono::steady_clock::now();
    std::cout << "\nAdding two days in v2: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
    std::cout << "\nmjd= " << mjd2;
    std::cout<<"\nDifference (in days) = " << mjd1 - mjd2 << " = " << (mjd1 - mjd2)*86400.0 << " seconds";
    std::cout<<" = " << (mjd1 - mjd2)*86400000.0 << " milliseconds";

    std::cout << "\n";
    return 0;   
}
