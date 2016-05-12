#include <iostream>
#include <chrono>
// #include "datetime.hpp"
#include "datetime_v2.hpp"

using ngpt::datev2;

int main()
{
    // datetime<ngpt::datetime_clock::milli_seconds> d1 (2015, 12, 30);
    datev2<ngpt::seconds> d2(ngpt::year(2015), ngpt::month(12), ngpt::day_of_month(30));
    d2.add_seconds(ngpt::seconds(10));
    
    // nice, the following won't compile
    //datev2<ngpt::year> d3(ngpt::year(2015), ngpt::month(12), ngpt::day_of_month(30));
 
    // this is fine; microseconds to seconds is allowed
    datev2<ngpt::seconds> d24(ngpt::year(2015), ngpt::month(12),
        ngpt::day_of_month(30), ngpt::microseconds(100));

    // the opposite however id not allowed!
    //datev2<ngpt::microseconds> d5(ngpt::year(2015), ngpt::month(12),
    //    ngpt::day_of_month(30), ngpt::seconds(100));
            
    // std::cout<<"\nSize of v1 class: "<< sizeof(d1);
    std::cout<<"\nSize of v2 class: "<< sizeof(d2);

    std::chrono::steady_clock::time_point begin, end;
    double mjd1=0, mjd2;
/*
    begin = std::chrono::steady_clock::now();
    for (int i=0; i<86400 * 2.5; ++i) {
        d1.add_seconds( 1.0e0 );
        mjd1 = d1.as_mjd();
    }
*/
    end = std::chrono::steady_clock::now();

    //std::cout << "\nAdding two days in v1: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
    //std::cout << "\nmjd= " << mjd1;
    
    // if no user input, the noexcept may bias the computation!
    //long secs;
    //std::cout<<"\nSeconds to add: ";
    //std::cin >> secs;

    begin = std::chrono::steady_clock::now();
    for (int i=0; i<86400 * 2.5; ++i) {
        d2.add_seconds( ngpt::seconds(1L) );
        mjd2 = d2.as_mjd();
    }
    end = std::chrono::steady_clock::now();
    std::cout << "\nAdding two days in v2: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
    std::cout << "\nmjd= " << mjd2;
    std::cout<<"\nDifference (in days) = " << mjd1 - mjd2 << " = " << (mjd1 - mjd2)*86400.0 << " seconds";
    std::cout<<" = " << (mjd1 - mjd2)*86400000.0 << " milliseconds";
  
    datev2<ngpt::milliseconds> d3(ngpt::year(2015), ngpt::month(12), ngpt::day_of_month(30));
    begin = std::chrono::steady_clock::now();
    for (int i=0; i<86400 * 2.5; ++i) {
        d3.add_seconds( ngpt::milliseconds(1000L) );
        mjd2 = d3.as_mjd();
    }
    end = std::chrono::steady_clock::now();
    std::cout << "\nAdding two days in v2: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
    std::cout << "\nmjd= " << mjd2;
    std::cout<<"\nDifference (in days) = " << mjd1 - mjd2 << " = " << (mjd1 - mjd2)*86400.0 << " seconds";
    std::cout<<" = " << (mjd1 - mjd2)*86400000.0 << " milliseconds";
   
    datev2<ngpt::microseconds> d4(ngpt::year(2015), ngpt::month(12), ngpt::day_of_month(30));
    begin = std::chrono::steady_clock::now();
    for (int i=0; i<86400 * 2.5; ++i) {
        d4.add_seconds( ngpt::microseconds(1000000L) );
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
